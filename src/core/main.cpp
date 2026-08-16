#include "camera.h"
#include "bvh.h"
#include "../gpu/memory.h"
#include "../gpu/intersect.h"
#include "../gpu/generateRaysKernel.h"
#include "../gpu/soa.h"
#include "../gpu/testShading.h"
#include "../math/math_types.h"
#include "../math/transform.h"
#include "../math/nurbs.h"
#include "../math/bezier_patch.h"
#include <cstdio>
#include <cuda_runtime_api.h>
#include <driver_types.h>


//remove later
static float patchDelta(const bezierPatch& P){
      auto deh = [&](int idx){                       // dehomogenize one control point
          MathTypes::vector4f c = P.control_point[idx];
          float iw = 1.0f / c.v[3];
          return MathTypes::vector3f{ c.v[0]*iw, c.v[1]*iw, c.v[2]*iw };
      };
      int p = P.p, q = P.q, row = q + 1;
      MathTypes::vector3f c00 = deh(0);              // (u,v)=(0,0)
      MathTypes::vector3f c10 = deh(p*row + 0);      // (1,0)
      MathTypes::vector3f c01 = deh(q);              // (0,1)
      MathTypes::vector3f c11 = deh(p*row + q);      // (1,1)

      float d = 0.0f;
      for(int i = 0; i <= p; ++i)
        for(int j = 0; j <= q; ++j){
          float u = (float)i / p, v = (float)j / q;                  // this point's node
          MathTypes::vector3f B = {                                   // bilinear at (u,v)
              (1-u)*(1-v)*c00.v[0] + u*(1-v)*c10.v[0] + (1-u)*v*c01.v[0] + u*v*c11.v[0],
              (1-u)*(1-v)*c00.v[1] + u*(1-v)*c10.v[1] + (1-u)*v*c01.v[1] + u*v*c11.v[1],
              (1-u)*(1-v)*c00.v[2] + u*(1-v)*c10.v[2] + (1-u)*v*c01.v[2] + u*v*c11.v[2] };
          MathTypes::vector3f Pij = deh(i*row + j);
          float dx=Pij.v[0]-B.v[0], dy=Pij.v[1]-B.v[1], dz=Pij.v[2]-B.v[2];
          d = fmaxf(d, sqrtf(dx*dx + dy*dy + dz*dz));                 // corners give 0
        }
      return d;
}

int main(){
	MathTypes::matrix4x4f world_from_camera = Transform::lookAt(	MathTypes::point3f{3,3,0}, 
																	MathTypes::point3f{0,0,6},
																	MathTypes::vector3f{0,0.5,0.5});
	int res_x = 1920;
	int res_y = 1080;

	
	camera_config main_cam = CameraOps::init_camera(world_from_camera, res_x, res_y);

	int max_pixels = res_x * res_y;
	
	
	//main_cam = CameraOps::getOrthographicCamera(main_cam, 0, 1000);
	float d_fov = 60;
	main_cam = CameraOps::getPerspectiveCamera(main_cam, 60, 0.01f, 1000.0f);

	SoAQueue::ray_queue ray_queue = {	CudaDevAlloc<float>(max_pixels),
										CudaDevAlloc<float>(max_pixels),
										CudaDevAlloc<float>(max_pixels),
										CudaDevAlloc<float>(max_pixels),
										CudaDevAlloc<float>(max_pixels),
										CudaDevAlloc<float>(max_pixels),
										CudaDevAlloc<int>(max_pixels)};
	
	// ---------- NURBS torus (biquadratic × biquadratic, rational) ----------
	const int p = 2, q = 2;
	int n = 8, m = 8;              // 9 × 9 control points
	const float R  = 1.0f;              // ring (center-of-tube) radius
	const float r  = 0.4f;              // tube radius
	const float cz = 6.0f;              // translate in front of the +Z camera
	const float S  = 0.70710678f;       // sqrt(2)/2  (circle mid-point weight)

	// 9-point rational unit circle: control points (cx,cy) and weights (cw)
	const float cx[9] = { 1, 1, 0,-1,-1,-1, 0, 1, 1};
	const float cy[9] = { 0, 1, 1, 1, 0,-1,-1,-1, 0};
	const float cw[9] = { 1, S, 1, S, 1, S, 1, S, 1};

	// full-circle knot vectors (4 quarter-arcs), length n+p+2 = 12
	float* U = new float[12]{0,0,0, 0.25f,0.25f, 0.5f,0.5f, 0.75f,0.75f, 1,1,1};
	float* V = new float[12]{0,0,0, 0.25f,0.25f, 0.5f,0.5f, 0.75f,0.75f, 1,1,1};

	MathTypes::vector4f* control_points = new MathTypes::vector4f[9*9];
	for(int i = 0; i < 9; ++i){                 // u : revolution around Z
		for(int j = 0; j < 9; ++j){             // v : tube circle in (radial, z)
			float radial = R + r * cx[j];       // distance from the Z axis
			float zz     = r * cy[j];           // tube height
			float w      = cw[i] * cw[j];       // tensor-product weight
			float x = radial * cx[i];
			float y = radial * cy[i];
			float z = zz + cz;                  // Euclidean point, pushed to +Z
			control_points[i*9 + j] = { w*x, w*y, w*z, w };  // store weighted-homogeneous
		}
	}
	
	const int C = 4; 
	NurbsOps::curvatureRefinementU(U, p, n, m, C, control_points);
	NurbsOps::curvatureRefinementV(V, q, m, n, C, control_points);

	ParametricSurf::nurbsSurf nurbs_prim = {U, V, control_points, n, m, p, q};	



	int nbU, nbV; 
	
	
	bezierPatch* patches = BezierPatchOps::getPatch(nurbs_prim, nbU, nbV);
	int primitive_size = nbU*nbV;

	int psz = (patches[0].p + 1)*(patches[0].q +1);
	int block_size = primitive_size * psz;

	MathTypes::vector4f* gpu_control_points = CudaAllocManaged<MathTypes::vector4f>(block_size);

	BVHNode* nodes = CudaAllocManaged<BVHNode>(2*(nbU*nbV)-1);
	bezierPatch* gpu_Patches = CudaAllocManaged<bezierPatch>(nbU*nbV);


	BvhOps::buildBVH(nodes, 0, patches, primitive_size);
	
	MathTypes::vector4f* base_pointer = patches[0].heap_pointer;


	memcpy(gpu_control_points, base_pointer, block_size*sizeof(MathTypes::vector4f));

	for(int i=0; i<primitive_size; i++){
		gpu_Patches[i] = patches[i];
		gpu_Patches[i].control_point	= gpu_control_points + (patches[i].control_point - base_pointer);
		gpu_Patches[i].heap_pointer 	= patches[i].heap_pointer;
	}

	for(int i=0; i<2*primitive_size-1; i++)
		printf("%f %f %f %f %f %f\n", nodes[i].node_bound.p_min_x, nodes[i].node_bound.p_min_y, nodes[i].node_bound.p_min_z, nodes[i].node_bound.p_max_x, nodes[i].node_bound.p_max_y, nodes[i].node_bound.p_max_z);

	BezierPatchOps::freePatches(patches);

	calculateRaysWrapper(ray_queue, perspective, main_cam);
	
	surfaceInteraction* surfaceInteract = CudaAllocManaged<surfaceInteraction>(max_pixels);	

	float* delta = CudaAllocManaged<float>(nbU*nbV);

	for(int i=0; i<nbU*nbV; i++){
		delta[i] = patchDelta(gpu_Patches[i]);
	}

	
	float flat_threshold = 0;	
	calculateIntersectWrapper(ray_queue, nodes, gpu_Patches, surfaceInteract, res_x, res_y, flat_threshold, MathUtils::radians(d_fov), delta);

	CudaFree(nodes);

	SoAQueue::frame_buffer screen = {	CudaAllocManaged<float>(res_x*res_y),
										CudaAllocManaged<float>(res_x*res_y),
										CudaAllocManaged<float>(res_x*res_y)};

	collisionTestShadingWrapper(ray_queue, screen, res_x, res_y);
	// Collision test with red pixels	
	FILE* ppmFile = fopen("CollisionPpmTest.ppm", "w");
	// ppm file header
	fprintf(ppmFile,"P3 \n%d %d \n255\n", res_x, res_y);
	for(int i=0; i<res_x*res_y; i++){
		
		fprintf(ppmFile,"%d %d %d\n", (int)(255.999f * screen.r[i]), 
											(int)(255.999f * screen.g[i]), 
											(int)(255.999f * screen.b[i]));
	}
	fclose(ppmFile);
	

	normalsTestShadingWrapper(ray_queue, screen, surfaceInteract, res_x, res_y);
	// Surface normals test
	ppmFile = fopen("SurfaceNormPpmTest.ppm", "w");
	// ppm file header
	fprintf(ppmFile,"P3 \n%d %d \n255\n", res_x, res_y);
	for(int i=0; i<res_x*res_y; i++){
		
		fprintf(ppmFile,"%d %d %d\n",	(int)(255.999f * screen.r[i]), 
										(int)(255.999f * screen.g[i]), 
										(int)(255.999f * screen.b[i]));
	}
	fclose(ppmFile);

	CudaFree(gpu_control_points);
	CudaFree(gpu_Patches);
	CudaFree(delta);
	CudaFree(ray_queue.direction_x);	
	CudaFree(ray_queue.direction_y);
	CudaFree(ray_queue.direction_z);	
	CudaFree(ray_queue.origin_x);	
	CudaFree(ray_queue.origin_z);	
	CudaFree(ray_queue.origin_y);	
	CudaFree(ray_queue.idx_interact);
	CudaFree(screen.r);
	CudaFree(screen.g);
	CudaFree(screen.b);
	CudaFree(surfaceInteract);


}
