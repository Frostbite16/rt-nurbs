#include "../math/math_types.h"
#include "../math/math_utils.h"
#include "../core/core_defines.h"
#include "../core/bvh.h"
#include "../math/bounds_ops.h"
#include "../math/bezier_patch.h"
#include "soa.h"
#include "intersect.h"
#include "interact.h"
#include <cmath>

using namespace VectorOps;


struct minimalIntersctPayload{
	MathTypes::vector3f hit_p;
	MathTypes::vector3f normal;
};

CUDA_GLOB 
void calculateIntersect(	BVHNode* nodes, 
							SoAQueue::ray_queue ray_queue, 
							bezierPatch* patches,
							surfaceInteraction* interactions,
							int screen_res_x,
							int screen_res_y,
							const float flat_threshold,
							const float cam_k,
							const float* delta,
							unsigned long long* newton_count,
							unsigned long long* skip_count){
	
	int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
	int idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	if(idx_x >= screen_res_x || idx_y >= screen_res_y)
		return;
	
	int idx = idx_y * screen_res_x + idx_x;
	MathExtra::ray ray = MathExtra::ray{MathTypes::vector3f{ 	ray_queue.origin_x[idx],
																ray_queue.origin_y[idx],
																ray_queue.origin_z[idx]},
										MathTypes::vector3f{	ray_queue.direction_x[idx],
																ray_queue.direction_y[idx],
																ray_queue.direction_z[idx]}};
	
	bool do_newton = true;
	unsigned long long nc_local = 0, sc_local = 0;

	// For ray-bound intersection calculation
	// you do not need a struct here, but just in the unlikely case the compiler
	// send this to the global memory, it is best to have this data continuously
	// rather than spread out
	MathTypes::vector3f ray_d_rec = MathTypes::vector3f{1/ray.direction.v[0],
														1/ray.direction.v[1],
														1/ray.direction.v[2]};

	MathExtra::ray_plane ray_plane = MathUtils::getRayPlanes(ray);

	int nodes_to_visit[32];
	int to_visit_offset = 0;
	float t_max = INFINITY;

	//intersectionPayload t_payload = {INFINITY, -1};
	
	BVHNode cur_node = nodes[0];

	// it calculates the intersection with the root node before the loop starts
	// this way we can already discard quickly the thread ray that do not hit the root
	const MathExtra::bounds3f* bound = &cur_node.node_bound;
	float intersect_t = BoundsOps::intersectBound(	bound->p_min_x, bound->p_max_x,
													bound->p_min_y, bound->p_max_y,
													bound->p_min_z, bound->p_max_z,
													ray.origin.v[0], ray.origin.v[1],
													ray.origin.v[2], ray_d_rec.v[0],
													ray_d_rec.v[1],  ray_d_rec.v[2]);
	int initial_test = intersect_t < t_max;	

	
	newtonRhapsonPayload nr_payload;
	// This algorithm was taken in parts from the pbrt bvh traversal implementation but also 
	// from https://jacco.ompf2.com/2019/07/18/wavefront-path-tracing/https://jacco.ompf2.com/2019/07/18/wavefront-path-tracing/
	// the loop will be ignored if the ray misses the initial box
	while(initial_test){
		if(cur_node.prim_count>0){
			intersect_t = BoundsOps::intersectBound(	cur_node.node_bound.p_min_x, cur_node.node_bound.p_max_x,
														cur_node.node_bound.p_min_y, cur_node.node_bound.p_max_y,
														cur_node.node_bound.p_min_z, cur_node.node_bound.p_max_z,
														ray.origin.v[0], ray.origin.v[1],
														ray.origin.v[2], ray_d_rec.v[0],
														ray_d_rec.v[1],  ray_d_rec.v[2]);


			for(int i=0; i<cur_node.prim_count; i++){
				do_newton = true;
				bezierPatch patch = patches[cur_node.left_first + i];
				MathTypes::vector3f q00;
				MathTypes::vector3f q10;
				MathTypes::vector3f q11;
				MathTypes::vector3f q01;
			
				corners(patch, q00, q10, q11, q01);
				MathTypes::vector2f g = BezierPatchOps::biliniearRay(q00, q10, q11, q01, ray.origin, ray.direction);	
				
				if(delta[cur_node.left_first + i]*cam_k / intersect_t < flat_threshold ){
					MathTypes::vector3f bl_S, S_u, S_v;
					BezierPatchOps::evalPointNormal(patch, g.v[0], g.v[1], ray.direction, bl_S, S_u, S_v);				
					float resid = 	fabsf(dotProduct(ray_plane.n1, bl_S) + ray_plane.d1) +
									fabsf(dotProduct(ray_plane.n2, bl_S) + ray_plane.d2);
					
					if(resid*cam_k/intersect_t < flat_threshold){
						sc_local++;
						float t = dotProduct(subtract(bl_S, ray.origin),ray.direction);
						if(0<t && t<t_max){
							t_max = t;
							interactions[idx].t_hit = t;
							interactions[idx].uv = {g.v[0], g.v[1]};
							interactions[idx].dpdu = S_u;
							interactions[idx].dpdv = S_v;
							ray_queue.idx_interact[idx] = idx;
							do_newton = false;
						}	
					}
				}

				if(do_newton){
					nc_local++;
					nr_payload = BezierPatchOps::newtonRhapson(ray_plane, patch, g.v[0], g.v[1], 1e-5f,1e-4f, idx);
					if(nr_payload.hit){
						float t=(nr_payload.s.v[0]-ray.origin.v[0])*ray.direction.v[0]+
								(nr_payload.s.v[1]-ray.origin.v[1])*ray.direction.v[1]+
								(nr_payload.s.v[2]-ray.origin.v[2])*ray.direction.v[2];
						if(0<t && t < t_max){
							interactions[idx].dpdu = nr_payload.su;
							interactions[idx].dpdv = nr_payload.sv;
							interactions[idx].uv = {nr_payload.u, nr_payload.v};
							interactions[idx].t_hit = t;
							ray_queue.idx_interact[idx] = idx;
							t_max = t;
						}
					}
				}
								
			}

			if(to_visit_offset == 0) break;
			cur_node = nodes[nodes_to_visit[--to_visit_offset]];

		}
		else{
			MathExtra::bounds3f left_bound = nodes[cur_node.left_first].node_bound;
			MathExtra::bounds3f right_bound = nodes[cur_node.left_first+1].node_bound;
			float left_t = BoundsOps::intersectBound(	left_bound.p_min_x, left_bound.p_max_x,
														left_bound.p_min_y, left_bound.p_max_y,
														left_bound.p_min_z, left_bound.p_max_z,
														ray.origin.v[0], ray.origin.v[1],
														ray.origin.v[2], ray_d_rec.v[0],
														ray_d_rec.v[1],  ray_d_rec.v[2]);

			float right_t = BoundsOps::intersectBound(	right_bound.p_min_x, right_bound.p_max_x,
														right_bound.p_min_y, right_bound.p_max_y,
														right_bound.p_min_z, right_bound.p_max_z,
														ray.origin.v[0], ray.origin.v[1],
														ray.origin.v[2], ray_d_rec.v[0],
														ray_d_rec.v[1],  ray_d_rec.v[2]);
			
			// Intersect_t is the t value from the closest intersection between the child nodes 
			intersect_t = fminf(left_t, right_t);

			if(intersect_t >= t_max){
				if(to_visit_offset==0) break;
				cur_node = nodes[nodes_to_visit[--to_visit_offset]];
			}
			else{
				if(left_t <= intersect_t){
					if(right_t < t_max)
						nodes_to_visit[to_visit_offset++] = cur_node.left_first+1;
					cur_node = nodes[cur_node.left_first];
				}
				else{
					if(left_t < t_max)
						nodes_to_visit[to_visit_offset++] = cur_node.left_first;
					cur_node = nodes[cur_node.left_first+1];
				}
			}
		}
	}


	if(newton_count) atomicAdd(newton_count, nc_local);
	if(skip_count)   atomicAdd(skip_count,   sc_local);

	if(ray_queue.idx_interact[idx] != -1)
		interactions[idx].wo = VectorOps::multiply(-1, ray.direction);
}




CORE_H
void calculateIntersectWrapper(	SoAQueue::ray_queue& ray_queue, BVHNode* nodes, bezierPatch* patches,
								surfaceInteraction* interact, int screen_res_x, int screen_res_y,
								float flat_threshold, float v_fov, float* delta,
								unsigned long long* newton_count, unsigned long long* skip_count){

	int thread_x_quant = 16;
	int thread_y_quant = 16;
	dim3 thread_2D(thread_x_quant,thread_y_quant);
	dim3 block_2D(
		(screen_res_x + thread_x_quant - 1) / thread_x_quant,
		(screen_res_y + thread_y_quant - 1) / thread_y_quant);

	float cam_k = 0.5 * (float)screen_res_y / tan(v_fov/2.0f);

	printf("camk: %f\n", cam_k);

	calculateIntersect<<<block_2D,thread_2D>>>(	nodes, ray_queue, patches, interact , screen_res_x, screen_res_y,
												flat_threshold, cam_k, delta, newton_count, skip_count);
	cudaDeviceSynchronize();

}



