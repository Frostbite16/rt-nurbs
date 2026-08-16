#ifndef BEZIER_PATCH_H
#define BEZIER_PATCH_H
#include "math_types.h"
#include "vector_ops.h"
#include "math_utils.h"
#include "matrix_ops.h"
#include "nurbs.h"
#include <cmath>

struct bezierPatch{
	MathTypes::vector4f* control_point;
	float u0, u1, v0, v1;
	int p, q;

	MathTypes::vector4f* heap_pointer; // So we can free the entire array of control points;
};

struct beziePatchInteraction{


};

// minimal structure for the newton-rhapson iterator return 
struct newtonRhapsonPayload{
	bool hit;
	MathTypes::vector3f s;
	MathTypes::vector3f su;
	MathTypes::vector3f sv;
	float u,v;
};

CORE_D CUDA_FINL 
MathTypes::vector3f dehom(const MathTypes::vector4f& c){ 
	float iw=1.f/c.v[3]; return {c.v[0]*iw,c.v[1]*iw,c.v[2]*iw};
}

CORE_D CUDA_FINL
void corners(	const bezierPatch& p, MathTypes::vector3f& q00, 
				MathTypes::vector3f& q10, MathTypes::vector3f& q11,
				MathTypes::vector3f& q01)
{
    const MathTypes::vector4f* cp = p.control_point; int P=p.p, Q=p.q;
    q00=dehom(cp[0]); q10=dehom(cp[P*(Q+1)]); q11=dehom(cp[P*(Q+1)+Q]); q01=dehom(cp[Q]);
}

inline int segs(const float* U, int last, int p, float* out){
	int nb = 0;
	out[nb++] = U[p];
	for(int i=p+1; i<=last-p-1; ++i)
		if(U[i] != U[i-1]) out[nb++] = U[i];
	out[nb] = U[last-p];
	return nb;
}

// Algorithm A1.5 (The Nurbs Book)
CORE_HD inline
MathTypes::vector4f deCasteljau1(const MathTypes::vector4f* P, int n, float u){
	MathTypes::vector4f Q[NURBS_MAX_DEGREE+1]; 
	for(int i=0; i<=n; i++)
		Q[i] = P[i];
	for(int k=1; k<=n; k++)
		for(int i=0; i<=n-k; i++)
			Q[i] = VectorOps::sum(VectorOps::multiply((1.0f-u), Q[i]), VectorOps::multiply(u, Q[i+1]));
	return Q[0];	
} 




// Algorithm A1.5 (The Nurbs Book)
CORE_HD inline
MathTypes::vector4f deCasteljauDers(const MathTypes::vector4f* P, int n, float u){
	MathTypes::vector4f Q[NURBS_MAX_DEGREE+1]; 
	for(int i=0; i<n; i++)
		Q[i] = VectorOps::multiply(n,VectorOps::subtract(P[i+1],P[i]));
	for(int k=1; k<n; k++)
		for(int i=0; i<=n-k; i++)
			Q[i] = VectorOps::sum(VectorOps::multiply((1.0f-u), Q[i]), VectorOps::multiply(u, Q[i+1]));
	return Q[0];	
}

namespace BezierPatchOps{

	// Caller must delete[] new patch	
	// You need to dealocate all the bezierPatches with freePatch()
	inline CORE_H
	bezierPatch* getPatch(const ParametricSurf::nurbsSurf& sel_nurbs, int& nbU, int& nbV){
		float* U_break = new float[sel_nurbs.n-sel_nurbs.p+2];
		float* V_break = new float[sel_nurbs.m-sel_nurbs.q+2];

		nbU = segs(sel_nurbs.u_knots, sel_nurbs.n+sel_nurbs.p+1, sel_nurbs.p, U_break);	
		nbV = segs(sel_nurbs.v_knots, sel_nurbs.m+sel_nurbs.q+1, sel_nurbs.q, V_break);	
		
		MathTypes::vector4f* strips = new MathTypes::vector4f[nbU * (sel_nurbs.p+1) * (sel_nurbs.m+1)];
		MathTypes::vector4f* patches = new MathTypes::vector4f[nbU * nbV * (sel_nurbs.p+1) * (sel_nurbs.q+1)];

		int got_U = 0;
		NurbsOps::decomposeSurface(	sel_nurbs.n, sel_nurbs.p, sel_nurbs.u_knots, 
									sel_nurbs.m, sel_nurbs.q, sel_nurbs.v_knots, 
									U_DIRECTION, got_U,
									sel_nurbs.control_point, strips);

		const int strip_size = (sel_nurbs.p+1) * (sel_nurbs.m+1);
		const int patch_size = (sel_nurbs.p+1) * (sel_nurbs.q+1);

		for(int s = 0; s<nbU; ++s){
			int got_V = 0;
			NurbsOps::decomposeSurface(	sel_nurbs.p, sel_nurbs.p, sel_nurbs.u_knots, 
									sel_nurbs.m, sel_nurbs.q, sel_nurbs.v_knots, 
									V_DIRECTION, got_V,
									strips + s*strip_size,
									patches + s*nbV*patch_size);

		}
		
		bezierPatch* new_patch = new bezierPatch[nbU*nbV];
		for(int i=0; i<nbU; i++){
			for(int j=0; j<nbV; j++){
				new_patch[i*nbV+j].p = sel_nurbs.p;
				new_patch[i*nbV+j].q = sel_nurbs.q;
				new_patch[i*nbV+j].u0 = U_break[i];
				new_patch[i*nbV+j].u1 = U_break[i+1];
				
				new_patch[i*nbV+j].v0 = V_break[j];
				new_patch[i*nbV+j].v1 = V_break[j+1];

				new_patch[i*nbV+j].control_point = patches + (i*nbV+j) * patch_size;
				new_patch[i*nbV+j].heap_pointer = patches;
			}

		}
		
		delete[] strips;
		delete[] U_break;
		delete[] V_break;
		
		return new_patch;

	}
	
	// Free all the patches
	CORE_H inline
	void freePatches(bezierPatch* patch){
		if(!patch) return;
		
		if(patch[0].heap_pointer)
			delete[] patch[0].heap_pointer;
		delete[] patch;
	}

	// From the formula of the Efficient Bounding of displaced Bezier Patches
	// Returned the bounds of the patch
	CORE_H inline
	MathExtra::bounds3f bounds(bezierPatch patch){
		float min[3] = {INFINITY, INFINITY, INFINITY};
		float max[3] = {-INFINITY, -INFINITY, -INFINITY};

		for(int i=0; i< (patch.p+1)*(patch.q+1); i++){
			MathTypes::vector4f sel_crp = patch.control_point[i];

			for(int j=0; j<3; j++){
				min[j] = fminf(min[j], sel_crp.v[j]/sel_crp.v[3]);
				max[j] = fmaxf(max[j], sel_crp.v[j]/sel_crp.v[3]);
			}
		}
			
		return {	min[0], min[1], min[2],
					max[0], max[1], max[2],};
	}
	
	// Evaluate bezier patch and calculates the partial derivatives
	CORE_HD inline 
	void BezierPatchEval(const MathTypes::vector4f* P, int p, int q, float u, float v, MathTypes::vector4f* ders, MathTypes::vector4f& point){

		// Calculating B_u
		MathTypes::vector4f curve_cr_p[NURBS_MAX_DEGREE+1];
		for(int i=0; i<=p; i++){
			curve_cr_p[i] = deCasteljau1(&P[i*(q+1)], q, v); // P is already contiguous so no need to isolate curve variable in here
		}
		ders[0] = deCasteljauDers(curve_cr_p, p, u);
		point = deCasteljau1(curve_cr_p, p, u);
		
		for(int j=0; j<=q; j++){
			
			MathTypes::vector4f internal_cr_p[NURBS_MAX_DEGREE+1];
			for(int i=0; i<=p; i++)
				internal_cr_p[i] = P[i*(q+1)+j];
			curve_cr_p[j] = deCasteljau1(internal_cr_p, p, u); // P is already contiguous so no need to isolate curve variable in here
		}
		ders[1] = deCasteljauDers(curve_cr_p, q, v);
	}
	

	// Inspired by the pseudo-code on Ray Tracing NURBS Surfaces using CUDA
	// The best initial guess for u,v in the case of the patches is 
	// (u_1 - u_0)/2 which is the center of the patch, repeat for v
	CUDA_FINL CORE_D 
	newtonRhapsonPayload newtonRhapson(MathExtra::ray_plane& ray_plane, bezierPatch& patch, float initial_u, float initial_v, float threshold1, float threshold2, uint seed){
		float u_it = initial_u;
		float v_it = initial_v;
		float error_prev = INFINITY;

		using namespace VectorOps;
	
		for(int i=1; i<=NEWTON_MAX_ITERATIONS; i++){
			MathTypes::vector4f ders[2];
			MathTypes::vector4f point;
			
			BezierPatchEval(patch.control_point, patch.p, patch.q, u_it, v_it, ders, point);
			
			float invW = 1.0f/point.v[3];
			MathTypes::vector3f S = {	point.v[0]*invW,
										point.v[1]*invW,
										point.v[2]*invW}; 

			MathTypes::vector3f Su = {	(ders[0].v[0] - S.v[0]*ders[0].v[3]) * invW,
										(ders[0].v[1] - S.v[1]*ders[0].v[3]) * invW,
										(ders[0].v[2] - S.v[2]*ders[0].v[3]) * invW};

			MathTypes::vector3f Sv = {	(ders[1].v[0] - S.v[0]*ders[1].v[3]) * invW,
										(ders[1].v[1] - S.v[1]*ders[1].v[3]) * invW,
										(ders[1].v[2] - S.v[2]*ders[1].v[3]) * invW};

			MathTypes::vector2f roots = {	dotProduct(ray_plane.n1, S)+ray_plane.d1,
											dotProduct(ray_plane.n2, S)+ray_plane.d2};
			
			float error = fabsf(roots.v[0]) + fabsf(roots.v[1]);
			if(error < threshold2)
				return {1,S,Su,Sv,u_it,v_it};
			
//			if(error > error_prev)
//				return {0};


			error_prev = error;
			MathTypes::matrix2x2f jacobian = {	dotProduct(ray_plane.n1, Su),
												dotProduct(ray_plane.n2, Su),
												dotProduct(ray_plane.n1, Sv),
												dotProduct(ray_plane.n2, Sv)};

			// These parts are conditional, but they are kind of heavy calculation
			// So every thread will calculate them to lower the chance of warp diverg
			float rand1 = MathUtils::randf(seed);
			float rand2 = MathUtils::randf(seed);
			
			// I took this way of calculating if it is singular from Martin et Al.
			if(fabsf(MatrixOps::determinant(jacobian)) < threshold1){
				u_it = u_it+0.1f*((initial_u-u_it)*rand1);
				v_it = v_it+0.1f*((initial_v-v_it)*rand2);
			}
			else{
				// Reusing the same variable to save register space
				jacobian = MatrixOps::inverse(jacobian);
				MathTypes::vector2f j_inv_rt = MatrixOps::multiply(jacobian, roots);
				u_it = u_it - j_inv_rt.v[0];
				v_it = v_it - j_inv_rt.v[1];

			}
			//if (u_it < 0.0f || u_it > 1.0f || v_it < 0.0f || v_it > 1.0f)
			//	return {0};
			// Clamping achieve much better results than simply exiting
			u_it = fminf(fmaxf(u_it, 0.0f), 1.0f);
			v_it = fminf(fmaxf(v_it, 0.0f), 1.0f);
		}
		return {0};	
	}

	// Fix Shadow acne 
	// From Reshetov's ray-bilinear-patch intersection
	// Because of the refiniment of the patches, a lot of the patches are "flat" which means we can calculate 
	// this function to provide a better initial u and v for the newton-rhapson
	CUDA_FINL CORE_D 
	 MathTypes::vector2f biliniearRay(	MathTypes::vector3f q00, MathTypes::vector3f q10,
											MathTypes::vector3f q11, MathTypes::vector3f q01,
											MathTypes::vector3f orig, MathTypes::vector3f dir)
	{
		MathTypes::vector3f e10 = VectorOps::subtract(q10, q00);
		MathTypes::vector3f e11 = VectorOps::subtract(q11, q10);
		MathTypes::vector3f e00  = VectorOps::subtract(q01, q00);
		
		MathTypes::vector3f qn = VectorOps::cross(e10, VectorOps::subtract(q01, q11));

		q00 = VectorOps::subtract(q00, orig);
		q10 = VectorOps::subtract(q10, orig);
		
		float a = VectorOps::dotProduct(VectorOps::cross(q00, dir), e00);
		float c = VectorOps::dotProduct(qn, dir);
		float b = VectorOps::dotProduct(VectorOps::cross(q10, dir), e11);

		b -= a+c;

		float det = b*b - 4*a*c;
		if (det < 0) return {0.5f, 0.5f};
		det = sqrtf(det);
		float ux[2];

		float u=0.5f, v=0.5f;
		if (c==0){
			ux[0] = -a/b; ux[1] = -1;
		}
		else{
			ux[0] = (-b - copysignf(det, b))/2;
			ux[1] = a/ux[0];
			ux[0] /= c;
		}
		
		float best_t = INFINITY;
		for(int i=0; i<2; i++){
			if(ux[i] < 0 || ux[i] > 1 ) continue;
			MathTypes::vector3f pa = MathUtils::lerp<MathTypes::vector3f>(q00, q10, ux[i]);
			MathTypes::vector3f pb = MathUtils::lerp(e00, e11, ux[i]);

			MathTypes::vector3f n = VectorOps::cross(dir, pb);

			det = VectorOps::dotProduct(n,n);
			float tx = VectorOps::dotProduct(n, pb);
			float vx = VectorOps::dotProduct(n, dir);

			if(tx > 0 && 0 <= vx && vx <= det && tx < best_t){
				u = ux[i]; 
				v = vx/det;
				best_t = tx;
			}
		}
		return {u,v};

	}


	CUDA_FINL CORE_D
	void evalPointNormal(	const bezierPatch& patch, float u, float v, 
							MathTypes::vector3f dir, MathTypes::vector3f& S,
							MathTypes::vector3f& Su, MathTypes::vector3f& Sv){
		MathTypes::vector4f ders[2], point;
		BezierPatchOps::BezierPatchEval(patch.control_point, patch.p, patch.q, u, v, ders, point);
		float iw=1.f/point.v[3];
		S = {point.v[0]*iw, point.v[1]*iw, point.v[2]*iw};
		Su={(ders[0].v[0]-S.v[0]*ders[0].v[3])*iw,(ders[0].v[1]-S.v[1]*ders[0].v[3])*iw,(ders[0].v[2]-S.v[2]*ders[0].v[3])*iw};
		Sv={(ders[1].v[0]-S.v[0]*ders[1].v[3])*iw,(ders[1].v[1]-S.v[1]*ders[1].v[3])*iw,(ders[1].v[2]-S.v[2]*ders[1].v[3])*iw};
	}
}

#endif
