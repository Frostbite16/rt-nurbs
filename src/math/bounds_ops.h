#ifndef BOUNDS_OPS_H
#define BOUNDS_OPS_H

#include <cmath>
#include "math_types.h"
#include "vector_ops.h"
#include "../core/core_defines.h"

namespace BoundsOps{
	
	// Initialize bounds2f 
	inline CORE_HD
	MathExtra::bounds2f create(float x1, float y1, float x2, float y2){
		MathExtra::bounds2f bound;
		bound.p_min_x = fminf(x1, x2);
		bound.p_min_y = fminf(y1, y2);
		bound.p_max_x = fmaxf(x1, x2);
		bound.p_max_y = fmaxf(y1, y2);
		return bound;
	}

	
	// Initialize bounds3f
	// No branching
	inline CORE_HD
	MathExtra::bounds3f create(float x1, float y1, float z1, float x2, float y2, float z2){
		MathExtra::bounds3f bound;
		bound.p_min_x = fminf(x1, x2);
		bound.p_min_y = fminf(y1, y2);
		bound.p_min_z = fminf(z1, z2);
		
		bound.p_max_x = fmaxf(x1, x2);
		bound.p_max_y = fmaxf(y1, y2);
		bound.p_max_z = fmaxf(z1, z2);	
		return bound;
	}

	// Returns the centroid of an bound
	inline constexpr CORE_HD
	MathTypes::point3f centroid(MathExtra::bounds3f bd){
		return MathTypes::point3f{ 	.5f*(bd.p_min_x + bd.p_max_x),
									.5f*(bd.p_min_y + bd.p_max_y),
									.5f*(bd.p_min_z + bd.p_max_z)
									 };
	}

	// Returns a new bound that encapsulates the two given bounds
	inline CORE_HD
	MathExtra::bounds3f boundUnion(const MathExtra::bounds3f& b1, const MathExtra::bounds3f& b2){
		MathExtra::bounds3f new_bound;
		new_bound.p_min_x = fminf(b1.p_min_x, b2.p_min_x);
		new_bound.p_min_y = fminf(b1.p_min_y, b2.p_min_y);
		new_bound.p_min_z = fminf(b1.p_min_z, b2.p_min_z);

		new_bound.p_max_x = fmaxf(b1.p_max_x, b2.p_max_x);
		new_bound.p_max_y = fmaxf(b1.p_max_y, b2.p_max_y);
		new_bound.p_max_z = fmaxf(b1.p_max_z, b2.p_max_z);	
		return new_bound;
	}
	
	// Grow 3 dimensional bound to parameter point
	inline CORE_HD
	MathExtra::bounds3f grow(const MathTypes::point3f point, MathExtra::bounds3f& bound){
		return MathExtra::bounds3f{ fminf(bound.p_min_x, point.v[0]),
									fminf(bound.p_min_y, point.v[1]),
									fminf(bound.p_min_z, point.v[2]),
									fmaxf(bound.p_max_x, point.v[0]),
                                    fmaxf(bound.p_max_y, point.v[1]),                               
									fmaxf(bound.p_max_z, point.v[2])};
	}

	inline CORE_HD
	float bound_area(MathExtra::bounds3f bound){
		MathTypes::vector3f e = VectorOps::subtract(	MathTypes::point3f{bound.p_max_x, bound.p_max_y, bound.p_max_z},
														MathTypes::point3f{bound.p_min_x, bound.p_min_y, bound.p_min_z});	
		return e.v[0] * e.v[1] + e.v[1] * e.v[2] + e.v[2] * e.v[0];
	}

	// Given an 3 dimensional bound and an ray, calculate the intersection 
	// of the given ray and return the t-value of that intersection (INFINITY if there is not hit)
	// Inline function, so the compiler just puts this into the kernel instead of calling an function
	// DEVICE ONLY FUNCTION
	CUDA_FINL CORE_D
	float intersectBound(float bounds_min_x, float bounds_max_x,
						 float bounds_min_y, float bounds_max_y,
						 float bounds_min_z, float bounds_max_z,
						 float ray_o_x, float ray_o_y, float ray_o_z,
						 float ray_d_rec_x, float ray_d_rec_y, float ray_d_rec_z){
			
		float tmin = 0, tmax = INFINITY;

		// Calculating the reciprocals of the ray
		// if ray->dir = 0 reciprocal goes to infinity, that is intended behaviour

		float t_min_x = (bounds_min_x - ray_o_x) * ray_d_rec_x;
		float t_max_x = (bounds_max_x - ray_o_x) * ray_d_rec_x; 
		float t_min_y = (bounds_min_y - ray_o_y) * ray_d_rec_y;
		float t_max_y = (bounds_max_y - ray_o_y) * ray_d_rec_y; 
		float t_min_z = (bounds_min_z - ray_o_z) * ray_d_rec_z;
		float t_max_z = (bounds_max_z - ray_o_z) * ray_d_rec_z; 
	
		float near_x = fminf(t_min_x, t_max_x);	
		float far_x = fmaxf(t_min_x, t_max_x);	
		float near_y = fminf(t_min_y, t_max_y);	
		float far_y = fmaxf(t_min_y, t_max_y);	
		float near_z = fminf(t_min_z, t_max_z);	
		float far_z = fmaxf(t_min_z, t_max_z);	

		tmin = fmaxf(near_z, fmaxf(near_y, fmaxf(near_x, tmin)));
		tmax = fminf(far_z, fminf(far_y, fminf(far_x, tmax)));

		return (tmin <= tmax && tmax > 0) ? tmin : INFINITY;

	}	
}

#endif

