#ifndef MATH_UTILS
#define MATH_UTILS

#include "../constants.h"
#include "../core/core_defines.h"
#include "math_types.h"
#include "vector_ops.h"
#include <cmath>
#include <cstdint>


CUDA_FINL CORE_D  
uint32_t pcg_hash(uint32_t x) {
   	x = x * 747796405u + 2891336453u;
  	uint32_t w = ((x >> ((x >> 28u) + 4u)) ^ x) * 277803737u;
   	return (w >> 22u) ^ w;
}
namespace MathUtils{
	
	CUDA_FINL constexpr CORE_HD
	float radians(float degree){
		return (PI / 180) * degree;
	} 

	// Implementation from:
	// https://github.com/dangets/cuda_examples/blob/master/clamp_function.cu
	CUDA_FINL constexpr CORE_HD
	float fclamp(float val, float v_max, float v_min){
		return fminf(fmaxf(val, v_min), v_max);
	}
	
	// Get the plane representation of the ray
	CUDA_FINL CORE_HD
	MathExtra::ray_plane getRayPlanes(MathExtra::ray ray){
		const float ax = fabsf(ray.direction.v[0]);
		const float ay = fabsf(ray.direction.v[1]);
		const float az = fabsf(ray.direction.v[2]);
		const bool x_dominant = (ax > ay) && (ax > az);
		

		MathTypes::vector3f base = (x_dominant) ? 	MathTypes::vector3f{ray.direction.v[0], ray.direction.v[1], 0} : 
													MathTypes::vector3f{ray.direction.v[1], ray.direction.v[2], 0};
		
		float den = VectorOps::reciprocal_length(base);
		
		MathTypes::vector3f ort_vec_n1 = (x_dominant) ? 	VectorOps::multiply(den, MathTypes::vector3f{ray.direction.v[1], -ray.direction.v[0],0}) :
															VectorOps::multiply(den, MathTypes::vector3f{0,ray.direction.v[2], -ray.direction.v[1]});

		MathTypes::vector3f ort_vec_n2 = VectorOps::cross(ort_vec_n1, VectorOps::multiply(ray.direction, 1.0f/VectorOps::length(ray.direction)));

		float orig_dist_1 = -VectorOps::dotProduct(ort_vec_n1, ray.origin);
		float orig_dist_2 = -VectorOps::dotProduct(ort_vec_n2, ray.origin);
		return {ort_vec_n1, ort_vec_n2, orig_dist_1, orig_dist_2};
	}		
	

	template <typename T, typename F>
	CUDA_FINL CORE_HD
	T lerp(const T& v0, const T& v1, F t);

	template <>
	CUDA_FINL CORE_HD
	MathTypes::vector3f lerp(const MathTypes::vector3f& v0, const MathTypes::vector3f& v1, float t){
		MathTypes::vector3f r1 = VectorOps::multiply((1-t),v0); 
		MathTypes::vector3f r2 = VectorOps::multiply(t,v1); 

		return VectorOps::sum(r1,r2);
	}

	CUDA_FINL CORE_D
	float randf(uint32_t& seed) {
    	seed = pcg_hash(seed);
    	return seed * 2.3283064e-10f;   // / 2^32, gives [0,1)
	}
}

#endif

