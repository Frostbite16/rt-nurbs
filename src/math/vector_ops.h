#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include "math_types.h"
#include "../core/core_defines.h"
#include <cmath>

namespace VectorOps{
	

	// Vector3 ops
	
	// Returns sum of 3 dimensional float vectors 
	inline constexpr CORE_HD 
	MathTypes::vector3f sum(const MathTypes::vector3f& vec1, const MathTypes::vector3f& vec2){
		return MathTypes::vector3f{vec1.v[0]+vec2.v[0], vec1.v[1]+vec2.v[1], vec1.v[2]+vec2.v[2]};
	}


	// Returns sum of 3 dimensional float vector an point 
	inline constexpr CORE_HD 
	MathTypes::vector3f sum(const MathTypes::point3f& p1, const MathTypes::vector3f& vec2){
		return MathTypes::vector3f{p1.v[0]+vec2.v[0], p1.v[1]+vec2.v[1], p1.v[2]+vec2.v[2]};
	}

	// Returns subtraction of 3 dimensional float vectors 
	inline constexpr CORE_HD 
	MathTypes::vector3f subtract(const MathTypes::vector3f& vec1, const MathTypes::vector3f& vec2){
		return MathTypes::vector3f{vec1.v[0]-vec2.v[0], vec1.v[1]-vec2.v[1], vec1.v[2]-vec2.v[2]};
	}

	// Returns subtraction of 3 dimensional float points (direction)
	inline constexpr CORE_HD
	MathTypes::vector3f subtract(const MathTypes::point3f& vec1, const MathTypes::point3f& vec2){
		return MathTypes::vector3f{vec1.v[0]-vec2.v[0], vec1.v[1]-vec2.v[1], vec1.v[2]-vec2.v[2]};
	}

	// Returns product of 3 dimensional float vector and float scalar
	inline constexpr CORE_HD
	MathTypes::vector3f multiply(const MathTypes::vector3f& vec, const float s){
		return MathTypes::vector3f{vec.v[0]*s, vec.v[1]*s, vec.v[2]*s};
	}

	inline constexpr CORE_HD
	MathTypes::vector3f multiply(const float s, const MathTypes::vector3f& vec){
		return multiply(vec, s);
	}

	// Returns the squared lenght of an 3 dimensional float vector
	inline constexpr CORE_HD
	float lengthSquare(const MathTypes::vector3f& vec){
		return vec.v[0]*vec.v[0] + vec.v[1]*vec.v[1] + vec.v[2]*vec.v[2];
	}
	
	// Returns the dot product of two 3 dimensional float vectors
	inline constexpr CORE_HD
	float dotProduct(const MathTypes::vector3f& vec1, const MathTypes::vector3f& vec2){
		return vec1.v[0]*vec2.v[0] + vec1.v[1]*vec2.v[1] + vec1.v[2]*vec2.v[2];
	}
	
	// Returns the dot product of two 3 dimensional float vectors
	inline constexpr CORE_HD
	float dotProduct(const MathTypes::vector3f& vec1, const MathTypes::vector4f& vec2){
		return vec1.v[0]*vec2.v[0] + vec1.v[1]*vec2.v[1] + vec1.v[2]*vec2.v[2];
	}

	// Returns the dot product of two 3 dimensional float points
	inline constexpr CORE_HD
	float dotProduct(const MathTypes::point3f& p1, const MathTypes::point3f& p2){
		return p1.v[0]*p2.v[0] + p1.v[1]*p2.v[1] + p1.v[2]*p2.v[2];
	}
	
	// Return lenght of the vector
	inline constexpr CORE_HD
	float length(const MathTypes::vector3f& vec){
		return sqrt(vec.v[0]*vec.v[0] + vec.v[1]*vec.v[1] + vec.v[2]*vec.v[2]);
	}

	// Returns sum of 4 dimensional float vector an point 
	inline constexpr CORE_HD 
	MathTypes::vector4f sum(const MathTypes::vector4f& vec1, const MathTypes::vector4f& vec2){
		return MathTypes::vector4f{vec1.v[0]+vec2.v[0], vec1.v[1]+vec2.v[1], vec1.v[2]+vec2.v[2], vec1.v[3]+vec2.v[3]};
	}

	// Returns subtraction of 3 dimensional float vectors 
	inline constexpr CORE_HD 
	MathTypes::vector4f subtract(const MathTypes::vector4f& vec1, const MathTypes::vector4f& vec2){
		return MathTypes::vector4f{vec1.v[0]-vec2.v[0], vec1.v[1]-vec2.v[1], vec1.v[2]-vec2.v[2], vec1.v[3] - vec2.v[3]};
	}


// Fix weird compile error
#ifdef __CUDA_ARCH__
		// Return the reciprocal lenght of the vector 
		inline constexpr CORE_D 
		float reciprocal_length(const MathTypes::vector3f& vec){
			return rsqrtf(vec.v[0]*vec.v[0] + vec.v[1]*vec.v[1] + vec.v[2]*vec.v[2]);
		}
	
#else 

	// Return the reciprocal lenght of the vector 
	inline constexpr CORE_H 
	float reciprocal_length(const MathTypes::vector3f& vec){
		return 1.0f/sqrtf(vec.v[0]*vec.v[0] + vec.v[1]*vec.v[1] + vec.v[2]*vec.v[2]);
	}

#endif
	
	// Return normalized 3 dimensional vector
	inline constexpr CORE_HD
	MathTypes::vector3f normalize(const MathTypes::vector3f& vec){
		return multiply(vec, reciprocal_length(vec));
	}

	// Return cross product vector	
	inline constexpr CORE_HD
	MathTypes::vector3f cross(const MathTypes::vector3f& vec1, const MathTypes::vector3f& vec2){
		return MathTypes::vector3f{	vec1.v[1]*vec2.v[2] - vec1.v[2]*vec2.v[1],
									vec1.v[2]*vec2.v[0] - vec1.v[0]*vec2.v[2], 
									vec1.v[0]*vec2.v[1] - vec1.v[1]*vec2.v[0]};										
	}
	
	// Returns product of 4 dimensional float vector and float scalar
	inline constexpr CORE_HD
	MathTypes::vector4f multiply(const MathTypes::vector4f& vec, const float s){
		return MathTypes::vector4f{vec.v[0]*s, vec.v[1]*s, vec.v[2]*s, vec.v[3]*s};
	}

	inline constexpr CORE_HD
	MathTypes::vector4f multiply(const float s, const MathTypes::vector4f& vec){
		return multiply(vec, s);
	}
}



#endif


