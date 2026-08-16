#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "math_types.h"
#include "vector_ops.h"
#include "../core/core_defines.h"

// Namespace for transform matrices
namespace Transform{

	// Returns an usable scale matrix with the scalar arguments
	// Function can be run on Both host and the device
	// can be evaluated at compile time
	inline constexpr CORE_HD 
	MathTypes::matrix4x4f scale(float sx, float sy, float sz){
		MathTypes::matrix4x4f mat = {0};
		mat.m[0] = sx;
		mat.m[5] = sy;
		mat.m[10] = sz;
		mat.m[15] = 1.0f;

		return mat;
	} 

	// Returns an usable tranlation matrix with the scalar arguments
	// Function can be run on Both host and the device
	// can be evaluated at compile time
	inline constexpr CORE_HD 
	MathTypes::matrix4x4f translate(float tx, float ty, float tz){
		MathTypes::matrix4x4f mat = {0};
		mat.m[0] = 1.0f;
		mat.m[5] = 1.0f;
		mat.m[10] = 1.0f;
		mat.m[15] = 1.0f;

		// Translation part
		mat.m[12] = tx;
		mat.m[13] = ty;
		mat.m[14] = tz;

		return mat;
	}

	// Transformation exclusive for placing the camera
	// Returns an transform matrix that transform any point from camera space to world space
	inline constexpr CORE_HD 
	MathTypes::matrix4x4f lookAt(const MathTypes::point3f& pos,
								 const MathTypes::point3f& look, 
								 const MathTypes::vector3f& up){
		using namespace VectorOps;

		MathTypes::vector3f dir = normalize(subtract(look, pos));
		MathTypes::vector3f right = normalize(cross(up, dir));
		MathTypes::vector3f newUp = cross(dir, right);

		return MathTypes::matrix4x4f{right.v[0], right.v[1], right.v[2], 0.0f,
									 newUp.v[0], newUp.v[1], newUp.v[2], 0.0f,
									 dir.v[0],	 dir.v[1],	 dir.v[2], 	 0.0f,
									 pos.v[0], 	 pos.v[1], 	 pos.v[2],	 1.0f};
	}	
}

#endif
