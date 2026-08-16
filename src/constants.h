#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "math/math_types.h"
#include "core/core_defines.h"
#include <cmath>

// Compile time definitions
// Inline is needed to force compiler to share memory adress across files

// Identity of an 4 by 4 matrix
CORE_HD inline constexpr 
MathTypes::matrix4x4f IDENTITY_MATRIX_4X4 = {	1,0,0,0,
				  												0,1,0,0,
																0,0,1,0,
																0,0,0,1	};
// PI
CORE_HD inline constexpr 
float PI = 3.14159265;


// Maximum value of float
CORE_HD inline constexpr 
float MAX_FLOAT = MAXFLOAT;

// Maximum value of float
CORE_HD inline constexpr 
float INF_FLOAT = INFINITY;

CORE_HD inline constexpr
int NURBS_MAX_DEGREE = 4;

CORE_H inline constexpr
int U_DIRECTION = 0;

CORE_H inline constexpr
int V_DIRECTION = 1;

CORE_HD inline constexpr
int NEWTON_MAX_ITERATIONS = 7;

CORE_HD inline constexpr
int NEWTON_MAX_STALL = 2;

#endif
