#ifndef MATHTYPES_H
#define MATHTYPES_H


namespace MathTypes {
	
	// Define a 4 by 4 matrix in column major order
	struct matrix4x4f{
		float m[16];
	};
	
	// Define a 4 by 4 matrix in column major order
	struct matrix4x4i{
		int m[16];
	};
	
	// Define a 3 by 3 matrix in column major order
	struct matrix3x3f{
		float m[9];
	};

	// Define a 2 by 2 matrix in column major order
	struct matrix2x2f{
		float m[4];
	};
	
	// Define a 4 dimensional vector of int 
	struct vector4i{int v[4];};

	// Define a 4 dimensional POINT of int
	// it is important that vectors and points remain two different structures
	struct point4i{int v[4];};


	// Define a 4 dimensional vector of floats
	struct vector4f{float v[4];};

	// Define a 3 dimensional vector of floats
	struct vector3f{float v[3];};

	// Define a 3 dimensional vector of int 
	struct vector3i{int v[3];};
	
	// Define a 2 dimensional vector of floats
	struct vector2f{float v[2];};

	// Define a 4 dimensional point of floats
	struct point4f{float v[4];};

	// Define a 3 dimensional point of floats
	struct point3f{float v[3];};

	// Define a 2 dimensional point of floats
	struct point2f{float v[2];};

}


// Not inheritenly a mathematical structure 
namespace MathExtra{
	
	// Define an 2 dimensional box defined by 2 tuples of values;
	struct bounds2f{
		float p_min_x, p_min_y;
		float p_max_x, p_max_y; 
	};

	// Define an 3 dimensional box defined by 3 tuples of values;
	struct bounds3f{
		float p_min_x, p_min_y, p_min_z;
		float p_max_x, p_max_y, p_max_z;
	};
	
	// Define an ray with origin and direction 
	struct ray{
		MathTypes::vector3f origin;
		MathTypes::vector3f direction;
	};

	//plane definition of an ray
	struct ray_plane{
		MathTypes::vector3f n1;
		MathTypes::vector3f n2;
		float d1;
		float d2;
	};

	struct plane{
		MathTypes::vector3f n;
		float d;
	};


}

#endif

