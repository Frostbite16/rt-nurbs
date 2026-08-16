#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "math_types.h"
#include "../constants.h"
#include "../core/core_defines.h"

namespace MatrixOps{
	

	// 2x2 operations
	
	//	determinant of an 2x2 matrix 
	inline constexpr CORE_HD
	float determinant(const MathTypes::matrix2x2f& mat){
		return mat.m[0]*mat.m[3] - mat.m[2]*mat.m[1]; 
	}
	
	// inverse of an 2x2 matrix 
	inline constexpr CORE_HD
	MathTypes::matrix2x2f inverse(MathTypes::matrix2x2f& mat){
		MathTypes::matrix2x2f adj_m = {	mat.m[3], -mat.m[1],
										-mat.m[2], mat.m[0]};
		float m_det = determinant(mat);
		return {adj_m.m[0]/m_det, adj_m.m[1]/m_det, adj_m.m[2]/m_det, adj_m.m[3]/m_det};
	}

	// multiplication of an 2x2 matrix with an 2 dimensional vector
	inline constexpr CORE_HD
	MathTypes::vector2f multiply(MathTypes::matrix2x2f& mat, MathTypes::vector2f& vec){
		return {mat.m[0]*vec.v[0] + mat.m[2]*vec.v[1], mat.m[1]*vec.v[0] + mat.m[3]*vec.v[1]};
	}

	// 3x3 operations

	//	determinant of an 3x3 matrix using leibniz formula
	inline constexpr CORE_HD
	float determinant(MathTypes::matrix3x3f mat){
		return  mat.m[0]*mat.m[4]*mat.m[8] +
				mat.m[3]*mat.m[7]*mat.m[2] +
				mat.m[6]*mat.m[1]*mat.m[5] -
				mat.m[6]*mat.m[4]*mat.m[2] -
				mat.m[3]*mat.m[1]*mat.m[8] - 
				mat.m[0]*mat.m[7]*mat.m[5];
	}

	// 4x4 operations
	
	// DEPRECATED
	// determinant of an 4x4 matrix using Laplace expansion
	inline constexpr CORE_HD
	float determinant(const MathTypes::matrix4x4f& mat){
		/*
		 * Index of each element
		 * |0 4 08 12|
		 * |1 5 09 13| 
		 * |2 6 10 14|
		 * |3 7 11 15|
		 * */

		// The first element of each iteration is (-1)^(i+j) where i is the row starting at 1 and j is the column starting at 1
		// Expanding along the first column
		return	mat.m[0] * determinant(MathTypes::matrix3x3f{mat.m[5], mat.m[6], mat.m[7], mat.m[9], mat.m[10], mat.m[11], mat.m[13], mat.m[14], mat.m[15]})-
				mat.m[1] * determinant(MathTypes::matrix3x3f{mat.m[4], mat.m[6], mat.m[7], mat.m[8], mat.m[10], mat.m[11], mat.m[12], mat.m[14], mat.m[15]})+
				mat.m[2] * determinant(MathTypes::matrix3x3f{mat.m[4], mat.m[5], mat.m[7], mat.m[8], mat.m[9], mat.m[11], mat.m[12], mat.m[13], mat.m[15]})-
				mat.m[3] * determinant(MathTypes::matrix3x3f{mat.m[4], mat.m[5], mat.m[6], mat.m[8], mat.m[9], mat.m[10], mat.m[12], mat.m[13], mat.m[14]}); 
	}

	// Adjugate of an 4x4 matrix
	inline constexpr CORE_HD 
	MathTypes::matrix4x4f adjugate(const MathTypes::matrix4x4f& mat){

		// Ugly but efficient code
		// If i find a more pretty implementation that has the same performance change the code

		// This is already transposing the cofactor matrix 
		MathTypes::matrix4x4f adj_matrix = {0};
		adj_matrix.m[0] 	=  determinant(MathTypes::matrix3x3f{mat.m[5], mat.m[6], mat.m[7], mat.m[9], mat.m[10], mat.m[11], mat.m[13], mat.m[14], mat.m[15]});
		adj_matrix.m[4] 	= -determinant(MathTypes::matrix3x3f{mat.m[4], mat.m[6], mat.m[7], mat.m[8], mat.m[10], mat.m[11], mat.m[12], mat.m[14], mat.m[15]}); 
		adj_matrix.m[8] 	=  determinant(MathTypes::matrix3x3f{mat.m[4], mat.m[5], mat.m[7], mat.m[8], mat.m[9], mat.m[11], mat.m[12], mat.m[13], mat.m[15]});
		adj_matrix.m[12] 	= -determinant(MathTypes::matrix3x3f{mat.m[4], mat.m[5], mat.m[6], mat.m[8], mat.m[9], mat.m[10], mat.m[12], mat.m[13], mat.m[14]}); 
		
		adj_matrix.m[1] 	= -determinant(MathTypes::matrix3x3f{mat.m[1], mat.m[2], mat.m[3], mat.m[9], mat.m[10], mat.m[11], mat.m[13], mat.m[14], mat.m[15]}); 
		adj_matrix.m[5]		=  determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[2], mat.m[3], mat.m[8], mat.m[10], mat.m[11], mat.m[12], mat.m[14], mat.m[15]}); 
		adj_matrix.m[9]		= -determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[1], mat.m[3], mat.m[8], mat.m[9], mat.m[11], mat.m[12], mat.m[13], mat.m[15]}); 
		adj_matrix.m[13]	=  determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[1], mat.m[2], mat.m[8], mat.m[9], mat.m[10], mat.m[12], mat.m[13], mat.m[14]}); 

		adj_matrix.m[2]		=  determinant(MathTypes::matrix3x3f{mat.m[1], mat.m[2], mat.m[3], mat.m[5], mat.m[6], mat.m[7], mat.m[13], mat.m[14], mat.m[15]}); 
		adj_matrix.m[6]		= -determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[2], mat.m[3], mat.m[4], mat.m[6], mat.m[7], mat.m[12], mat.m[14], mat.m[15]}); 
		adj_matrix.m[10]	=  determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[1], mat.m[3], mat.m[4], mat.m[5], mat.m[7], mat.m[12], mat.m[13], mat.m[15]}); 
		adj_matrix.m[14]	= -determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[1], mat.m[2], mat.m[4], mat.m[5], mat.m[6], mat.m[12], mat.m[13], mat.m[14]}); 

		adj_matrix.m[3]		= -determinant(MathTypes::matrix3x3f{mat.m[1], mat.m[2], mat.m[3], mat.m[5], mat.m[6], mat.m[7], mat.m[9], mat.m[10], mat.m[11]}); 
		adj_matrix.m[7]		=  determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[2], mat.m[3], mat.m[4], mat.m[6], mat.m[7], mat.m[8], mat.m[10], mat.m[11]}); 
		adj_matrix.m[11]	= -determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[1], mat.m[3], mat.m[4], mat.m[5], mat.m[7], mat.m[8], mat.m[9], mat.m[11]}); 
		adj_matrix.m[15]	=  determinant(MathTypes::matrix3x3f{mat.m[0], mat.m[1], mat.m[2], mat.m[4], mat.m[5], mat.m[6], mat.m[8], mat.m[9], mat.m[10]}); 
		return adj_matrix;
	}

	// 4x4 matrix multiplication (Column major)
	inline constexpr CORE_HD 
	MathTypes::matrix4x4f multiply(const MathTypes::matrix4x4f& mat1, const MathTypes::matrix4x4f& mat2){
		MathTypes::matrix4x4f result = {0};
		
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				float sum = 0;
				for(int k=0; k<4; k++){
					sum += mat1.m[k*4 + i] * mat2.m[j*4 + k];
				}
				result.m[j*4 + i] = sum;
			}
		}
		return result;
	}

	// 4x4 matrix multiplication with 4 dimensional point/vector (Column major)
	template <typename T> 
	inline constexpr CORE_HD 
	T multiply(const MathTypes::matrix4x4f& mat, const T& elem){
		static_assert(sizeof(elem.v) >= 4 * sizeof(float), "multiply: T must have 4 elements");

		return T{	mat.m[0]*elem.v[0] + mat.m[4]*elem.v[1] + mat.m[8]*elem.v[2] + mat.m[12]*elem.v[3],
				 	mat.m[1]*elem.v[0] + mat.m[5]*elem.v[1] + mat.m[9]*elem.v[2] + mat.m[13]*elem.v[3],
					mat.m[2]*elem.v[0] + mat.m[6]*elem.v[1] + mat.m[10]*elem.v[2] + mat.m[14]*elem.v[3],
					mat.m[3]*elem.v[0] + mat.m[7]*elem.v[1] + mat.m[11]*elem.v[2] + mat.m[15]*elem.v[3]};

	}

	inline constexpr CORE_D 
	MathTypes::point3f multiply(const MathTypes::matrix4x4f& mat, const MathTypes::point3f& elem){

		return MathTypes::point3f{	mat.m[0]*elem.v[0] + mat.m[4]*elem.v[1] + mat.m[8]*elem.v[2] + mat.m[12],
				 					mat.m[1]*elem.v[0] + mat.m[5]*elem.v[1] + mat.m[9]*elem.v[2] + mat.m[13],
									mat.m[2]*elem.v[0] + mat.m[6]*elem.v[1] + mat.m[10]*elem.v[2] + mat.m[14]};
	} 

	inline constexpr CORE_D 
	MathTypes::vector3f multiply(const MathTypes::matrix4x4f& mat, const MathTypes::vector3f& elem){

		return MathTypes::vector3f{	mat.m[0]*elem.v[0] + mat.m[4]*elem.v[1] + mat.m[8]*elem.v[2],
				 					mat.m[1]*elem.v[0] + mat.m[5]*elem.v[1] + mat.m[9]*elem.v[2],
									mat.m[2]*elem.v[0] + mat.m[6]*elem.v[1] + mat.m[10]*elem.v[2]};
	} 



	// 4x4 matrix inverse (1/Mat) using cramer's rule (Scary!)
	inline constexpr CORE_HD 
	MathTypes::matrix4x4f inverse(const MathTypes::matrix4x4f& mat){
		
		MathTypes::matrix4x4f result_matrix = adjugate(mat);
		// Calculating determinant
		// Using transposed adjunct matrix to get the cofactor matrix
		float matrix_det = mat.m[0]  * result_matrix.m[0] + 
							mat.m[1]  * result_matrix.m[4] +
							mat.m[2]  * result_matrix.m[8] + 
							mat.m[3] * result_matrix.m[12];
		float abs_det = (matrix_det < 0.0f) ? -matrix_det : matrix_det; 
		if(abs_det < 1e-8f)
			return IDENTITY_MATRIX_4X4;

		float det_inverse = 1.0f/matrix_det;
		// implementation with for
		// if i need matrix division with scalar later i switch it up
		// GPU and CPU can unroll this loop
		for(int i=0; i<16; i++){
			result_matrix.m[i]*=det_inverse;
		}
		return result_matrix;
	}
}

#endif
