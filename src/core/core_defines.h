#pragma once
// Without this the c++ compiler complaims about cuda code in header files
#ifdef __CUDACC__
	#define CORE_HD __host__ __device__
	#define CORE_H __host__ 
	#define CORE_D __device__ 
	#define CUDA_CONST __constant__
	#define CUDA_GLOB __global__
	#define CUDA_FINL __forceinline__
#else
	#define CORE_HD
	#define CORE_H 
	#define CORE_D 
	#define CUDA_CONST 
	#define CUDA_GLOB 
	#define CUDA_FINL inline
#endif


