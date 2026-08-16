#ifndef MEMORY_H
#define MEMORY_H

#include <cstdint>
#include <cstdio>
#include <cstdlib>

// Error checking CUDA
#define CUDA_CHECK(expr_to_check) do {            	\
    cudaError_t result  = expr_to_check;          	\
    if(result != cudaSuccess)                     	\
    {                                             	\
        fprintf(stderr,                           	\
                "CUDA Runtime Error: %s:%i:%d = %s\n", \
                __FILE__,                         	\
                __LINE__,                         	\
                result,								\
                cudaGetErrorString(result));      	\
    }                                             	\
} while(0)

#include <cuda_runtime.h>
// Alocate sizeof(T)*count bytes in the CPU
template <typename T>
T* CudaHostAlloc(size_t count){
	T* addr;
	CUDA_CHECK(cudaMallocHost(&addr, sizeof(T)*count));
	return addr;
}

// Alocate sizeof(T)*count bytes in the GPU 
template <typename T>
T* CudaDevAlloc(size_t count){
	if(count > SIZE_MAX / sizeof(size_t)) {
		abort();
	}

	T* addr;
	printf("Alloc %zu bytes\n", sizeof(T)*count);
	CUDA_CHECK(cudaMalloc(&addr, sizeof(T)*count));
	return addr;
}

// Alocate sizeof(T)*count bytes of managed memory between GPU and CPU 
template <typename T>
T* CudaAllocManaged(size_t count){
	T* addr;
	CUDA_CHECK(cudaMallocManaged(&addr, sizeof(T)*count));
	return addr;
}

inline void CudaFree(void* addr){
	if(addr){
		CUDA_CHECK(cudaFree(addr));
		addr = nullptr;
	}
}

inline void CudaCpuFree(void* addr){
	if(addr){
		CUDA_CHECK(cudaFreeHost(addr));
		addr = nullptr;
	}
}

#endif

