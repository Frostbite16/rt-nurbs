#include "soa.h"
#include "testShading.h"
#include "interact.h"
#include "../core/core_defines.h"
#include "../math/math_types.h"
#include "../math/vector_ops.h"
#include <cstdio>

CUDA_GLOB
void redTestShading(	SoAQueue::ray_queue ray_queue, SoAQueue::frame_buffer screen,
						int res_x, int res_y ){
	
	int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
	int idx_y = threadIdx.y + blockIdx.y*blockDim.y;

	if(idx_x >= res_x || idx_y >= res_y)
		return;
	
	int idx = idx_y * res_x + idx_x;
	
	screen.r[idx] = (ray_queue.idx_interact[idx] > -1);
	screen.g[idx] = 0.0f;
	screen.b[idx] = 0.0f;
}


CUDA_GLOB
void normalsTestShading(	SoAQueue::ray_queue ray_queue, SoAQueue::frame_buffer screen,
							surfaceInteraction* interact, int res_x, int res_y ){
	
	int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
	int idx_y = threadIdx.y + blockIdx.y*blockDim.y;

	if(idx_x >= res_x || idx_y >= res_y)
		return;
	
	int idx = idx_y * res_x + idx_x;
		
	screen.r[idx] = 0.0f;
	screen.g[idx] = 0.0f;
	screen.b[idx] = 0.0f;

	int interact_idx = ray_queue.idx_interact[idx];
	if(interact_idx > -1){

		MathTypes::vector3f normal_vector = VectorOps::normalize(VectorOps::cross(interact[interact_idx].dpdv, interact[interact_idx].dpdu));
		
		screen.r[idx] = (normal_vector.v[0] + 1)*.5f;
		screen.g[idx] = (normal_vector.v[1] + 1)*.5f;
		screen.b[idx] = (normal_vector.v[2] + 1)*.5f;
	}
	


}


CORE_H
void collisionTestShadingWrapper(	SoAQueue::ray_queue ray_queue, SoAQueue::frame_buffer screen,
									int res_x, int res_y){

	int thread_x_quant = 16;
	int thread_y_quant = 16;
	dim3 thread_2D(thread_x_quant,thread_y_quant);
	dim3 block_2D(
		(res_x + thread_x_quant - 1) / thread_x_quant,
		(res_y + thread_y_quant - 1) / thread_y_quant);

	
	redTestShading<<<block_2D, thread_2D>>>( ray_queue, screen, res_x, res_y);
	cudaDeviceSynchronize();
	
	printf("Test shading done\n");
}


CORE_H
void normalsTestShadingWrapper(	SoAQueue::ray_queue ray_queue, SoAQueue::frame_buffer screen,
								surfaceInteraction* interact, int res_x, int res_y){

	int thread_x_quant = 16;
	int thread_y_quant = 16;
	dim3 thread_2D(thread_x_quant,thread_y_quant);
	dim3 block_2D(
		(res_x + thread_x_quant - 1) / thread_x_quant,
		(res_y + thread_y_quant - 1) / thread_y_quant);

	
	normalsTestShading<<<block_2D, thread_2D>>>(ray_queue, screen, interact, res_x, res_y);
	cudaDeviceSynchronize();
	
	printf("Test shading done\n");
}



