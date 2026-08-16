#include "soa.h"
#include "../core/camera.h"
#include "../math/math_types.h"
#include "../core/core_defines.h"
#include "generateRaysKernel.h"


#include <cstdio>

CUDA_CONST camera_gpu_payload const_camera_config;


// GPU Kernel that calculates the origin and direction of every ray in ray queue
// Template is for choosing projection without branching
template<typename T>
CUDA_GLOB void calculateRays(SoAQueue::ray_queue ray_queue, int screen_res_x, int screen_res_y){
	
	int idx_x = threadIdx.x + blockIdx.x*blockDim.x;
	int idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	
	if(idx_x >= screen_res_x || idx_y >= screen_res_y)
		return;
	
	MathExtra::ray new_ray = CameraOps::calculateRay<T>(const_camera_config, idx_x, idx_y);

	ray_queue.origin_x[idx_y * screen_res_x + idx_x] = new_ray.origin.v[0];
	ray_queue.origin_y[idx_y * screen_res_x + idx_x] = new_ray.origin.v[1];
	ray_queue.origin_z[idx_y * screen_res_x + idx_x] = new_ray.origin.v[2];

	ray_queue.direction_x[idx_y * screen_res_x + idx_x] = new_ray.direction.v[0];
	ray_queue.direction_y[idx_y * screen_res_x + idx_x] = new_ray.direction.v[1];
	ray_queue.direction_z[idx_y * screen_res_x + idx_x] = new_ray.direction.v[2];
	
	ray_queue.idx_interact[idx_y * screen_res_x + idx_x] = -1;
}


// Wrapper for calculateRays kernel
void calculateRaysWrapper(SoAQueue::ray_queue& ray_queue, projection proj, camera_config& camera){

	// copy payload to const space in device memory
	camera_gpu_payload buffer_payload = {camera.world_from_camera, camera.camera_from_raster};
	cudaMemcpyToSymbol(const_camera_config, &buffer_payload, sizeof(camera_gpu_payload));

	// Make an 2D thread array in the gpu
	// this has no performance cost
	dim3 thread_2D(16,16);
	dim3 block_2D(
		(camera.screen_resolution_x + 15) / 16,
		(camera.screen_resolution_y + 15) / 16
			);	

	printf("block X,Y dim: %d,%d \n", block_2D.x, block_2D.y);	
	cudaGetLastError();

	// Run kernel per projection
	// cuda synchronize is needed to make the cpu wait for the gpu workload to finish
	if(proj == orthographic){
		calculateRays<orthographic_projection><<<block_2D,thread_2D>>>(ray_queue, 
				camera.screen_resolution_x, camera.screen_resolution_y);
		cudaDeviceSynchronize();
	}
	if(proj == perspective){
		calculateRays<perspective_projection><<<block_2D,thread_2D>>>(ray_queue, 
				camera.screen_resolution_x, camera.screen_resolution_y);	
		cudaDeviceSynchronize();
	}

}
