#ifndef GENERATE_RAYS_KERNEL_H
#define GENERATE_RAYS_KERNEL_H
#include "../core/camera.h"
#include "soa.h"

void calculateRaysWrapper(SoAQueue::ray_queue& ray_queue, projection proj, camera_config& camera);


#endif 
