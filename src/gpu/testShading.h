#ifndef TESTSHADING_H
#define TESTSHADING_H

#include "../core/core_defines.h"
#include "soa.h"
#include "interact.h"

CORE_H
void collisionTestShadingWrapper(	SoAQueue::ray_queue ray_queue, SoAQueue::frame_buffer screen,
									int res_x, int res_y);

CORE_H
void normalsTestShadingWrapper(	SoAQueue::ray_queue ray_queue, SoAQueue::frame_buffer screen,
								surfaceInteraction* patches, int res_x, int res_y);

#endif
