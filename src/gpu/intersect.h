#ifndef INTERSECT_H
#define INTERSECT_H

#include "../core/bvh.h"
#include "../math/bezier_patch.h"
#include "soa.h"
#include "interact.h"

void calculateIntersectWrapper(	SoAQueue::ray_queue& ray_queue, BVHNode* nodes, bezierPatch* patches,
								surfaceInteraction* interact, int screen_res_x, int screen_res_y,
								float flat_threshold, float cam_k, float* delta,
								unsigned long long* newton_count = nullptr,
								unsigned long long* skip_count = nullptr);


#endif


