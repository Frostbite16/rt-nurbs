#ifndef BVH_H
#define BVH_H

#include "../math/bezier_patch.h"
#include "../math/math_types.h"
#include "../math/bounds_ops.h"
#include <cmath>
#include <sys/types.h>

// Left first can have two definitions depending if the node is a leaf or not
// if it is a leaf it is the index of the primitive
// if it is not a leaf it is the index of the left child node
struct BVHNode{
	MathExtra::bounds3f node_bound;
	uint left_first;
	uint prim_count;
};

// Variables needed for the building process of the BVH
struct BVHBuild{
	MathExtra::bounds3f* bounds;
	uint* pri_idx;
	int primitive_size;
	int nodes_used;
};


namespace BvhOps{
	
	// Free the dynamic allocated arrays of build
	inline
	void freeBuild(BVHBuild& build){
		delete[] build.pri_idx;
		delete[] build.bounds;
	}
	
	// Allocate the arrays of build
	inline 
	BVHBuild initBuild(const BVHBuild& build, int primitive_size){
		MathExtra::bounds3f* pri_bounds = new MathExtra::bounds3f[primitive_size];
		uint* pri_idx = new uint[primitive_size];
		return BVHBuild{pri_bounds, pri_idx, 0, 0};
	}

	// Evaluate the cost of spliting the node in a certain axis_idx
	// Implementation from: https://jacco.ompf2.com/2022/04/18/how-to-build-a-bvh-part-2-faster-rays/
	inline
	float evaluateSAH( const BVHBuild& build, const BVHNode& node, int axis, float pos){
		MathExtra::bounds3f left_bound = {	INFINITY, INFINITY, INFINITY,
											-INFINITY, -INFINITY, -INFINITY};
		MathExtra::bounds3f right_bound = {	INFINITY, INFINITY, INFINITY,
											-INFINITY, -INFINITY, -INFINITY};
		int left_count = 0, right_count = 0;
		for(int i=0; i < node.prim_count; i++){
			MathExtra::bounds3f& prim_bound = build.bounds[build.pri_idx[node.left_first + i]]; 
			MathTypes::point3f bd_centroid = BoundsOps::centroid(prim_bound);
			if(bd_centroid.v[axis] < pos){
				left_count++;
				left_bound = BoundsOps::boundUnion(left_bound, prim_bound);
				
			}
			else{
				right_count++;
				right_bound =BoundsOps::boundUnion(right_bound, prim_bound);
			}
		}
		float cost = left_count * BoundsOps::bound_area(left_bound) + right_count * BoundsOps::bound_area(right_bound);
		return cost > 0 ? cost : 1e30f;
	}
	
	// Resize the bounds of the node to fit all their child bounds
	inline
	void updateNodeBounds(const BVHBuild& build, BVHNode* nodes, uint node_idx){
		BVHNode& cur_node = nodes[node_idx];
		cur_node.node_bound = MathExtra::bounds3f{	 INFINITY, 	INFINITY,  INFINITY, 
													-INFINITY, -INFINITY, -INFINITY};

		for(int i = cur_node.left_first; i<cur_node.left_first+cur_node.prim_count; i++){
			cur_node.node_bound = BoundsOps::boundUnion(cur_node.node_bound,
													 	build.bounds[build.pri_idx[i]]); 
		}
	}
	
	// Divide an bound between left and right child
	inline
	void subdivide(BVHBuild& build, BVHNode* nodes, uint node_idx){
		BVHNode& cur_node = nodes[node_idx];
		
		int best_axis = 1;
		float best_pos = 0, best_cost = 1e30f;
		
		for(int axis=0; axis < 3; axis++){
			for(int i=0; i<cur_node.prim_count; i++){
				MathTypes::point3f bd_centroid = BoundsOps::centroid(build.bounds[build.pri_idx[cur_node.left_first + i]]);
				float candidate_pos = bd_centroid.v[axis];
				float cost = evaluateSAH(build, cur_node, axis, candidate_pos);
				if(cost < best_cost){
					best_pos = candidate_pos; best_axis = axis, best_cost = cost;
				}
			}
		}

		int axis_idx = best_axis;
		float split_pos = best_pos;

		// Organizes primitives on the node
		// putting the primitives left of the axis position to the left of the node
		int i = cur_node.left_first;
		int j = i + cur_node.prim_count - 1;

		while(i <= j){
			MathTypes::point3f bd_centroid = BoundsOps::centroid(build.bounds[build.pri_idx[i]]);
			if(bd_centroid.v[axis_idx] < split_pos){
				i++;
			}
			else{
				int buffer = build.pri_idx[i];
				build.pri_idx[i] = build.pri_idx[j];
				build.pri_idx[j] = buffer;
				j--;
			} // Implement a swap
		}
		
		int left_count = i - cur_node.left_first;
		if(left_count == 0 || left_count == cur_node.prim_count) return;

		int left_child_idx = build.nodes_used++;
		int right_child_idx = build.nodes_used++;
		
		nodes[left_child_idx].left_first = cur_node.left_first;
		nodes[left_child_idx].prim_count = left_count;

		nodes[right_child_idx].left_first = i;
		nodes[right_child_idx].prim_count = cur_node.prim_count - left_count;
		cur_node.left_first = left_child_idx;
		cur_node.prim_count = 0;

		updateNodeBounds(build, nodes, left_child_idx);
		updateNodeBounds(build, nodes, right_child_idx);

		subdivide(build, nodes, left_child_idx);
		subdivide(build, nodes, right_child_idx);
	}



	// Build root node of the bvh tree
	inline
	void buildBVH(BVHNode* nodes, uint node_root_idx, bezierPatch*& patches, int primitive_size){
		BVHBuild build = initBuild(build, primitive_size);

		BVHNode& root = nodes[node_root_idx];
		root.left_first = 0;
		root.prim_count = primitive_size;
		root.node_bound = MathExtra::bounds3f{	 INFINITY,  INFINITY,  INFINITY, 
												-INFINITY, -INFINITY, -INFINITY};

		build.primitive_size = primitive_size;
		build.nodes_used = 1;
		for(int i=0; i<primitive_size; i++){
			build.bounds[i] = BezierPatchOps::bounds(patches[i]);
			build.pri_idx[i] = i;
			root.node_bound = BoundsOps::boundUnion(root.node_bound, build.bounds[i]);
		}
		subdivide(build, nodes, node_root_idx);
		
		// reordering spheres based on the pri_idx array
		bezierPatch* primitives_buffer = new bezierPatch[primitive_size];
		for(int i=0; i<build.primitive_size; i++){
			primitives_buffer[i] = patches[build.pri_idx[i]];
		}
	
		for(int i=0; i<build.primitive_size; i++){
			patches[i] = primitives_buffer[i];
		}

		delete[] primitives_buffer;
		freeBuild(build);
	}
	
	

	

	
}

#endif

