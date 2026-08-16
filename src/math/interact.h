#ifndef INTERACT_H
#define INTERACT_H

#include "../math/math_types.h"

struct surfaceInteraction{
	
	float t_hit;
	MathTypes::vector3f wo;
	MathTypes::point2f uv;
	MathTypes::vector3f dpdu;
	MathTypes::vector3f dpdv;
	MathTypes::vector3f dndu;
	MathTypes::vector3f dndv;

	// Removed for now, but can be used later, we will need this if we add mesh
	// it is used for texture mapping 
	//int face_index;
	

};






#endif
