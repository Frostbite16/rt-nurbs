#ifndef SOA_H
#define SOA_H

namespace SoAQueue{
	
	struct ray_queue{
		float* origin_x;
		float* origin_y;
		float* origin_z;
		float* direction_x;
		float* direction_y;
		float* direction_z;
		int* idx_interact;
	};

	struct ray_pixels{
		

	};


	struct frame_buffer{
		float* r;
		float* g;
		float* b;
	};
}

namespace SoAPrimitve{
	
	struct sphere{
		float* world_pos_x;
		float* world_pos_y;
		float* world_pos_z;
		float* radius;
		float* z_min;
		float* z_max;
	};

}

#endif
