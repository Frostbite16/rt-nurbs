#ifndef CAMERA_H
#define CAMERA_H

#include "../math/transform.h"
#include "../math/math_types.h"
#include "../math/matrix_ops.h"
#include "../math/math_utils.h"
#include <cmath>


// Define camera initial parameters configuration
// an render point is the same as an camera point
struct camera_config{
	//MathTypes::matrix4x4f render_from_camera;
	MathTypes::matrix4x4f world_from_camera;

	// Camera sub coordinates transform
	MathTypes::matrix4x4f screen_from_camera;
	MathTypes::matrix4x4f NDC_from_screen;
	MathTypes::matrix4x4f raster_from_NDC;
	MathTypes::matrix4x4f raster_from_screen;
	MathTypes::matrix4x4f screen_from_raster;
	MathTypes::matrix4x4f camera_from_raster;
	int screen_resolution_x;
	int screen_resolution_y;

	MathExtra::bounds2f screen_window;
};
struct orthographic_projection{};
struct perspective_projection{};

struct camera_gpu_payload{
	MathTypes::matrix4x4f world_from_camera;
	MathTypes::matrix4x4f camera_from_raster;
};

enum projection{orthographic, perspective};

namespace CameraOps{

	// Get orthographic camera from initial camera configuration
	inline constexpr
	camera_config getOrthographicCamera(const camera_config& initial_config, float z_near, float z_far){
		MathTypes::matrix4x4f screen_from_camera = MatrixOps::multiply(Transform::scale(1, 1, 1/(z_far - z_near)),
																	  Transform::translate(0, 0, -z_near));	
		MathTypes::matrix4x4f camera_from_raster = MatrixOps::multiply(MatrixOps::inverse(screen_from_camera), initial_config.screen_from_raster);

		return camera_config{
				initial_config.world_from_camera,
				screen_from_camera,
				initial_config.NDC_from_screen,
				initial_config.raster_from_NDC,
				initial_config.raster_from_screen,
				initial_config.screen_from_raster,
				camera_from_raster,
				initial_config.screen_resolution_x,
				initial_config.screen_resolution_y,
				initial_config.screen_window};
	}

	// Get perpective camera from initial camera configuration
	// using pbrt perpective transform matrix
	inline
	camera_config getPerspectiveCamera(const camera_config& initial_config, float fov, float n, float f){
		MathTypes::matrix4x4f p_transform = {	1,0,		0,			0,
											 	0,1,		0,			0,
												0,0, 	 f/(f-n),		1,
												0,0,	-f*n/(f-n),		0};
		float inv_tang_ang = 1.0f / tanf(MathUtils::radians(fov) / 2);
		MathTypes::matrix4x4f screen_from_camera = MatrixOps::multiply(Transform::scale(inv_tang_ang, inv_tang_ang, 1), p_transform);
		MathTypes::matrix4x4f camera_from_raster = MatrixOps::multiply(MatrixOps::inverse(screen_from_camera), initial_config.screen_from_raster);

		return camera_config{
				initial_config.world_from_camera,
				screen_from_camera,
				initial_config.NDC_from_screen,
				initial_config.raster_from_NDC,
				initial_config.raster_from_screen,
				initial_config.screen_from_raster,
				camera_from_raster,
				initial_config.screen_resolution_x,
				initial_config.screen_resolution_y,
				initial_config.screen_window};
	
	}

	// Defines camera coordinates tranform matrices
	// Except screen from raster matrix because of expecific projection necessity
	inline
	camera_config initCameraTransforms(const camera_config& initial_config){
			
		MathExtra::bounds2f screen_window(initial_config.screen_window);
		
		int screen_res_x = initial_config.screen_resolution_x;
		int screen_res_y = initial_config.screen_resolution_y;

		MathTypes::matrix4x4f NDC_from_screen = 
			MatrixOps::multiply(	Transform::scale(	1 / (screen_window.p_max_x - screen_window.p_min_x),
														1 / (screen_window.p_max_y - screen_window.p_min_y), 1),
									Transform::translate( -screen_window.p_min_x, -screen_window.p_max_y, 0));
		MathTypes::matrix4x4f raster_from_NDC = 
			Transform::scale(screen_res_x, -screen_res_y, 1);

		MathTypes::matrix4x4f raster_from_screen = MatrixOps::multiply(raster_from_NDC, NDC_from_screen);
		MathTypes::matrix4x4f screen_from_raster = MatrixOps::inverse(raster_from_screen);
	
		//MathTypes::matrix4x4f camera_from_raster = MatrixOps::multiply(MatrixOps::inverse(initial_config.screen_from_camera), screen_from_raster);

		return camera_config{
				initial_config.world_from_camera,
				initial_config.screen_from_camera,
				NDC_from_screen,
				raster_from_NDC,
				raster_from_screen,
				screen_from_raster,
				initial_config.camera_from_raster,
				initial_config.screen_resolution_x,
				initial_config.screen_resolution_y,
				initial_config.screen_window};
	}

	// Initiate the camera
	inline
	camera_config init_camera(const MathTypes::matrix4x4f world_from_camera, int res_x, int res_y){
		camera_config clean_config = {0};
			
		clean_config.world_from_camera = world_from_camera;	
		clean_config.screen_resolution_x = res_x;
		clean_config.screen_resolution_y = res_y;	

		float aspect = float(res_x) / float(res_y);
		if(aspect > 1.0f)
			clean_config.screen_window = {-aspect, -1.0f, aspect, 1.0f}; 	
		else
			clean_config.screen_window = {-1.0f, -1.0f/aspect, 1.0f, 1.0f/aspect}; 	
		clean_config = initCameraTransforms(clean_config);


		return clean_config; 
	
	}

	// Using specialized templates for static conditional
	template<typename T>
	CORE_D MathExtra::ray calculateRay(const camera_gpu_payload& camera, int pixel_x, int pixel_y);

	
	// For ortographic camera
	// returns ray with origin (pixel_x, pixel_y) and direction perpenticular to the screen in world coordinates
	template <>
	inline CORE_D MathExtra::ray calculateRay<orthographic_projection>(const camera_gpu_payload& camera, int pixel_x, int pixel_y){
		// I'm using casting from cuda library here
		MathTypes::point4f point_in_camera = MatrixOps::multiply(camera.camera_from_raster, 
																 MathTypes::point4f{float(pixel_x) + 0.5f, float(pixel_y) + 0.5f, 0.0f, 1.0f});
		MathTypes::point4f ray_origin = MatrixOps::multiply(camera.world_from_camera, point_in_camera);
		MathTypes::vector3f ray_direction = MathTypes::vector3f{ camera.world_from_camera.m[8], camera.world_from_camera.m[9], camera.world_from_camera.m[10]};	

		return MathExtra::ray{MathTypes::vector3f{ray_origin.v[0], ray_origin.v[1], ray_origin.v[2]}, ray_direction};
	}

	// For perspective camera
	template <>
	inline CORE_D MathExtra::ray calculateRay<perspective_projection>(const camera_gpu_payload &camera, int pixel_x, int pixel_y){
		MathTypes::point4f point_in_camera = MatrixOps::multiply(camera.camera_from_raster, 
																 MathTypes::point4f{float(pixel_x) + 0.5f, float(pixel_y) + 0.5f, 0.0f, 1.0f});

		MathTypes::vector3f p_camera = VectorOps::normalize(MathTypes::vector3f{point_in_camera.v[0], point_in_camera.v[1], point_in_camera.v[2]});

		MathTypes::point3f ray_origin = MathTypes::point3f{camera.world_from_camera.m[12], camera.world_from_camera.m[13], camera.world_from_camera.m[14]};
		MathTypes::vector4f ray_direction = MatrixOps::multiply(camera.world_from_camera, MathTypes::vector4f{p_camera.v[0], p_camera.v[1], p_camera.v[2], 0.0f} );

		return MathExtra::ray{MathTypes::vector3f{ray_origin.v[0], ray_origin.v[1], ray_origin.v[2]}, 
							  MathTypes::vector3f{ray_direction.v[0], ray_direction.v[1], ray_direction.v[2]}};
	}

}

#endif
