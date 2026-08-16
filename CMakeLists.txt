cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

project(
	rt_nurbs 
	VERSION 0.0.1
	LANGUAGES CXX 
	CUDA
)
set(CMAKE_CXX_STANDARD 17)
set(FETCHCONTENT_QUIET FALSE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

find_package(CUDAToolkit REQUIRED)

add_executable(${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PRIVATE CUDA::cudart)

target_sources(${PROJECT_NAME}
	PRIVATE
		src/core/main.cpp
		src/core/camera.h
		src/core/core_defines.h
		src/math/transform.h
		src/math/bounds_ops.h
		src/math/math_types.h
		src/math/math_utils.h
		src/math/matrix_ops.h
		src/math/vector_ops.h
		src/gpu/soa.h
		src/gpu/memory.h
		src/gpu/generateRaysKernel.h
		src/gpu/generateRaysKernel.cu
		src/gpu/intersect.h
		src/gpu/intersect.cu
		src/gpu/testShading.h
		src/gpu/testShading.cu
		src/constants.h
		src/math/nurbs.h
		src/math/bezier_patch.h
)

set_target_properties( ${PROJECT_NAME}
						PROPERTIES
							CUDA_SEPARABLE_COMPILATION ON
							CUDA_ARCHITECTURES native)

# ---- throwaway data-gathering target (benchmark.cpp instead of main.cpp) ----
add_executable(benchmark)

target_link_libraries(benchmark PRIVATE CUDA::cudart)

target_sources(benchmark
	PRIVATE
		src/core/benchmark.cpp
		src/core/camera.h
		src/core/core_defines.h
		src/math/transform.h
		src/math/bounds_ops.h
		src/math/math_types.h
		src/math/math_utils.h
		src/math/matrix_ops.h
		src/math/vector_ops.h
		src/gpu/soa.h
		src/gpu/memory.h
		src/gpu/generateRaysKernel.h
		src/gpu/generateRaysKernel.cu
		src/gpu/intersect.h
		src/gpu/intersect.cu
		src/constants.h
		src/math/nurbs.h
		src/math/bezier_patch.h
)

set_target_properties( benchmark
						PROPERTIES
							CUDA_SEPARABLE_COMPILATION ON
							CUDA_ARCHITECTURES native)



