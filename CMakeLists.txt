cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

project(
	rt_nurbs 
	VERSION 0.0.1
	LANGUAGES CXX 
	CUDA
)
set(CMAKE_C_STANDARD 99)
set(FETCHCONTENT_QUIET FALSE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

FetchContent_Declare(
	raylib
	GIT_REPOSITORY "https://github.com/raysan5/raylib.git"
    GIT_TAG "master"
	GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(raylib)


add_executable(${PROJECT_NAME})


target_sources(${PROJECT_NAME}
	PRIVATE
		src/main.cu
)

target_link_libraries(${PROJECT_NAME} PRIVATE raylib)



