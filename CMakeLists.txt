cmake_minimum_required(VERSION 3.19)
project(RayTracing)

set(CMAKE_CXX_STANDARD 17)

add_executable(main main.cpp Scene.cpp Scene.h structs.h)