cmake_minimum_required(VERSION 3.24)
project(sph)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(D:/openglrepo/Include include)
link_directories(D:/openglrepo/Libs)

add_executable(sph main.cpp glad.c include/config.cpp include/particle.cpp)
target_link_libraries(sph glfw3dll)