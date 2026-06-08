cmake_minimum_required(VERSION 3.16)
project(openlb_droplet_particle_jump LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(OPENLB_DIR "" CACHE PATH "Path to an OpenLB source tree")
if(NOT OPENLB_DIR)
  message(FATAL_ERROR "Set -DOPENLB_DIR=/path/to/openlb")
endif()

add_executable(droplet_particle_jump3d src/droplet_particle_jump3d.cpp)
target_include_directories(droplet_particle_jump3d PRIVATE ${OPENLB_DIR} ${OPENLB_DIR}/src)
