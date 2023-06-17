cmake_minimum_required(VERSION 3.25)
project(untitled55 C)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_C_STANDARD 11)

add_executable(untitled55 main.c matrix.c)
target_link_libraries(untitled55 ${MPI_C_LIBRARIES})