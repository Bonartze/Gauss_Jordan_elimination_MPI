cmake_minimum_required(VERSION 3.25)
project(untitled55 C)

find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER "${MPI_C_COMPILER}")
add_executable(untitled55 main.c matrix.c)
target_link_libraries(untitled55 PUBLIC MPI::MPI_C)
