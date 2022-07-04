cmake_minimum_required(VERSION 3.21)
project(matrix_inversion)

set(CMAKE_CXX_STANDARD 17)

add_executable(matrix_inversion main.cpp matrix_inversion.cpp matrix_inversion.hpp)
