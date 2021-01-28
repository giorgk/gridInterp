cmake_minimum_required(VERSION 3.17)
project(gridInterp)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(gridInterp
        gridInterp.cpp
        gridInterp.h
        tests.h)
