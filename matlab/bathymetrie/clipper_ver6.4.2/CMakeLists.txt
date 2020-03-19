CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0)
PROJECT(polyintersec)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../..")

ADD_EXECUTABLE(polyintersec src/clipper.cpp src/poly_intersection.cpp)
