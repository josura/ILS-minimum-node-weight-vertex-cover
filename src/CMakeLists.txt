cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Debug)

#including values
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(vertex_cover)

add_executable(
    vertex_cover
    main.cpp
    WeightedVertexGraph.cpp
    utilities.cpp
    heuristics.cpp
    LocalSearch.cpp
    #IteratedLocalSearch.cpp
)