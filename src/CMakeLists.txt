cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Release)

#including values
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(vertex_cover)

# add_executable(
#     vertex_cover
#     main.cpp
#     WeightedVertexGraph.cpp
#     utilities.cpp
#     heuristics.cpp
#     LocalSearch.cpp
#     IteratedLocalSearch.cpp
# )

set(EXEC_MAIN "main.cpp")
set(EXPERIMENT_MAIN "experiments.cpp")
set(LIB_SRCS 
    WeightedVertexGraph.cpp
    utilities.cpp
    heuristics.cpp
    LocalSearch.cpp
    IteratedLocalSearch.cpp
)

add_library( mysharedlib SHARED ${LIB_SRCS} )

add_executable( vertex_cover ${EXEC_MAIN} )
target_link_libraries( vertex_cover mysharedlib )

add_executable( experiments ${EXPERIMENT_MAIN})
target_link_libraries( experiments mysharedlib )