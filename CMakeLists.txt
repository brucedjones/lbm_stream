cmake_minimum_required (VERSION 2.8)
project (lbm_stream)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(TWO_LATTICE "TWO_LATTICE" ON)
option(ORIGIN_STREAMING "ORIGIN_STREAMING" OFF)
option(OUTPUT "OUTPUT" OFF)

if(TWO_LATTICE)
    add_definitions(-DTWO_LATTICE)
endif()

if(ORIGIN_STREAMING)
    add_definitions(-DORIGIN_STREAMING)
endif()

if(OUTPUT)
    add_definitions(-DOUTPUT)
endif()

add_executable(lbm_stream lbm_stream.cpp)
