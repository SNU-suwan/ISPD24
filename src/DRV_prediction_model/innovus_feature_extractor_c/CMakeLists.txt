# Specify the minimum version for CMake
cmake_minimum_required(VERSION 2.8)

# Message
message(STATUS "Input Feature Extractor")
message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS ${PROJECT_SOURCE_DIR})

# Find includes in the build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Project's name
project(featurextract)

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_STANDARD 17)

set(PATH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)



# Warning/Errors
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")

# Source code
file(GLOB_RECURSE SRC_FILES_ ${PATH_SRC}/*.cpp)

# FLUTE
file(GLOB_RECURSE FLUTE_FILES_ ${PATH_SRC}/flute-3.1/*.c)
list(REMOVE_ITEM FLUTE_FILES_ ${PATH_SRC}/flute-3.1/flute-ckt.c)
list(REMOVE_ITEM FLUTE_FILES_ ${PATH_SRC}/flute-3.1/flute-net.c)
list(REMOVE_ITEM FLUTE_FILES_ ${PATH_SRC}/flute-3.1/rand-pts.c)

set(SRC_FILES ${SRC_FILES_} ${FLUTE_FILES_})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


# Target
add_executable(main ${SRC_FILES})

# Boost Library
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(main ${Boost_LIBRARIES})


##For Debugging with GDB
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-O0")
set(CMAKE_C_FLAGS "-O0")
