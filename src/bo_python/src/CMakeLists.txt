cmake_minimum_required(VERSION 3.0.0)
project(processing_input VERSION 0.1.0)

message(STATUS "Placement Refinement")

include(CTest)
enable_testing()

set(PATH_SRC ${CMAKE_CURRENT_SOURCE_DIR})

#file(GLOB_RECURSE PARSERS ${PATH_SRC}/parser/*.cpp)
#file(GLOB_RECURSE UTILS ${PATH_SRC}/*.cpp)
file(GLOB_RECURSE SRC_FILES main.cpp cell.cpp parser.cpp BBox.cpp utils.cpp drv_region.cpp)

#set(SRC_FILES ${PARSERS} ${UTILS})
#set(SRC_FILES ${PARSERS} main.cpp utils.cpp)
#add_executable(processing_input main.cpp)
#add_executable(main ${SRC_FILES} main.cpp)
add_executable(main ${SRC_FILES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-O0 -ggdb3")
set(CMAKE_C_FLAGS "-O0")
