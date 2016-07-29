cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_FLAGS "-g3")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wall")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

MESSAGE(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(
  /usr/local/include
  /usr/include
)

ADD_EXECUTABLE(scpxf main.cpp)
