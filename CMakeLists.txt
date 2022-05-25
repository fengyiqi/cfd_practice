cmake_minimum_required(VERSION 3.16)

project(CFD_Practice)

set(CMAKE_CXX_STANDARD 20)

include_directories(src)
file(GLOB_RECURSE SOURCE_FILES "src/*cpp")
add_executable(CFD_Practice ${SOURCE_FILES})