cmake_minimum_required(VERSION 4.0)
project(mls)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(.)

add_executable(cpp
        test.cpp)
