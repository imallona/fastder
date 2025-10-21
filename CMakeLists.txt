cmake_minimum_required(VERSION 4.0)
project(mls)

set(CMAKE_CXX_STANDARD 20)

include_directories(cpp)

add_executable(mls
        cpp/Averager.cpp
        cpp/Averager.h
        cpp/BedGraphRow.cpp
        cpp/BedGraphRow.h
        # cpp/integrate_sj.cpp
        cpp/main.cpp
        cpp/Parser.cpp
        cpp/Parser.h
        cpp/SJRow.cpp
        cpp/SJRow.h
        cpp/StitchedER.cpp
        cpp/StitchedER.h
        cpp/Integrator.cpp
        cpp/Integrator.h
        #cpp/test.cpp
)
