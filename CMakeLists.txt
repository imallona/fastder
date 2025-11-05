cmake_minimum_required(VERSION 4.0)
project(mls)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(cpp)

add_library(mls_lib
        cpp/Averager.cpp
        cpp/Averager.h
        cpp/BedGraphRow.cpp
        cpp/BedGraphRow.h
        # cpp/integrate_sj.cpp
        cpp/Parser.cpp
        cpp/Parser.h
        cpp/SJRow.cpp
        cpp/SJRow.h
        cpp/StitchedER.cpp
        cpp/StitchedER.h
        cpp/Integrator.cpp
        cpp/Integrator.h
        cpp/GTFRow.cpp
        cpp/GTFRow.h


)


add_executable(mls
        cpp/main.cpp
)
target_link_libraries(mls PRIVATE mls_lib)

enable_testing()
add_subdirectory(tests)