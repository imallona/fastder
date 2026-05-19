option(CONDA_BUILD "Building inside conda-build" OFF)
# Build with libBigWig support so Parser::read_bigwig can parse .bw files
# directly. Off by default to keep the build hermetic; turning it on requires
# libBigWig and its transitive deps (libcurl, zlib) to be available.
option(FASTDER_USE_LIBBIGWIG "Build with libBigWig support for direct .bw reading" OFF)

cmake_minimum_required(VERSION 4.0)
project(fastder)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)

include_directories(cpp)

add_library(fastder_lib
        cpp/Averager.cpp
        cpp/BedGraphRow.cpp
        cpp/Parser.cpp
        cpp/SJRow.cpp
        cpp/StitchedER.cpp
        cpp/Integrator.cpp
        cpp/GTFRow.cpp
)


add_executable(fastder
        cpp/main.cpp
)
target_link_libraries(fastder PRIVATE fastder_lib Threads::Threads)

if(FASTDER_USE_LIBBIGWIG)
    # libBigWig ships its own CMakeLists and exports a static target called
    # BigWig along with the headers under its source directory. Pull it in,
    # link against the target, and add the source dir to our include path so
    # the <bigWig.h> include in Parser.cpp resolves.
    include(FetchContent)
    FetchContent_Declare(
        libBigWig_src
        GIT_REPOSITORY https://github.com/dpryan79/libBigWig.git
        GIT_TAG 0.4.7
    )
    FetchContent_MakeAvailable(libBigWig_src)

    target_link_libraries(fastder_lib PRIVATE BigWig)
    target_include_directories(fastder_lib PRIVATE ${libbigwig_src_SOURCE_DIR})
    target_compile_definitions(fastder_lib PRIVATE FASTDER_USE_LIBBIGWIG)
endif()

enable_testing()

if(NOT CONDA_BUILD)
       add_subdirectory(tests)
endif()

# packaging
set(CPACK_PACKAGE_NAME "fastder")
set(CPACK_PACKAGE_VENDOR "Martina Lehmann")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "fastder: A Fast, Agnostic, Splice-Junction-Based Tool for the Detection of Expressed Regions")
set(CPACK_PACKAGE_VERSION "0.1.0")
set(CPACK_PACKAGE_CONTACT "martina.lavanya@gmail.com, izaskun.mallona.work@gmail.com")

# Debian-specific packaging settings
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Martina Lehmann")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6")

include(CPack)
