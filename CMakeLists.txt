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

# packaging
set(CPACK_PACKAGE_NAME "fastder")
set(CPACK_PACKAGE_VENDOR "Martina Lehmann")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "fastder, a fast expressed regions finder")
set(CPACK_PACKAGE_VERSION "0.0.1")
set(CPACK_PACKAGE_CONTACT "martina.lavanya@gmail.com")

# Debian-specific packaging settings
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Martina Lehmann")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6")

include(CPack)
