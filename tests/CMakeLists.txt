project(fastder-tests)

include(FetchContent)

# disable installation of googletest/gmock before fetching
set(INSTALL_GTEST OFF CACHE BOOL "Disable installation of googletest" FORCE)
set(INSTALL_GMOCK OFF CACHE BOOL "Disable installation of googlemock" FORCE)
set(gtest_install OFF CACHE BOOL "Disable installation of googletest (legacy var)" FORCE)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/heads/main.zip
)
FetchContent_MakeAvailable(googletest)

add_executable(fastder-tests
    UnitTests.cpp
)

target_link_libraries(fastder-tests
        PRIVATE
        fastder_lib
        GTest::gtest
        GTest::gtest_main
)

# When fastder is built with libBigWig, the BigWig round-trip and equivalence
# tests in UnitTests.cpp include <bigWig.h> and call libBigWig's writer API,
# so the test target needs both the include path and the static lib too.
if(FASTDER_USE_LIBBIGWIG)
    target_link_libraries(fastder-tests PRIVATE BigWig)
    target_include_directories(fastder-tests PRIVATE ${libbigwig_src_SOURCE_DIR})
    target_compile_definitions(fastder-tests PRIVATE FASTDER_USE_LIBBIGWIG)
endif()

include(GoogleTest)
gtest_discover_tests(fastder-tests)
