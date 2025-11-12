project(mls-tests)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/heads/main.zip
)
FetchContent_MakeAvailable(googletest)

add_executable(mls-tests
        UnitTests.cpp
)

target_link_libraries(mls-tests
        PRIVATE
        mls_lib
        GTest::gtest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(mls-tests)
