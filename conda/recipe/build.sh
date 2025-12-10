#!/bin/bash

set -ex

mkdir -p build
cd build
cmake ${CMAKE_ARGS} -DCONDA_BUILD=ON ..
make -j${CPU_COUNT}
