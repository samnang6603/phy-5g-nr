#!/bin/zsh

BUILD_TYPE=${1:-Debug}   # ./build.sh Release -> Release, otherwise defaults to Debug

cmake -S . -B build -DCMAKE_BUILD_TYPE=$BUILD_TYPE
cmake --build build