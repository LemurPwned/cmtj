#!/bin/bash

# builds the test suite
rm -rf ./build
cmake -S . -B build
cmake --build build
cd build && ctest --rerun-failed --output-on-failure
