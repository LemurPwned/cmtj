#!/bin/bash
find {../core,.}/*{.cpp,.hpp} CMakeLists.txt  | entr sh -c "./build_script.sh"
