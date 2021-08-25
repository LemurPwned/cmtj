#!/bin/bash

CORE_CODE_LOC="../../core/*.hpp" 
LAUNCH_CODE_LOC="./*.cpp"
find ${CORE_CODE_LOC} ${LAUNCH_CODE_LOC} | entr sh -c 'rm -f stack_exp && make stack && ./stack_exp && python3 plot_spectrum.py'
