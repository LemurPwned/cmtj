#!/bin/bash

CORE_CODE_LOC="../../core/*.hpp" 
LAUNCH_CODE_LOC="./*.cpp"
find ${CORE_CODE_LOC} ${LAUNCH_CODE_LOC} | entr sh -c 'rm -f vsd && make vsd && ./vsd && python3 vsd.py VSD.csv'