#!/bin/bash

CORE_CODE_LOC="../../core/*.hpp" 
LAUNCH_CODE_LOC="./*.cpp ./*.py"
find ${CORE_CODE_LOC} ${LAUNCH_CODE_LOC} | entr sh -c 'rm -f pimm && make pimm && ./pimm && python3 pimm.py PIMM.csv'