#!/bin/bash

CORE_CODE_LOC="../../core/*.hpp"
LAUNCH_CODE_LOC="./*.cpp"
EXEC="torque"
find ${CORE_CODE_LOC} ${LAUNCH_CODE_LOC} | entr sh -c "rm -f ${EXEC} && make ${EXEC} && ./${EXEC} && python3 plot_torque.py"
