#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSPenny -g plugins/libTestGeometry.so -gd macros/geo.mac -d macros/run.mac -ro GGSRootOutput.root
