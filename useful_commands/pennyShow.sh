#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

mv `pwd`/macros/vis.mac `pwd`
mv `pwd`/macros/rootlogon.C `pwd`

GGSPenny -g plugins/libTestGeometry.so -gd macros/geo.mac -d macros/run.mac -X
