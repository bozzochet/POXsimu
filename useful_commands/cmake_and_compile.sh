#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

cmake -DCMAKE_INSTALL_PREFIX=${PWD/%_build}_install ${PWD/%_build}

make
#make install

rm -rf plugins/libTestGeometry.gdml
GGSWolowitz -g plugins/libTestGeometry.so -o plugins/libTestGeometry.gdml -t gdml

