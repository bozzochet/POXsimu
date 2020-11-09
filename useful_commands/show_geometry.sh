#!/bin/bash

rm -rf plugins/libTestGeometry.gdml
GGSWolowitz -g plugins/libTestGeometry.so -t gdml -o plugins/libTestGeometry.gdml
GGSLeonard -g plugins/libTestGeometry.gdml

