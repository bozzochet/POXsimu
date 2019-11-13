# POXsimu

GGS Simulation of POX prototype detector\n

vis.mac:\n
-> gamma 100 MeV\n
-> 10000 events\n
geo.mac:\n
-> 10+4 layers\n
-> thickness 300 um\n
\n
```
GGSPenny -g plugins/libTestGeometry.so -gd geo.mac -d vis.mac -ro GGSRootOutput.root > GGSOut.txt
```
\n
produces the file GGSRootOutput.root, and the .wrl files of the first 100 events\n
see for example g4_03.wrl for conferted gamma\n
(the creation of the .wrl event display is suppressed commenting the\n
/vis/open VRML2FILE\n
line the vis.mac macro)\n 