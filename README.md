# POXsimu

GGS Simulation of POX prototype detector

vis.mac:  
-> gamma 100 MeV   
-> 10000 events  
geo.mac:  
-> 10+4 layers  
-> thickness 300 um  
  
```
GGSPenny -g plugins/libTestGeometry.so -gd geo.mac -d vis.mac -ro GGSRootOutput.root > GGSOut.txt
```
  
produces the file GGSRootOutput.root, and the .wrl files of the first 100 events see for example g4_03.wrl for conferted gamma  
(the creation of the .wrl event display is suppressed commenting the
```
/vis/open VRML2FILE
```
line the vis.mac macro)

Then, to analyze the GGS output:

root [0] .L Analysis/Analysis.C 
root [1] SimpleAnalysis("GGSRootOutput.root", "anaOut.root")