# POXsimu

GGS Simulation of POX prototype detector

Ingredients:

- macros/vis.mac: datacard, a la Geant, to set the simulation parameters 
- macros/geo.mac: datacard for the parametric geometry
- {src,include}/DetectorConstruction.{cc,hh}: definition, a la Geant, of the geometry 
- Analysis/Analysis.C: ROOT macro to read the GGS output file

Typical commands:

- simuation

```
GGSPenny -g plugins/libTestGeometry.so -gd geo.mac -d vis.mac -ro GGSRootOutput.root > GGSOut.txt
```
  
produces the file GGSRootOutput.root, and the .wrl files of the first 100 events see for example g4_03.wrl for conferted gamma  
(the creation of the .wrl event display is suppressed commenting the
```
/vis/open VRML2FILE
```
line in the vis.mac macro)

- conversion from GGS output to plain ROOT file
Then, to analyze the GGS output:
```
root [0] .L Analysis/Analysis.C 
root [1] SimpleAnalysis("GGSRootOutput.root", "anaOut.root")
```

