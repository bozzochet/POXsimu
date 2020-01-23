# POXsimu

GGS Simulation of POX prototype detector

Ingredients:

- macros/vis.mac: datacard, a la Geant, to set the simulation parameters 
- macros/geo.mac: datacard for the parametric geometry
- {src,include}/DetectorConstruction.{cc,hh}: definition, a la Geant, of the geometry 
- Analysis/Analysis.C: ROOT macro to read the GGS output file

Typical commands:

- compilation and installation:

```
cd <build_path>
cmake -DCMAKE_INSTALL_PREFIX=<installation_path> <source_path>
make
make install
```

- simulation:

```
GGSPenny -g plugins/libTestGeometry.so -gd macros/geo.mac -d macros/vis.mac -ro GGSRootOutput.root > GGSOut.txt
```
  
produces the file GGSRootOutput.root, and the .wrl files of the first 100 events, see for example g4_03.wrl for converted gamma  
(the creation of the .wrl event display is suppressed commenting the
```
/vis/open VRML2FILE
```
line in the vis.mac macro)

- conversion from parametric geometry to GDML
```
GGSWolowitz -g plugins/libTestGeometry.so -gd macros/geo.mac -t gdml -o plugins/libTestGeometry.gdml
```

- conversion from parametric geometry to VGM (http://ivana.home.cern.ch/ivana/VGM.html)
```
GGSWolowitz -g plugins/libTestGeometry.so -gd macros/geo.mac -t vgm -o plugins/libTestGeometry.vgm.root
```

- opening of the geometry display:
```
GGSLeonard -g plugins/libTestGeometry.vgm.root
```

- opening of the event display:
```
GGSLeonard -g plugins/libTestGeometry.vgm.root -i GGSRootOutput.root
```

- conversion from GGS output to plain ROOT file:
```
root [0] .L Analysis/Analysis.C 
root [1] SimpleAnalysis("GGSRootOutput.root", "anaOut.root")
```

