# POXsimu

GGS Simulation of POX prototype detector

## Simulation:

### Main ingredients:

- macros/run.mac: datacard, a la Geant, to set the simulation parameters 
- macros/geo.mac: datacard for the parametric geometry
- {src,include}/DetectorConstruction.{cc,hh}: definition, a la Geant, of the geometry 
- Analysis/Analysis.C: ROOT macro to read the GGS output file

### Typical commands:

- compilation and installation:

```
cd <build_path>
cmake -DCMAKE_INSTALL_PREFIX=<installation_path> <source_path>
make
make install
```

- simulation:

```
GGSPenny -g plugins/libTestGeometry.so -gd macros/geo.mac -d macros/run.mac -ro GGSRootOutput.root
```

produces the file `GGSRootOutput.root` that can be later analyzed (for example with the `Analysis.C` macro, see below) or inspected with `GGSLeonard` (see below).

- conversion from parametric geometry to GDML

```
GGSWolowitz -g plugins/libTestGeometry.so -gd macros/geo.mac -t gdml -o plugins/libTestGeometry.gdml
```

- conversion from parametric geometry to VGM (http://ivana.home.cern.ch/ivana/VGM.html)

```
GGSWolowitz -g plugins/libTestGeometry.so -gd macros/geo.mac -t vgm -o plugins/libTestGeometry.vgm.root
```

- open the geometry with `GGSLeonard`:

```
GGSLeonard -g plugins/libTestGeometry.vgm.root
```

- display the events with `GGSLeonard`:

```
GGSLeonard -g plugins/libTestGeometry.vgm.root -i GGSRootOutput.root
```

- run `GGSPenny` interactively and look the event visualization (it requires Geant4 to be compiled with Qt and OpenGL support):

```
GGSPenny -g plugins/libTestGeometry.so -gd macros/geo.mac -d macros/run.mac -X
```
(adding the `-X` flag). If the `macros/run.mac` contains the usual  `/run/beamOn NNNN` the full simulation is performed and only after it the visualizazion is opened. Commenting that line just open the visualization. To simulate one event is enough, in the visualizer window, to click on the upper right green arrow or to send a `/run/beamOn 1` in the "Session:" input form.

The look of the visualization is customized by the `/macros/vic.mac` file from the "source". This file in the "install" is moved to the main directory since `GGSPenny` is looking in the current working dir from which it is executed (usually the main directory of the "install").

- during the `GGSPenny` run, one can save each event as `.wrl` file. This needs the inclusion of some `/vis/XXXX` instructions. An axample of this is done inside the `macros/run_savewrl.mac` file:

```
GGSPenny -g plugins/libTestGeometry.so -gd macros/geo.mac -d macros/run_savewrl.mac -ro GGSRootOutput.root
```

## Data Analysis

```
root [0] .L Analysis/Analysis.C 
root [1] SimpleAnalysis("GGSRootOutput.root", "anaOut.root")
```

