/*
 * DetectorConstruction.h
 *
 *  Created on: 22 Feb 2017
 *      Authors: Junjing Wang, Ming Xu, Zheng Quan
 *      Adapted by: Nicola Mori
 */

#ifndef DETECTORCONSTRUCTION_HH_
#define DETECTORCONSTRUCTION_HH_

// GGS headers
#include "geometry/GGSVGeometryConstruction.h"

// Geant4 headers
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"

//#include "DetectorGeometry.hh"
//#include "DetectorMessenger.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

class G4GenericMessenger;
//class DetectorMessenger;

//class G4PVReplica;
//class G4VPVParameterisation;

class DetectorConstruction: public GGSVGeometryConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();
  
  bool ExportParameters();
  const std::string GetVersion();
  
public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
public:
  void updateGeometry();
  
  /*! @brief Return the pointer to physical world */
  G4VPhysicalVolume* GetVolume() {
    return fPhysicalWorld;
  }
  
private:
  void DefineMaterials();
  
  G4bool fCheckOverlaps;
  
  G4VPhysicalVolume* fPhysicalWorld;
  
  G4UniformMagField* fMagField;
  G4LogicalVolume* magnetLog;
  
  //  DetectorMessenger* detMessenger;
  G4GenericMessenger *_messenger;
  
  
  /// GEOMETRY DEFAULT VALUES	
  
  
  G4String fAlign;
  
  
  G4double fLayersOffsetX = 2. * cm;
  G4double fLayersOffsetY = 3.52 * cm;
  
  /// target tracker
  // target tracker z pos
  G4double fTargetPosZ = 0. * cm;
  
  // target tracker box
  
  G4double fTargetSizeX = 50.  * cm;
  G4double fTargetSizeY = 50.  * cm;
  G4double fTargetSizeZ = 50.  * cm;
  
  // offset position of first layer	
  G4double fTargetLayerOffsetZ = 0. * cm;	
  // number of target layers	
  G4int fTargetLayerNo = 10;
  // MAYBE WITH A TBITSARRAY
  // target alignment (X=0 Y=1)
  // G4int tAlign[targetLayerNo]={0,0,0,0,0};

  // number of mini ("new") layers
  G4int fTargetLayerMiniNo = 4;
  // number of short layers
  G4int fTargetLayerShortNo = 1;
    // number of long layers
  G4int fTargetLayerLongNo = 0;
  // number of Dampe layers
  G4int fTargetLayerDampeNo = 6;
  // inter layer distance
  G4double fTargetLayerDistance = 3. * cm;

  // AMS like sensors dimensions
  G4double fAMSTileX = 4.  * cm;
  G4double fAMSTileY = 7.04 * cm;
  G4double fAMSTileThickness = 0.300 * mm;
  // sensor SPitch	
   G4double fAMSTileSPitch = 0.110 * mm;
  // sensor KPitch
  G4double fAMSTileKPitch = 0.208 * mm;

  // DAMPE like sensors dimensions
  G4double fDAMPETileX = 9.5  * cm;
  G4double fDAMPETileY = 9.5 * cm;
  G4double fDAMPETileThickness = 0.320 * mm;
  // sensor SPitch	
   G4double fDAMPETileSPitch = 0.240 * mm;
  // sensor KPitch
  G4double fDAMPETileKPitch = 9.5 * cm;

  // Mini ("new") like sensors dimensions
  G4double fTileX = 9.5  * cm;
  G4double fTileY = 9.5 * cm;
  G4double fTileThickness = 0.150 * mm;
  // sensor SPitch	
   G4double fTileSPitch = 0.150 * mm;
  // sensor KPitch
  G4double fTileKPitch = 9.5 * cm;
  
  // number of sensors per ladder
  G4int fLayerMiniTileNo = 1;
  G4int fLayerShortTileNo = 4; 
  G4int fLayerLongTileNo = 12;
  G4int fLayerDampeTileNo = 4; 
	
  /// spectrometer and magnet
  
  // spectrometer Z offset
  G4double fSMOffsetZ = 0. * cm;
  // spectrometer box	
  G4double fSMSizeX = 50.  * cm;
  G4double fSMSizeY = 50.  * cm;
  G4double fSMSizeZ = 50.  * cm;
  
  // magnetic vol dimensions
  G4double fMagVolR = 7. * cm;
  G4double fMagVolZ = 20. * cm;
  G4double fMagFieldVal = 0.05 *tesla;
  
  // number of spectrometer layers	
  G4int fSMLayerNo = 4;
  // FRONT inter layer distance
  G4double fLayerFDistance = 10. * cm;
  // END inter layer distance
  G4double fLayerEDistance = 10. * cm;
	
  // FRONT layer-magfield gap
  G4double fLayerFGap = 0. * cm; // distance from magnetic vol
  // END layer-magfield gap
  G4double fLayerEGap = 0. * cm; // distance from magnetic vol

    
  //  G4VPVParameterization * siLayerParam;
    /*
  G4double pad_x;
  G4double pad_y;
  G4double pad_z; //thickness
  G4double zfirst;
  G4double paddist;
  G4Int npads;
  */


};

#endif /* DETECTORCONSTRUCTION_HH_ */
