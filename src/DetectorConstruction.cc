/*
 *  DetectorConstruction.cc
 *
 *  Created on: 25 Apr 2018
 *  Author: Viviana Scherini
*/


#include "DetectorConstruction.hh"
// parameterisation
//#include "SiLayerParameterisation.hh"

// GGS headers
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

#include "G4GenericMessenger.hh"

// Geant4 headers
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
//#include "G4VPVParameterisation.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4UnionSolid.hh"

#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
//#include "G4Mag_UsualEqRhs.hh"
//#include "G4MagIntegratorStepper.hh"
//#include "G4ChordFinder.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make DetectorConstruction a GGS plugin
GeometryPlugin (DetectorConstruction)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :
GGSVGeometryConstruction(), fCheckOverlaps(false), fPhysicalWorld(NULL) {
  DefineMaterials();
  
  _messenger = new G4GenericMessenger(this, "/Detector/");
  _messenger->DeclareProperty("layerAlignment", fAlign, "Set alignment of layers (X=0, Y=1)");
  _messenger->DeclareProperty("layersOffsetX", fLayersOffsetX, "Set layers X offset").SetUnit("cm");
  _messenger->DeclareProperty("layersOffsetY", fLayersOffsetY, "Set layers Y offset").SetUnit("cm");
  _messenger->DeclareProperty("targetPosZ", fTargetPosZ, "Set target box Z position").SetUnit("cm");
  _messenger->DeclareProperty("targetSizeX", fTargetSizeX, "Set target box X size").SetUnit("cm");
  _messenger->DeclareProperty("targetSizeY", fTargetSizeY, "Set target box Y size").SetUnit("cm");
  _messenger->DeclareProperty("targetSizeZ", fTargetSizeZ, "Set target box Z size").SetUnit("cm");
  _messenger->DeclareProperty("targetLayerOffsetZ", fTargetLayerOffsetZ, "Set target layer Z offset").SetUnit("cm");
  _messenger->DeclareProperty("targetLayerDistance", fTargetLayerDistance, "Set target inter-layer distance").SetUnit("cm");
  _messenger->DeclareProperty("targetLayerNo", fTargetLayerNo, "Set target layer number");
 
  _messenger->DeclareProperty("tileX", fTileX, "Set sensor tile X size").SetUnit("cm");
  _messenger->DeclareProperty("tileY", fTileY, "Set sensor tile Y size").SetUnit("cm");
  _messenger->DeclareProperty("tileThickness", fTileThickness, "Set sensor tile thickness").SetUnit("mm");
  _messenger->DeclareProperty("tileSPitch", fTileSPitch, "Set sensor tile S Pitch").SetUnit("mm");
  _messenger->DeclareProperty("tileKPitch", fTileKPitch, "Set sensor tile K Pitch").SetUnit("mm");
   
  _messenger->DeclareProperty("targetLayerMiniNo", fTargetLayerMiniNo, "Set target mini-layer number");
  _messenger->DeclareProperty("targetLayerShortNo", fTargetLayerShortNo, "Set target short-layer number");
  _messenger->DeclareProperty("targetLayerDampeNo", fTargetLayerDampeNo, "Set target Dampe-layer number");
  _messenger->DeclareProperty("layerMiniTileNo", fLayerMiniTileNo, "Set mini-layer tiles number");

  _messenger->DeclareProperty("layerShortTileNo", fLayerShortTileNo, "Set short-layer tiles number");
  _messenger->DeclareProperty("layerLongTileNo", fLayerLongTileNo, "Set long-layer tiles number");
  _messenger->DeclareProperty("layerDampeTileNo", fLayerDampeTileNo, "Set dampe-layer tiles number");
  
    _messenger->DeclareProperty("smOffsetZ", fSMOffsetZ, "Set spectrometer box Z offset position").SetUnit("cm");
  _messenger->DeclareProperty("smSizeX", fSMSizeX, "Set spectrometer box X size").SetUnit("cm");
  _messenger->DeclareProperty("smSizeY", fSMSizeY, "Set spectrometer box Y size").SetUnit("cm");
  _messenger->DeclareProperty("smSizeZ", fSMSizeZ, "Set spectrometer box Z size").SetUnit("cm");

  _messenger->DeclareProperty("magVolR", fMagVolR, "Set magnetic volume Radius").SetUnit("cm");
  _messenger->DeclareProperty("magVolZ", fMagVolZ, "Set magnetic volume Z size").SetUnit("cm");
  
  _messenger->DeclareProperty("magFieldVal", fMagFieldVal, "Set Constant Magnetic Field Value").SetUnit("tesla");

  _messenger->DeclareProperty("smLayerNo", fSMLayerNo, "Set spectrometer layers number");
  _messenger->DeclareProperty("layerFDist", fLayerFDistance, "Set sm Front inter-layer distance").SetUnit("cm");
  _messenger->DeclareProperty("layerEDist", fLayerEDistance, "Set sm Exit inter-layer distance").SetUnit("cm");
  _messenger->DeclareProperty("layerFGap", fLayerFGap, "Set Front layer-magvol distance").SetUnit("cm");
  _messenger->DeclareProperty("layerEGap", fLayerEGap, "Set Exit layer-magvol distance").SetUnit("cm");

  
  //  detMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  //delete detMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() {
  //G4NistManager* man = G4NistManager::Instance();

  //  G4bool isotopes = false;
  //  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4cout << "begin of Detector Construction" << G4endl;
  //Run geometry configuration script before building
  if (_geoDataCard != "") {
    G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute " + _geoDataCard));
  }

  // Delete the messenger so that the commands for configuring the geometry won't be available anymore
  delete _messenger;
  _messenger = NULL;
  
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_Galactic");

  G4double world_diameter = 1000. * cm;
  G4Sphere* solidWorld = new G4Sphere("World", 0., 0.5 * world_diameter, 0., 2*pi, 0., pi); //For compatibility with VGM [V.F.]

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,          //its solid
      default_mat,         //its material
      "World");            //its name

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  fPhysicalWorld = new G4PVPlacement(0,                     //no rotation
      G4ThreeVector(),       //at (0,0,0)
      logicWorld,            //its logical volume
      "World",               //its name
      0,                     //its mother  volume
      false,                 //no boolean operation
      0,                     //copy number
      fCheckOverlaps);       // checking overlaps


  //  G4double l = 0.0001 * mm; //tolerance  

  //// HARDCODED ALIGNMENT -> to refine
  //// maybe with a tbitsarray like stuff?
  
  // G4cout<< "**** warning-> HARDCODED DETECTOR ALIGNMENT (X=0 Y=1): ";
  // layers alignment (X=0 Y=1)

  
  G4cout<< "**** STRING DETECTOR ALIGNMENT: "<<fAlign<<G4endl;
  //int * al = new int[fAlign.size()];
  //std::copy(fAlign.begin(), fAlign.end(), al); 

  G4int tAlign[fTargetLayerNo];
  for(int it=0;it<fTargetLayerNo;it++){
    tAlign[it]=fAlign[it]-'0';
    //tAlign[it]=it%2==0?0:1;
    G4cout<<tAlign[it]<<" ";
   }

  // spectrometer alignment (X=0 Y=1)
  G4int sAlign[fSMLayerNo];
  sAlign[0]=1;
  G4cout<<sAlign[0]<<" ";
  for(int is=1;is<fSMLayerNo;is++){
    sAlign[is]=fAlign[fTargetLayerNo+is]-'0';
    G4cout<<sAlign[is]<<" ";
  }
  G4cout<<G4endl;

  
  G4RotationMatrix* myRotation = new G4RotationMatrix();
  myRotation->rotateX(0.*deg);
  myRotation->rotateY(0.*deg);
  myRotation->rotateZ(90.*deg);
  G4double layerMiniX = fTileX * fLayerMiniTileNo;  // x of mini ladder 
  G4double layerShortX = fTileX * fLayerShortTileNo; // x of short ladder
  G4double layerLongX = fTileX * fLayerLongTileNo; // x of long ladder
  G4double layerDampeX = fTileX * fLayerDampeTileNo; // x of dampe ladder 
      
  G4double siStripX = fTileX;
  G4double siStripY = fTileSPitch;
  G4double siStripZ = fTileThickness;

  G4double targetLayerFirstZ = fTargetLayerOffsetZ-0.5*fTargetSizeZ;

   
  G4Box* targetMother = new G4Box("target", 0.5 * fTargetSizeX, 0.5 * (fTargetSizeY), 0.5 * (fTargetSizeZ));
  G4LogicalVolume* targetLog = new G4LogicalVolume(targetMother, default_mat, "target");

  new G4PVPlacement(0,                     //no rotation
		    G4ThreeVector(0.,0.,fTargetPosZ),       //at (0,0,0)
		    targetLog,              //its logical volume
		    "target",                 //its name
		    logicWorld,            //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    fCheckOverlaps);       // checking overlaps
  

  
  // check if needed
  /*
  G4Box* padMother = new G4Box("pad", 0.5 * (pad_x + l), 0.5 * (pad_y + l), 0.5 * (pad_z + l));
  G4LogicalVolume* padLogic = new G4LogicalVolume(padMother, default_mat, "pad");
  new G4PVPlacement(0,                     //no rotation
    G4ThreeVector(),       //at (0,0,0)
    padLogic,              //its logical volume
    "pad",                 //its name
    logicWorld,            //its mother  volume
    false,                 //no boolean operation
    0,                     //copy number
    fCheckOverlaps);       // checking overlaps
  */
  
  G4Box* siLayerMini = new G4Box("siLayerMini", 0.5 * layerMiniX, 0.5 * fTileY, 0.5 * fTileThickness);
  G4Box* siLayerShort = new G4Box("siLayerShort", 0.5 * layerShortX, 0.5 * fTileY, 0.5 * fTileThickness);
  G4Box* siLayerLong = new G4Box("siLayerLong", 0.5 * layerLongX, 0.5 * fTileY, 0.5 * fTileThickness);
 G4Box* siLayerDampe = new G4Box("siLayerDampe", 0.5 * layerDampeX, 0.5 * fTileY, 0.5 * fTileThickness);
  

  G4LogicalVolume* siLayerMiniLog = new G4LogicalVolume(siLayerMini, silicon, "siLayerMiniLog");
  G4LogicalVolume* siLayerShortLog = new G4LogicalVolume(siLayerShort, silicon, "siLayerShortLog");
  G4LogicalVolume* siLayerLongLog = new G4LogicalVolume(siLayerLong, silicon, "siLayerLongLog");
  G4LogicalVolume* siLayerDampeLog = new G4LogicalVolume(siLayerDampe, silicon, "siLayerDampeLog");
  
  for (G4int iln = 0; iln <fTargetLayerNo ; iln++)
    {
      if (iln<fTargetLayerMiniNo)
	new G4PVPlacement(!tAlign[iln]?myRotation:0,G4ThreeVector(0.,0.,targetLayerFirstZ
								  +fTileThickness/2 
								  +(iln)*fTargetLayerDistance),
			  siLayerMiniLog,
			  "siLayerMiniPhys",
			  targetLog,
			  false,	
			  iln,
			  false);
      else if (iln-fTargetLayerMiniNo<fTargetLayerShortNo)
	new G4PVPlacement(!tAlign[iln]?myRotation:0,G4ThreeVector(tAlign[iln]*fLayersOffsetY,(!tAlign[iln])*fLayersOffsetX,targetLayerFirstZ
								  +fTileThickness/2 
								  +(iln)*fTargetLayerDistance),
			  siLayerShortLog,
			  "siLayerShortPhys",
			  targetLog,
			  false,	
			  iln,
			  false);
      
      else
	new G4PVPlacement(!tAlign[iln]?myRotation:0,G4ThreeVector(tAlign[iln]*fLayersOffsetY,(!tAlign[iln])*fLayersOffsetX,targetLayerFirstZ
								  +fTileThickness/2 
								  +(iln)*fTargetLayerDistance),
			  siLayerDampeLog,
			  "siLayerDampePhys",
			  targetLog,
			  false,	
			  iln,
			  false);
    }
  

  G4Box* siTile = new G4Box("siTile", 0.5 * fTileX, 0.5 * fTileY, 0.5 * fTileThickness);
  G4LogicalVolume* siTileLog = new G4LogicalVolume(siTile, silicon, "siTileLog");
  
  G4PVReplica * layerMiniReplica = new G4PVReplica("miniReplica", siTileLog,siLayerMiniLog,kXAxis, fLayerMiniTileNo, fTileX, 0);
  
  G4PVReplica * layerShortReplica = new G4PVReplica("shortReplica", siTileLog,siLayerShortLog,kXAxis, fLayerShortTileNo, fTileX, 0);

  G4PVReplica * layerLongReplica = new G4PVReplica("longReplica", siTileLog,siLayerLongLog,kXAxis, fLayerLongTileNo, fTileX, 0);

  G4PVReplica * layerDampeReplica = new G4PVReplica("dampeReplica", siTileLog,siLayerDampeLog,kXAxis, fLayerDampeTileNo, fTileX, 0);

  
  G4double spectrometerPosZ = fTargetSizeZ/2. + fSMSizeZ/2. + fSMOffsetZ;
  
  G4Box* spectrometerMother = new G4Box("spectrometer", 0.5*fSMSizeX, 0.5*fSMSizeY, 0.5*fSMSizeZ);
  G4LogicalVolume* spectrometerLog = new G4LogicalVolume(spectrometerMother, default_mat, "spectrometer");

  new G4PVPlacement(0,                     //no rotation
		    G4ThreeVector(0.,0.,spectrometerPosZ),       //at (0,0,0)
		    spectrometerLog,              //its logical volume
		    "spectrometer",                 //its name
		    logicWorld,            //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    fCheckOverlaps);       // checking overlaps
  

  
  G4double spectrometerFPosZ = - fMagVolZ/2. - fLayerFGap - fTileThickness/2.- fLayerFDistance;
  G4double spectrometerEPosZ = fMagVolZ/2. + fLayerEGap + fTileThickness/2. ;

    
  // HARDCODED LAYER NUMBER -->to change
  // Front layers PhysVol 
  
  new G4PVPlacement(!sAlign[0]?myRotation:0,G4ThreeVector((sAlign[0])*fLayersOffsetY,(!sAlign[0])*fLayersOffsetX,spectrometerFPosZ),siLayerLongLog,"siLayerF1Phys",spectrometerLog,false,fTargetLayerNo,false);
  
  new G4PVPlacement(!sAlign[1]?myRotation:0,G4ThreeVector((sAlign[1])*fLayersOffsetY,(!sAlign[1])*fLayersOffsetX,spectrometerFPosZ+fLayerFDistance),siLayerLongLog,"siLayerF2Phys",spectrometerLog,false,fTargetLayerNo+1,false);
  
  // End layers PhysVol 
 
  new G4PVPlacement(!sAlign[2]?myRotation:0,G4ThreeVector((sAlign[2])*fLayersOffsetY,(!sAlign[2])*fLayersOffsetX,spectrometerEPosZ),siLayerLongLog,"siLayerE1Phys",spectrometerLog,false,fTargetLayerNo+2,false);
  
  new G4PVPlacement(!sAlign[3]?myRotation:0,G4ThreeVector((sAlign[3])*fLayersOffsetY,(!sAlign[3])*fLayersOffsetX,spectrometerEPosZ+fLayerEDistance),siLayerLongLog,"siLayerE2Phys",spectrometerLog,false,fTargetLayerNo+3,false);
  // 


  G4int nSStrips = int(fTileY/fTileSPitch); // = siStripY
  G4int nKStrips = int(fTileY/fTileKPitch);
  
  G4Box* siStrip = new G4Box("siStrip", 0.5 * siStripX, 0.5 * siStripY, 0.5 * siStripZ);
  G4LogicalVolume* siStripLog = new G4LogicalVolume(siStrip, silicon, "siStripLog");

  // G4PVReplica * siStripReplica = new G4PVReplica("stripReplica", siStripLog, siTileLog, kYAxis, nSStrips, siStripY, 0);
  

  // magnet vol
  
  G4Tubs* magnetMother = new G4Tubs("magnet", 0., fMagVolR, 0.5*fMagVolZ,0., 2*pi);
  G4LogicalVolume* magnetLog = new G4LogicalVolume(magnetMother, default_mat, "magnet");
  
  new G4PVPlacement(0,                     //no rotation
		    G4ThreeVector(0.,0.,0.),       //at (0,0,0)
		    magnetLog,              //its logical volume
		    "magnet",                 //its name
		    spectrometerLog,            //its mother  volume
		    false,                 //no boolean operation
		    0,                     //copy number
		    fCheckOverlaps);       // checking overlaps
  

  //  G4FieldManager* fieldMgr
  //  = G4TransportationManager::GetTransportationManager()
  //  ->GetFieldManager();
  //fieldMgr->SetDetectorField(0);

  fMagField = new G4UniformMagField(G4ThreeVector(fMagFieldVal,0.,0.));
  
  G4FieldManager* localFieldMgr= new G4FieldManager(); 
  localFieldMgr->SetDetectorField(fMagField);
  localFieldMgr->CreateChordFinder(fMagField);
  magnetLog->SetFieldManager(localFieldMgr,false);

  //                                        
  // Visualization attributes
  //  
  
  // Some visualization styles

  G4VisAttributes* VisAtt1= new G4VisAttributes(G4Colour(0.3,0.8,0.1));
  VisAtt1->SetVisibility(true);
  VisAtt1->SetForceSolid(TRUE);

  G4VisAttributes* VisAtt2= new G4VisAttributes(G4Colour(0.2,0.3,0.8));
  VisAtt2->SetVisibility(true);
  VisAtt2->SetForceSolid(FALSE);

  G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.8,0.2,0.3));
  VisAtt3->SetVisibility(true);
  VisAtt3->SetForceWireframe(TRUE);  

  siLayerMiniLog->SetVisAttributes(VisAtt1);
  siLayerShortLog->SetVisAttributes(VisAtt1);
  siLayerLongLog->SetVisAttributes(VisAtt1);

  magnetLog->SetVisAttributes(VisAtt2);
  
  G4VisAttributes * Yellow     = new G4VisAttributes(G4Colour(255/255., 255/255.,0/255.));     
  //  siStripLog->SetVisAttributes(Yellow);
  siTileLog->SetVisAttributes(Yellow);

  
  // try parameterization
  //  siLayerParam = new SiLayerParameterisation (npads, zfirst, paddist, pad_z);
  
  // G4PhysicalVolume * siLayerPhys =
  //    new G4PVParameterized("target", siLayerLog, targetLog, kZAxis, npads, siLayerParam);
  
  
  //always return the physical World
  return fPhysicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void DetectorConstruction::ConstructSDandField() {
  //G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  /*
   G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
   G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
   cryst->RegisterPrimitive(primitiv1);
   SetSensitiveDetector("Crystal",cryst);
   */
  //G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  //G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  //cryst->RegisterPrimitive(primitiv1);
  //SetSensitiveDetector("Crystal", cryst);
// }

void DetectorConstruction::updateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}


bool DetectorConstruction::ExportParameters() {
  bool result = true;
  result = result && ExportStringParameter("layerAlignment", fAlign);
  result = result && ExportRealParameter("layersOffsetX", fLayersOffsetX / cm);
  result = result && ExportRealParameter("layersOffsetY", fLayersOffsetY / cm);
  result = result && ExportRealParameter("targetPosZ", fTargetPosZ / cm);
  result = result && ExportRealParameter("targetSizeX", fTargetSizeX / cm);
  result = result && ExportRealParameter("targetSizeY", fTargetSizeY / cm);
  result = result && ExportRealParameter("targetSizeZ", fTargetSizeZ / cm);
  result = result && ExportRealParameter("targetLayerOffsetZ", fTargetLayerOffsetZ / cm);
  result = result && ExportIntParameter("targetLayerNo", fTargetLayerNo);
  result = result && ExportIntParameter("targetLayerMiniNo", fTargetLayerMiniNo);
  result = result && ExportIntParameter("targetLayerShortNo", fTargetLayerShortNo);
  result = result && ExportIntParameter("targetLayerDampeNo", fTargetLayerDampeNo);
  result = result && ExportRealParameter("targetLayerDistance", fTargetLayerDistance / cm);
  result = result && ExportRealParameter("tileX", fTileX / cm);
  result = result && ExportRealParameter("tileY", fTileY / cm);
  result = result && ExportRealParameter("tileThickness", fTileThickness / mm);
  result = result && ExportRealParameter("tileSPitch", fTileSPitch / mm);
  result = result && ExportRealParameter("tileKPitch", fTileKPitch / mm);
  result = result && ExportIntParameter("layerMiniTileNo", fLayerMiniTileNo);
  result = result && ExportIntParameter("layerShortTileNo", fLayerShortTileNo);
  result = result && ExportIntParameter("layerLongTileNo", fLayerLongTileNo);
  result = result && ExportIntParameter("layerDampeTileNo", fLayerDampeTileNo);
  result = result && ExportRealParameter("smOffsetZ", fSMOffsetZ / cm);
  result = result && ExportRealParameter("smSizeX", fSMSizeX / cm);
  result = result && ExportRealParameter("smSizeY", fSMSizeY / cm);
  result = result && ExportRealParameter("smSizeZ", fSMSizeZ / cm);
  result = result && ExportRealParameter("magVolR", fMagVolR / cm);
  result = result && ExportRealParameter("magVolZ", fMagVolZ / cm);
  result = result && ExportRealParameter("magFieldVal", fMagFieldVal / tesla);
  result = result && ExportIntParameter("smLayerNo", fSMLayerNo);
  result = result && ExportRealParameter("layerFDist", fLayerFDistance / cm);
  result = result && ExportRealParameter("layerEDist", fLayerEDistance / cm);
  result = result && ExportRealParameter("layerFGap", fLayerFGap / cm);
  result = result && ExportRealParameter("layerEGap", fLayerEGap / cm);
  return result;
}

const std::string DetectorConstruction::GetVersion() {
  return "1.0";
}


/*
void DetectorConstruction::SetMagFieldVal(G4double v)
{
  printf("Setting constant value of magnetic field to %f Tesla \n",v/tesla);
  fMagField->SetFieldValue(G4ThreeVector(v,0.,0.));
}
*/
/*void DetectorConstruction::SetTargetX(G4double tx)
{
  printf("Setting target size X to %f cm \n",tx);
  det->SetTargetX(tx);
}
*/
