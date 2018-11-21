#ifndef DVCS_DETECTOR_CONSTRUCTION_HH
#define DVCS_DETECTOR_CONSTRUCTION_HH 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class dvcsDetectorConstruction: public G4VUserDetectorConstruction
{
public:
  dvcsDetectorConstruction();
  ~dvcsDetectorConstruction();
  
  G4VPhysicalVolume *Construct();

  G4VPhysicalVolume *GetHRS_window() {return HRS_window_phys;};
  
private:
  
  //Logical Volumes
  G4LogicalVolume *world_log;
  G4LogicalVolume *LH2_target_log;
  G4LogicalVolume *target_container_log;
  G4LogicalVolume *Al_scat_chamber_log;
  G4LogicalVolume *Vac_scat_chamber_log;
  G4LogicalVolume *air_window_log;
  G4LogicalVolume *kapton_window_log;
  G4LogicalVolume *nouse_shield_log;
  G4LogicalVolume *beamline_shield_log;
  G4LogicalVolume *HRS_window_log;

  //Physical Volumes
  G4VPhysicalVolume *world_phys;  // Experimental Hall
  G4VPhysicalVolume *LH2_target_phys; // target Liquid Hidrogen
  G4VPhysicalVolume *target_container_phys;  // target Container Aluminum
  G4VPhysicalVolume *Al_scat_chamber_phys; // Aluminium walls of Scattering Chamber
  G4VPhysicalVolume *Vac_scat_chamber_phys; // Vacuum inside of Scattering Chamber
  G4VPhysicalVolume *air_window_phys; // Kapton Window
  G4VPhysicalVolume *kapton_window_phys; // Kapton Window
  G4VPhysicalVolume *nouse_shield_phys; // nouse_shield
  G4VPhysicalVolume *beamline_shield_phys; // beamline_shield
  G4VPhysicalVolume *HRS_window_phys;       // HRS front window, in order to detect scettered electron, and it's energy loss
};


#endif
