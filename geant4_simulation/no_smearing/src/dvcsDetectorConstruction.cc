#include "dvcsDetectorConstruction.hh"
#include "DVCSCaloConstruction.hh"
#include "DVCSCaloGeom.hh"

//#include "TMath.h"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "dvcsHRS_windowSD.hh"
#include "dvcsGlobals.hh"
#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


using namespace std;

dvcsDetectorConstruction::dvcsDetectorConstruction():
  world_log(0), LH2_target_log(0), target_container_log(0), Al_scat_chamber_log(0), Vac_scat_chamber_log(0), air_window_log(0), kapton_window_log(0), 
  nouse_shield_log(0), beamline_shield_log(0), HRS_window_log(0),
  world_phys(0), LH2_target_phys(0), target_container_phys(0), Al_scat_chamber_phys(0), Vac_scat_chamber_phys(0), air_window_phys(0), 
  kapton_window_phys(0), nouse_shield_phys(0), beamline_shield_phys(0), HRS_window_phys(0)
{
  ;
}

dvcsDetectorConstruction::~dvcsDetectorConstruction()
{
  ;
}

G4VPhysicalVolume* dvcsDetectorConstruction::Construct()
{
  //=================Define Materials==============================
  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel;
  
  G4double vac_atomicNumber = 1.;
  G4double vac_massOfMole = 1.008*g/mole;
  G4double vac_density = 1.e-25*g/cm3;
  G4double vac_temperature = 2.73*kelvin;
  G4double vac_pressure = 3.e-18*pascal;
  G4Material* Vac_mat = new G4Material("interGalactic", vac_atomicNumber, vac_massOfMole, vac_density, kStateGas, vac_temperature, vac_pressure);

  G4NistManager *nist_man = G4NistManager::Instance();

  //Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  

  G4Material *Al_mat = nist_man->FindOrBuildMaterial("G4_Al");
  G4Material *LH2_mat = nist_man->FindOrBuildMaterial("G4_lH2");
  G4Material *Kapton_mat = nist_man->FindOrBuildMaterial("G4_KAPTON");
  G4Material *Wolfram_mat = nist_man->FindOrBuildMaterial("G4_W");
  
  //===============Some quantities====================================

  G4double inch = 2.54*cm;
  G4double radian = 57.2957795130823229;
  
  G4double world_lx = 10*m;
  G4double world_ly = 10*m;
  G4double world_lz = 10*m;

  G4double sc_radius=63.*cm;
  G4double sc_or = 24.51*inch;
  G4double sc_ir = 24.13*inch;
  G4double sc_thickness = sc_or - sc_ir;
  //G4double sc_ir=sc_radius-sc_thickness/2.;
  //G4double sc_or=sc_radius+sc_thickness/2.;
  
  G4double dt_thickness=0.28*inch;

  G4double dt_or = 0.5*6.625*inch;
  G4double dt_ir=dt_or-dt_thickness;
  G4double dt_length=100.*cm;

  G4double ut_radius=3.*cm;
  G4double ut_thickness=0.3*cm;
  G4double ut_ir=ut_radius-ut_thickness/2.;
  G4double ut_or=ut_radius+ut_thickness/2.;
  G4double ut_length=100.*cm;

  G4double target_LH2_length = dvcsGlobals::target_length*cm;
  //G4double target_LH2_length = 15*cm;
  G4double target_wall_thickness = 0.141*mm;
  G4double entrance_window_thick = 0.128*mm;
  G4double exit_window_thick = 0.207*mm;
  G4double target_container_radius = 0.5*1.6*2.54*cm;
  G4double target_container_length = target_LH2_length + entrance_window_thick + exit_window_thick;
  G4double target_tube_length = target_container_length - target_container_radius;
  
  G4double target_LH2_radius = target_container_radius - target_wall_thickness;
  G4double target_LH2_tube_length = target_LH2_length - target_LH2_radius;

  G4double target_offset = dvcsGlobals::target_offset*cm;
  G4double target_container_tube_center = target_offset - target_container_radius/2.;
  G4double target_LH2_tube_center = -target_tube_length/2. + entrance_window_thick + target_LH2_tube_length/2.;
  
  G4double dt_shift = sqrt( sc_ir*sc_ir - dt_ir*dt_ir );
  G4double ut_shift = sqrt( sc_ir*sc_ir - ut_ir*ut_ir );


  G4double kapton_thickness = 0.02*inch;
  G4double kapton_height = 4.24*inch;
  G4double kapton_vert_angle = 2*asin(0.5*kapton_height/sc_ir)*radian*degree;
  G4double kapton_horiz_angle = 36*degree;
  G4double kapton_horiz_shift = 28*degree;

  G4double nouse_shield_lz = 3.94*inch; // length along z direction
  G4double nouse_shield_height = 13.39*inch;
  G4double nouse_shield_l_thick = 0.84*inch; // large thickness
  G4double nouse_shield_angle = 7.97*degree;
  G4double nouse_shield_survay_angle = 9.5*degree;
  G4double nouse_shield_s_thick =  nouse_shield_l_thick - nouse_shield_lz*tan(nouse_shield_angle);
  G4double nouse_shield_z_shift = 68.3*cm + nouse_shield_lz/2.;
  //G4double nouse_shield_x_shift = -nouse_shield_z_shift*tan(nouse_shield_survay_angle) + 0.5*0.5*(nouse_shield_l_thick + nouse_shield_s_thick);
  G4double nouse_shield_x_shift = -1.01*dt_or - 0.5*0.5*(nouse_shield_l_thick + nouse_shield_s_thick);
  
  G4double nouse_shield_roll_angle = -0.3*degree; //clockwise about the beam's axis
  G4double nouse_shield_xrot = 90*degree;
  G4double nouse_shield_yrot = 180*degree;
  
  
  G4double beamline_shield_lz = 10.236*inch;
  G4double beamline_shield_height = 7.087*inch;
  G4double beamline_shield_thick = 0.591*inch;

  G4double beamline_shield_z_shift = nouse_shield_z_shift + 0.5*(nouse_shield_lz + beamline_shield_lz);
  //G4double beamline_shiels_x_shift = nouse_shield_x_shift + 0.5*0.5*(nouse_shield_l_thick + nouse_shield_s_thick) - 0.5*beamline_shield_thick;
  G4double beamline_shiels_x_shift = -1.01*dt_or - 0.5*beamline_shield_thick;
  
  
  G4double HRS_window_horiz_spanning = 2*60*milliradian;
  G4double HRS_window_vert_spanning = 2*70*milliradian;
  //G4double HRS_window_horiz_shift = 21.49*degree;
  G4double HRS_window_horiz_shift = dvcsGlobals::HRS_angle*rad;
  G4double HRS_window_vert_shift = 0*degree;
  G4double HRS_window_radius = 1.1086*m;

  //  G4double dt_shift = sc_ir - sqrt( sc_ir*sc_ir - dt_ir*dt_ir ); //down tube should be shifted by this quantity from z = sc_ir point 
                                                                 //to the upstream.
  //G4double ut_shift = sc_ir - sqrt( sc_ir*sc_ir - ut_ir*ut_ir ); //up tube should be shifted by this quantity from z = sc_ir point 
                                                                 //to the downstream.
 
  //==================world====================
  
  G4Box *world_box = new G4Box("world_box", world_lx, world_ly, world_lz);
  world_log = new G4LogicalVolume(world_box, Air, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0., 0.,0. ), world_log, 
				 "world_phys", 0, 0, false);
  world_log->SetVisAttributes(G4VisAttributes::Invisible);

  
  //==================Scattering Chamber with up and down tubes and Capton/Air windows=========================

  G4ThreeVector dt_pos(0, 0, +dt_shift + dt_length/2.);
  G4RotationMatrix *dt_rot = new G4RotationMatrix();

  G4ThreeVector ut_pos(0, 0, -ut_shift - ut_length/2.);
  G4RotationMatrix *ut_rot = new G4RotationMatrix();
  
  //==========================================================================================
  //-------------------A newer Method for construction of Scattering chamber------------------
  
  //---------Aluminum Scattering Chamber----------
  G4VSolid *Al_orb = new G4Orb("Al_orb", sc_or);
  G4VSolid *Al_down_tube = new G4Tubs("Al_down_tube", 0, dt_or, dt_length/2, 0, 360*degree);
  G4VSolid *Al_uper_tube = new G4Tubs("Al_uper_tube", 0, ut_or, ut_length/2, 0, 360*degree);

  G4VSolid *Al_intermed_chamber = new G4UnionSolid("Al_intermed_chamber", Al_orb, Al_down_tube, dt_rot, dt_pos);
  G4VSolid *Al_scat_chamber = new G4UnionSolid("Al_scat_chamber", Al_intermed_chamber, Al_uper_tube, ut_rot, ut_pos);


  //---------Vacuum Scattering Chamber----------
  G4VSolid *Vac_orb = new G4Orb("Vac_orb", sc_ir);
  G4VSolid *Vac_down_tube = new G4Tubs("Vac_down_tube", 0, dt_ir, dt_length/2, 0, 360*degree);
  G4VSolid *Vac_uper_tube = new G4Tubs("Vac_uper_tube", 0, ut_ir, ut_length/2, 0, 360*degree);

  G4VSolid *Vac_intermed_chamber = new G4UnionSolid("Vac_intermed_chamber", Vac_orb, Vac_down_tube, dt_rot, dt_pos);
  G4VSolid *Vac_scat_chamber = new G4UnionSolid("Vac_scat_chamber", Vac_intermed_chamber, Vac_uper_tube, ut_rot, ut_pos);

  Al_scat_chamber_log = new G4LogicalVolume(Al_scat_chamber, Al_mat, "Al_scat_chamber");

  Al_scat_chamber_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Al_scat_chamber_log, "Al_scat_chamber_phys", world_log, true, 0, true);
  
  Vac_scat_chamber_log = new G4LogicalVolume(Vac_scat_chamber, Vac_mat, "Vac_scat_chamber");
  Vac_scat_chamber_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Vac_scat_chamber_log, 
					    "Vac_scat_chamber_phys", Al_scat_chamber_log, true, 0, true);

  
  //==========================================================================================
  
  /*
  G4Sphere *scat_sphere = new G4Sphere("scat_sphere", sc_ir, sc_or, 0, 360.*degree, 0, 180.*degree);
  G4Tubs *down_tube = new G4Tubs("down_tube", 0, dt_or, dt_length/2., 0, 360.*degree);
  G4Tubs *down_vac_tube = new G4Tubs("down_vac_tube", 0, dt_ir, dt_length/2., 0, 360.*degree); // down Vacum tube
  G4Tubs *up_tube = new G4Tubs("up_tube", 0, ut_or, ut_length/2., 0, 360.*degree);
  G4Tubs *up_vac_tube = new G4Tubs("up_vac_tube", 0, ut_ir, ut_length/2., 0, 360.*degree); // down Vacum tube
  
  G4UnionSolid *sphere_downtube_intermed = new G4UnionSolid("sphere_downtube_intermed", scat_sphere, down_tube, dt_rot, dt_pos );
  
  G4SubtractionSolid *sphere_downtube = new G4SubtractionSolid("sphere_downtube", sphere_downtube_intermed, down_vac_tube, dt_rot, dt_pos);
  
  G4UnionSolid *sphere_uptube_intermed = new G4UnionSolid("sphere_uptube_intermed", sphere_downtube, up_tube, ut_rot, ut_pos);
  
  G4SubtractionSolid *scat_spher = new G4SubtractionSolid("scat_spher", sphere_uptube_intermed, up_vac_tube, ut_rot, ut_pos);
  
  scat_chamber_log = new G4LogicalVolume(scat_spher, Al_mat, "scat_chamber");
  scat_chamber_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), scat_chamber_log, "scat_chamber_phys", world_log, true, 0, true);
  */
  // G4VisAttributes *scat_chamber_vis = new G4VisAttributes(G4Colour(0.1, 0.2, 0.9));
  // scat_chamber_vis->SetForceWireframe(0);
  // scat_chamber_vis->SetForceSolid(1);
  // Al_scat_chamber_log->SetVisAttributes(scat_chamber_vis);
  Al_scat_chamber_log->SetVisAttributes(G4VisAttributes::Invisible);

  // G4VisAttributes *vac_scat_chamber_vis = new G4VisAttributes(G4Colour(0.99, 0.99, 0.02));
  // vac_scat_chamber_vis->SetForceWireframe(0);
  // vac_scat_chamber_vis->SetForceSolid(1);
  // Vac_scat_chamber_log->SetVisAttributes(vac_scat_chamber_vis);
  Vac_scat_chamber_log->SetVisAttributes(G4VisAttributes::Invisible);
  //===============================Air Window=============================

  G4Sphere *air_window = new G4Sphere("air_window", sc_ir, sc_ir + sc_thickness, -kapton_horiz_angle/2, kapton_horiz_angle, 
					 90*degree - kapton_vert_angle/2., kapton_vert_angle);


  //G4Sphere *air_window = new G4Sphere("air_window", sc_ir, sc_ir + sc_thickness, -kapton_vert_angle/2, kapton_vert_angle, 
  //				      kapton_horiz_shift, kapton_horiz_angle);
  air_window_log = new G4LogicalVolume(air_window, Air, "air_window_log");
  
  G4RotationMatrix *air_rot = new G4RotationMatrix();
  air_rot->rotateX(90*degree);
  air_rot->rotateZ(90*degree);
  air_rot->rotateZ(-kapton_horiz_shift);
  air_window_phys = new G4PVPlacement(air_rot, G4ThreeVector(0., 0., 0.), air_window_log, "air_window_phys", Al_scat_chamber_log, true, 0, true);

  G4VisAttributes *air_window_vis = new G4VisAttributes(G4Colour(0.01, 0.99, 0.99));
  air_window_vis->SetForceWireframe(0);
  air_window_vis->SetForceSolid(1);
  air_window_log->SetVisAttributes(air_window_vis);
  //air_window_log->SetVisAttributes(G4VisAttributes::Invisible);

  //=============================Capton Window============================
  
  G4Sphere *kapton_window = new G4Sphere("kapton_window", sc_ir, sc_ir + kapton_thickness, -kapton_horiz_angle/2, kapton_horiz_angle, 
					 90*degree - kapton_vert_angle/2., kapton_vert_angle);
  //G4Sphere *kapton_window = new G4Sphere("kapton_window", sc_ir, sc_ir + kapton_thickness, -kapton_vert_angle/2, kapton_vert_angle, 
  //kapton_horiz_shift, kapton_horiz_angle);
  
  kapton_window_log = new G4LogicalVolume(kapton_window, Kapton_mat, "kapton_window_log");
  
  G4RotationMatrix *kapton_rot = new G4RotationMatrix();
  //kapton_rot->rotateX(90*degree);
  //kapton_rot->rotateY(90*degree);
  //kapton_rot->rotateY(-kapton_horiz_shift);
  

  kapton_window_phys = new G4PVPlacement(kapton_rot, G4ThreeVector(0., 0., 0.), kapton_window_log, "kapton_window_phys", air_window_log, 
					 true, 0, true);
  G4VisAttributes *kapton_window_vis = new G4VisAttributes(G4Colour(0.99, 0.02, 0.99));
  kapton_window_vis->SetForceWireframe(0);
  kapton_window_vis->SetForceSolid(1);
  kapton_window_log->SetVisAttributes(kapton_window_vis);


  //====================target====================
  
  //=======container=======
  G4Tubs *target_container_tube = new G4Tubs("target_container_tube", 0, target_container_radius, target_tube_length/2., 0, 360.*degree);
  G4Sphere *target_container_endcup = new G4Sphere("target_container_endcup", 0, target_container_radius, 0.*degree, 360.*degree, 0., 90*degree);
  G4UnionSolid *target_container_solid = new G4UnionSolid("target_container_solid", target_container_tube, target_container_endcup, 
							  0, G4ThreeVector(0., 0., 0.999*(target_tube_length/2.)));
  target_container_log = new G4LogicalVolume(target_container_solid, Al_mat, "target_container_log");
  target_container_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., target_container_tube_center), 
  target_container_log, "target_container_phys", Vac_scat_chamber_log, true, 0, true);
  //=======LH2=======
  G4Tubs *target_LH2_tube = new G4Tubs("target_LH2_tube", 0, target_LH2_radius, target_LH2_tube_length/2., 0, 360.*degree);
  G4Sphere *target_LH2_endcup = new G4Sphere("target_LH2_endcup", 0, target_LH2_radius, 0., 360.*degree, 0., 90.*degree);
  G4UnionSolid *target_LH2_solid = new G4UnionSolid("target_LH2_solid", target_LH2_tube, target_LH2_endcup, 0, 
						    G4ThreeVector(0, 0, 0.999*(target_LH2_tube_length/2.)));
  LH2_target_log = new G4LogicalVolume(target_LH2_solid, LH2_mat, "LH2_target_log");
  LH2_target_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., target_LH2_tube_center), LH2_target_log, "LH2_target_phys", target_container_log,
				      true, 0, true);
  
  G4VisAttributes *target_container_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  target_container_vis->SetForceWireframe(1);
  //target_container_vis->SetForceSolid(1);
  target_container_log->SetVisAttributes(target_container_vis);
  //target_container_log->SetVisAttributes(G4VisAttributes::Invisible);  
  G4VisAttributes *target_LH2_vis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.7));
  target_LH2_vis->SetForceWireframe(1);
  //target_LH2_vis->SetForceSolid(1);
  LH2_target_log->SetVisAttributes(target_LH2_vis);
  //LH2_target_log->SetVisAttributes(G4VisAttributes::Invisible);  
  //===========================Nouse Shielding=========================
  G4RotationMatrix *nouse_shield_rot = new G4RotationMatrix();
  nouse_shield_rot->rotateX(nouse_shield_xrot);
  nouse_shield_rot->rotateY(nouse_shield_yrot);
  nouse_shield_rot->rotateY(nouse_shield_roll_angle);
  
  G4ThreeVector nouse_shield_pos(nouse_shield_x_shift, 0*cm, nouse_shield_z_shift);
  G4Trap *nouse_shield = new G4Trap("nouse_shield", nouse_shield_height, nouse_shield_lz, nouse_shield_l_thick, nouse_shield_s_thick);
  nouse_shield_log = new G4LogicalVolume(nouse_shield, Wolfram_mat, "nouse_shield_log");
  nouse_shield_phys = new G4PVPlacement(nouse_shield_rot, nouse_shield_pos, nouse_shield_log, "nouse_shield_phys", world_log, true, 0, true);
  G4VisAttributes *nose_shield_vis = new G4VisAttributes(G4Colour(0., 0.99, 0.));
  nose_shield_vis->SetForceWireframe(1);
  nose_shield_vis->SetForceSolid(1);
  nouse_shield_log->SetVisAttributes(nose_shield_vis);
  

  //========================Beamline Shielading========================
  
  G4ThreeVector beamline_shield_pos(beamline_shiels_x_shift, 0*cm, beamline_shield_z_shift);
  
  G4Box *beamline_shield = new G4Box("beamline_shield", beamline_shield_thick/2., beamline_shield_height/2.,  beamline_shield_lz/2. );
  beamline_shield_log = new G4LogicalVolume(beamline_shield, Wolfram_mat, "beamlibe_shield_log");
  beamline_shield_phys = new G4PVPlacement(0, beamline_shield_pos, beamline_shield_log, "beamline_shield_phys", world_log, true, 0, true);
  G4VisAttributes *beamline_shield_vis = new G4VisAttributes(G4Colour(0., 0.99, 0.));
  beamline_shield_vis->SetForceWireframe(1);
  beamline_shield_vis->SetForceSolid(1);
  beamline_shield_log->SetVisAttributes(beamline_shield_vis);
  


  //===================================HRS Window==========================================

  G4RotationMatrix *HRS_window_rot = new G4RotationMatrix();
  HRS_window_rot->rotateY(90*degree - HRS_window_horiz_shift);
  
  G4Sphere *HRS_window = new G4Sphere("HRS_window", HRS_window_radius, HRS_window_radius + 0.1*mm, -HRS_window_vert_spanning/2., 
				      HRS_window_vert_spanning, 90*degree - HRS_window_horiz_spanning/2., HRS_window_horiz_spanning);
  
  HRS_window_log = new G4LogicalVolume(HRS_window, Al_mat, "HRS_window_log");
  HRS_window_phys = new G4PVPlacement(HRS_window_rot, G4ThreeVector(0., 0., 0.), HRS_window_log, "HRS_window_phys", world_log, true, 0, true);
  
   G4VisAttributes *HRS_window_vis = new G4VisAttributes(G4Colour(0.25, 0.99, 0.75));
   HRS_window_vis->SetForceSolid(1);
   HRS_window_log->SetVisAttributes(HRS_window_vis);

  //=================Calorimeter By Maxime===================
   //G4double calo_dist=110.*cm;
   //G4double calo_dist=0*cm;
   //G4double calo_angle = -14.78*deg; // This comes from run 8000
   //G4double calo_angle=(-14.8071 - 2.344)*deg; // first number comes from run 8000, the scond comes from Calo center shift
   //G4double calo_angle=0*deg; // this is just for test to see where is calo at 0 degree

  double calo_angle = dvcsGlobals::Calo_angle*rad;
  double calo_dist = dvcsGlobals::Calo_distance*mm;

  cout<<"Calo angle = "<<calo_angle<<endl;
  cout<<"Calo dist = "<<calo_dist<<endl;

  DVCSCaloConstruction* Calo;
  Calo->Construction(world_log, calo_dist, calo_angle);

  //-----------------------------------------------------------------------
  //                        Sensitive Detectors
  //-----------------------------------------------------------------------


  G4SDManager *SDman = G4SDManager::GetSDMpointer();

  //HRSwindowSD *HRS_wind_SD = new HRSwindowSD("HRS_wind_SD");
  //SDman->AddNewDetector(HRS_wind_SD);
  //HRS_window_log->SetSensitiveDetector(HRS_wind_SD);

  return world_phys;
}
