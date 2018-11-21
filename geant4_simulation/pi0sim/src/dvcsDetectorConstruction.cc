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
#include "G4Polycone.hh"
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
  world_log(0), LH2_target_log(0), target_container_log(0), Al_scat_chamber_log(0), Vac_scat_chamber_log(0), air_window_log(0), kapton_window_log(0), calo_al_window_log(0),
  nouse_shield_log(0), beamline_shield_log(0), HRS_window_log(0), pb_pipe_log(0),
  world_phys(0), LH2_target_phys(0), target_container_phys(0), Al_scat_chamber_phys(0), Vac_scat_chamber_phys(0), air_window_phys(0), kapton_window_phys(0), calo_al_window_phys(0), nouse_shield_phys(0), beamline_shield_phys(0), HRS_window_phys(0),pb_pipe_phys(0)
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
  G4double density;
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
  G4Material *Pb_mat = nist_man->FindOrBuildMaterial("G4_Pb");
  
  //===============Some quantities====================================
  G4double inch = 2.54*cm;
  G4double world_lx = 10*m;
  G4double world_ly = 10*m;
  G4double world_lz = 10*m;

  //scattering chamber
  G4double sc_or = 30.315*inch;//from A06114-03-03-0201 
  G4double sc_ir = 29.315*inch;//from A06114-03-03-0201 
  G4double sc_thickness = sc_or - sc_ir;
  G4double sc_height = 12.5*inch;//from A06114-03-03-0201 
  
  //beamline tubes
  G4double dt_thickness=0.28*inch;//from A06114-03-03-0202 
  G4double dt_or = 0.5*6.625*inch;//from A06114-03-03-0202 
  G4double dt_ir=dt_or-dt_thickness;
  G4double dt_length=150.*cm;
  G4double ut_radius=3.*cm;
  G4double ut_thickness=0.3*cm;
  G4double ut_ir=ut_radius-ut_thickness/2.;
  G4double ut_or=ut_radius+ut_thickness/2.;
  G4double ut_length=100.*cm;
  G4double dt_shift = sqrt( sc_ir*sc_ir - dt_ir*dt_ir );
  G4double ut_shift = sqrt( sc_ir*sc_ir - ut_ir*ut_ir );

  //LH2 Target
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
  
  //HRS window **Aluminum or Kapton??
  G4double kapton_thickness = 0.016*inch;//from A06114-03-03-0105 
  G4double kapton_height = 5*inch;//from A06114-03-03-0102
  G4double kapton_horiz_angle = 47.6*degree;//from A06114-03-03-0102 
  G4double kapton_horiz_shift = (90-55.6)*degree;//from A06114-03-03-0100

  //Scattering chamber calorimeter window
  G4double al_window_height = 27.9*cm;//from survey
  G4double al_window_thickness = .394*inch;//from A06114-03-03-0201
  G4double calo_al_horiz_angle = 21.57*degree;//from survey
  G4double calo_al_horiz_shift = (21.57/2+7.3)*degree;//from survey
  
  //nouse shield
  G4double nouse_shield_lz = 9.979*cm; // length along z direction
  G4double nouse_shield_height = 13.39*inch;//from A08025-03-01-0001
  G4double nouse_shield_l_thick = .84*inch; //from A08025-03-01-0001
  G4double nouse_shield_angle = 9.3*degree;//from A08025-03-01-0001
  G4double nouse_shield_s_thick =  nouse_shield_l_thick - nouse_shield_lz*tan(nouse_shield_angle);
  G4double nouse_shield_m_thick = 0.5*(nouse_shield_l_thick + nouse_shield_s_thick); 
  G4double nouse_shield_z_shift =  77.813*cm + nouse_shield_lz/2.;//from survey
  G4double nouse_shield_x_shift = -11.114*cm + nouse_shield_l_thick - 0.5*nouse_shield_m_thick;//from survey 
  G4double nouse_shield_roll_angle = 0.0*degree; //clockwise about the beam's axis
  G4double nouse_shield_xrot = 90*degree;
  G4double nouse_shield_yrot = 180*degree;
  
  //beamline shield
  G4double beamline_shield_lz = 10.236*inch;
  G4double beamline_shield_height = 7.087*inch;
  G4double beamline_shield_thick = 0.591*inch;
  G4double beamline_shield_z_shift = nouse_shield_z_shift + 0.5*nouse_shield_lz + 0.5*( beamline_shield_lz);
  G4double beamline_shield_x_shift = nouse_shield_x_shift + nouse_shield_m_thick/2 - beamline_shield_thick/2;//beam side of two shields are flush
  
  //Lead shielding pipe
  G4double pb_pipe_length = 39.435*inch;
  G4double pb_pipe_x = 0*cm;
  G4double pb_pipe_y = 0*cm;
  G4double pb_pipe_z = beamline_shield_z_shift+ 0.5*( beamline_shield_lz)+ pb_pipe_length/2+1*inch;//looks like a 1" gap in drawing
  G4double pb_pipe_ir = 3.44*inch;
  G4double pb_pipe_or = 4.44*inch;

  G4double HRS_window_horiz_spanning = 2*60*milliradian;
  G4double HRS_window_vert_spanning = 2*70*milliradian;
  G4double HRS_window_horiz_shift = dvcsGlobals::HRS_angle*rad;
  G4double HRS_window_vert_shift = 0*degree;
  G4double HRS_window_radius = 1.1086*m;

 
  //==================world====================
  
  G4Box *world_box = new G4Box("world_box", world_lx, world_ly, world_lz);
  world_log = new G4LogicalVolume(world_box, Air, "world_log");
  world_phys= new G4PVPlacement(0, G4ThreeVector(0., 0.,0. ), world_log, "world_phys", 0, 0, false);
  world_log->SetVisAttributes(G4VisAttributes::Invisible);

  
  //==================Scattering Chamber with up and down tubes and Capton/Air windows=========================

  G4ThreeVector dt_pos( 0, +dt_shift + dt_length/2.,0);
  G4RotationMatrix *dt_rot = new G4RotationMatrix();
  dt_rot->rotateX(90*degree);
  G4ThreeVector ut_pos( 0, -ut_shift - ut_length/2.,0);
  G4RotationMatrix *ut_rot = new G4RotationMatrix();
  ut_rot->rotateX(90*degree);
  
  //==========================================================================================
  //-------------------A newer Method for construction of Scattering chamber------------------
  
  G4double Al_zplane[6]={-3*sc_height,-sc_height/2-1*inch,-sc_height/2-1*inch,sc_height/2+1*inch,sc_height/2+1*inch,3*sc_height};
  G4double Al_rInner[6]={0,0,0,0,0,0};
  G4double Al_rOuter[6]={sc_or-10*cm,sc_or-10*cm,sc_or,sc_or,sc_or-10*cm,sc_or-10*cm};

  G4VSolid *Al_scat_tube = new G4Polycone("Al_scat_tube",0*degree,360*degree,6,Al_zplane,Al_rInner,Al_rOuter);
  G4VSolid *Al_down_tube = new G4Tubs("Al_down_tube", 0, dt_or, dt_length/2, 0, 360*degree);
  G4VSolid *Al_uper_tube = new G4Tubs("Al_uper_tube", 0, ut_or, ut_length/2, 0, 360*degree);
  G4VSolid *Al_intermed_chamber = new G4UnionSolid("Al_intermed_chamber", Al_scat_tube, Al_down_tube, dt_rot, dt_pos);
  G4VSolid *Al_scat_chamber = new G4UnionSolid("Al_scat_chamber", Al_intermed_chamber, Al_uper_tube, ut_rot, ut_pos);
  
  //---------Vacuum Scattering Chamber----------

  G4double Vac_zplane[6]={-3*sc_height+1*inch,-sc_height/2,-sc_height/2,sc_height/2,sc_height/2,3*sc_height-1*inch};
  G4double Vac_rInner[6]={0,0,0,0,0,0};
  G4double Vac_rOuter[6]={sc_or-11*cm,sc_or-11*cm,sc_or-1*inch,sc_or-1*inch,sc_or-11*cm,sc_or-11*cm};

  G4VSolid *Vac_scat_tube = new G4Polycone("Vac_scat_tube",0*degree,360*degree,6,Vac_zplane,Vac_rInner,Vac_rOuter);
  G4VSolid *Vac_down_tube = new G4Tubs("Vac_down_tube", 0, dt_ir, dt_length/2, 0, 360*degree);
  G4VSolid *Vac_uper_tube = new G4Tubs("Vac_uper_tube", 0, ut_ir, ut_length/2, 0, 360*degree);
  G4VSolid *Vac_intermed_chamber = new G4UnionSolid("Vac_intermed_chamber", Vac_scat_tube, Vac_down_tube, dt_rot, dt_pos);
  G4VSolid *Vac_scat_chamber = new G4UnionSolid("Vac_scat_chamber", Vac_intermed_chamber, Vac_uper_tube, ut_rot, ut_pos);

  Al_scat_chamber_log = new G4LogicalVolume(Al_scat_chamber, Al_mat, "Al_scat_chamber");

  G4RotationMatrix *al_scat_rot = new G4RotationMatrix();
  al_scat_rot->rotateX(90*degree);
  al_scat_rot->rotateY(0*degree);
  al_scat_rot->rotateZ(180*degree);
  Al_scat_chamber_phys = new G4PVPlacement(al_scat_rot, G4ThreeVector(0., 0., 0.), Al_scat_chamber_log, "Al_scat_chamber_phys", world_log, true, 0, true);

  G4RotationMatrix *vac_scat_rot = new G4RotationMatrix();
  vac_scat_rot->rotateX(0*degree);
  vac_scat_rot->rotateY(0*degree);
  vac_scat_rot->rotateZ(0*degree);
  Vac_scat_chamber_log = new G4LogicalVolume(Vac_scat_chamber, Vac_mat, "Vac_scat_chamber");
  Vac_scat_chamber_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Vac_scat_chamber_log,"Vac_scat_chamber_phys", Al_scat_chamber_log, true, 0, true);

  G4VisAttributes *scat_chamber_vis = new G4VisAttributes(G4Colour(0.1, 0.2, 0.9));
  scat_chamber_vis->SetForceWireframe(1);
  // scat_chamber_vis->SetForceSolid(0);
  Al_scat_chamber_log->SetVisAttributes(scat_chamber_vis);
  //  Al_scat_chamber_log->SetVisAttributes(G4VisAttributes::Invisible);

   G4VisAttributes *vac_scat_chamber_vis = new G4VisAttributes(G4Colour(0.99, 0.99, 0.02));
   //   vac_scat_chamber_vis->SetForceWireframe(1);
   vac_scat_chamber_vis->SetForceSolid(1);
   Vac_scat_chamber_log->SetVisAttributes(vac_scat_chamber_vis);
   //  Vac_scat_chamber_log->SetVisAttributes(G4VisAttributes::Invisible);

  //===============================Air Window=============================

  G4Tubs *air_window = new G4Tubs("air_window",sc_ir, sc_ir + sc_thickness, kapton_height, -kapton_horiz_angle/2, kapton_horiz_angle);
  air_window_log = new G4LogicalVolume(air_window, Air, "air_window_log");
  G4RotationMatrix *air_rot = new G4RotationMatrix();
  air_rot->rotateX(0*degree);
  air_rot->rotateZ(0*degree);
  air_rot->rotateZ(-kapton_horiz_shift-90*degree);
  air_window_phys = new G4PVPlacement(air_rot, G4ThreeVector(0., 0., 0.), air_window_log, "air_window_phys", Al_scat_chamber_log, true, 0, true);
  G4VisAttributes *air_window_vis = new G4VisAttributes(G4Colour(0.01, 0.99, 0.99));
  air_window_vis->SetForceWireframe(0);
  air_window_vis->SetForceSolid(1);
  air_window_log->SetVisAttributes(air_window_vis);
  //air_window_log->SetVisAttributes(G4VisAttributes::Invisible);

  //=============================Capton Window============================

  G4Tubs *kapton_window = new G4Tubs("kapton_window",sc_ir, sc_ir + kapton_thickness, kapton_height, -kapton_horiz_angle/2, kapton_horiz_angle);
  kapton_window_log = new G4LogicalVolume(kapton_window, Kapton_mat, "kapton_window_log");
  G4RotationMatrix *kapton_rot = new G4RotationMatrix();
  //kapton_rot->rotateX(90*degree);
  //kapton_rot->rotateY(90*degree);
  //kapton_rot->rotateY(-kapton_horiz_shift);
  kapton_window_phys = new G4PVPlacement(kapton_rot, G4ThreeVector(0., 0., 0.), kapton_window_log, "kapton_window_phys", air_window_log,true, 0, true);
  G4VisAttributes *kapton_window_vis = new G4VisAttributes(G4Colour(0.99, 0.02, 0.99));
  kapton_window_vis->SetForceWireframe(0);
  kapton_window_vis->SetForceSolid(1);
  kapton_window_log->SetVisAttributes(kapton_window_vis);

  //===============================Calo Alum Window=============================

  G4Tubs *calo_al_window = new G4Tubs("calo_al_window", sc_ir + al_window_thickness,sc_or, al_window_height/2, -calo_al_horiz_angle/2, calo_al_horiz_angle);
  calo_al_window_log = new G4LogicalVolume(calo_al_window, Air, "calo_al_window_log");
  G4RotationMatrix *calo_al_rot = new G4RotationMatrix();
  calo_al_rot->rotateX(0*degree);
  calo_al_rot->rotateY(0*degree);
  calo_al_rot->rotateZ(calo_al_horiz_shift - 90*degree);
  calo_al_window_phys = new G4PVPlacement(calo_al_rot, G4ThreeVector(0., 0., 0.), calo_al_window_log, "calo_al_window_phys", Al_scat_chamber_log, true, 0, true);
  G4VisAttributes *calo_al_window_vis = new G4VisAttributes(G4Colour(0.01, 0.99, 0.99));
  // calo_al_window_vis->SetForceWireframe(0);
  calo_al_window_vis->SetForceSolid(1);
  calo_al_window_log->SetVisAttributes(calo_al_window_vis);
  //calo_al_window_log->SetVisAttributes(G4VisAttributes::Invisible);

  //====================target====================
  
  //=======container=======
  G4Tubs *target_container_tube = new G4Tubs("target_container_tube", 0, target_container_radius, target_tube_length/2., 0, 360.*degree);
  G4Sphere *target_container_endcup = new G4Sphere("target_container_endcup", 0, target_container_radius, 0.*degree, 360.*degree, 0., 90*degree);
  G4UnionSolid *target_container_solid = new G4UnionSolid("target_container_solid", target_container_tube, target_container_endcup, 0, G4ThreeVector(0., 0., 0.999*(target_tube_length/2.)));
  target_container_log = new G4LogicalVolume(target_container_solid, Al_mat, "target_container_log");

  G4RotationMatrix *target_rot = new G4RotationMatrix();
  target_rot->rotateX(90*degree);
  target_rot->rotateY(0*degree);
  target_rot->rotateZ(0*degree);

  target_container_phys = new G4PVPlacement(target_rot, G4ThreeVector(0., 0., target_container_tube_center), 
  target_container_log, "target_container_phys", Vac_scat_chamber_log, true, 0, true);

  //=======LH2=======
  G4Tubs *target_LH2_tube = new G4Tubs("target_LH2_tube", 0, target_LH2_radius, target_LH2_tube_length/2., 0, 360.*degree);
  G4Sphere *target_LH2_endcup = new G4Sphere("target_LH2_endcup", 0, target_LH2_radius, 0., 360.*degree, 0., 90.*degree);

  G4UnionSolid *target_LH2_solid = new G4UnionSolid("target_LH2_solid", target_LH2_tube, target_LH2_endcup, 0,G4ThreeVector(0, 0, 0.999*(target_LH2_tube_length/2.)));
  LH2_target_log = new G4LogicalVolume(target_LH2_solid, LH2_mat, "LH2_target_log");
  LH2_target_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., target_LH2_tube_center), LH2_target_log, "LH2_target_phys", target_container_log,true, 0, true);

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
  
  //========================Beamline Shielding========================
  
  G4ThreeVector beamline_shield_pos(beamline_shield_x_shift, 0*cm, beamline_shield_z_shift);
  G4Box *beamline_shield = new G4Box("beamline_shield", beamline_shield_thick/2., beamline_shield_height/2.,  beamline_shield_lz/2. );
  beamline_shield_log = new G4LogicalVolume(beamline_shield, Wolfram_mat, "beamlibe_shield_log");
  beamline_shield_phys = new G4PVPlacement(0, beamline_shield_pos, beamline_shield_log, "beamline_shield_phys", world_log, true, 0, true);

  G4VisAttributes *beamline_shield_vis = new G4VisAttributes(G4Colour(0., 0.99, 0.));
  beamline_shield_vis->SetForceWireframe(1);
  beamline_shield_vis->SetForceSolid(1);
  beamline_shield_log->SetVisAttributes(beamline_shield_vis);

  //========================Lead Pipe Beamline Shielding========================
  
  G4ThreeVector pb_pipe_pos(pb_pipe_x, pb_pipe_y, pb_pipe_z);
  G4Tubs *pb_pipe = new G4Tubs("pb_pipe", pb_pipe_ir, pb_pipe_or, pb_pipe_length/2,90*degree,180*degree );
  pb_pipe_log = new G4LogicalVolume(pb_pipe, Pb_mat, "pb_pipe_log");
  pb_pipe_phys = new G4PVPlacement(0, pb_pipe_pos, pb_pipe_log, "pb_pipe_phys", world_log, true, 0, true);

  G4VisAttributes *pb_pipe_vis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  pb_pipe_vis->SetForceWireframe(1);
  pb_pipe_vis->SetForceSolid(1);
  pb_pipe_log->SetVisAttributes(pb_pipe_vis);
  
  //===================================HRS Window==========================================

  G4RotationMatrix *HRS_window_rot = new G4RotationMatrix();
  HRS_window_rot->rotateY(90*degree - HRS_window_horiz_shift);

  G4Sphere *HRS_window = new G4Sphere("HRS_window", HRS_window_radius, HRS_window_radius + 0.1*mm, -HRS_window_vert_spanning/2.,HRS_window_vert_spanning, 90*degree - HRS_window_horiz_spanning/2., HRS_window_horiz_spanning);
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
