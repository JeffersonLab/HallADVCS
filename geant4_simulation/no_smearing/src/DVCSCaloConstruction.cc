// $Id: DVCSCaloConstruction.cc,v 1.9 2006/06/29 17:47:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-03-patch-01 $
//

#include "DVCSCaloConstruction.hh"
//#include "DVCSField.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "globals.hh"
//#include "DVCSCaloSD.hh"
#include "G4SDManager.hh"
#include "G4Orb.hh"
#include "G4VisAttributes.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "DVCSCaloGeom.hh"
#include "DVCSCalo.hh"
#include "TDVCSDB.h"
#include "dvcsGlobals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh""

DVCSCaloConstruction::DVCSCaloConstruction()
 : calorimeterBlock_log(0),  Laiton_frame_log(0), G10_log(0), Screw_log(0), Tyvek_log(0), Tedlar_log(0), OutAl_log(0), PVC_log(0), InAl_log(0),
     CaloBlock_phys(0), Laiton_frame_phys(0), G10_phys(0), Screw_phys(0), Tyvek_phys(0), Tedlar_phys(0), OutAl_phys(0), PVC_phys(0) , InAl_phys(0)
{;}

DVCSCaloConstruction::~DVCSCaloConstruction()
{
}

void DVCSCaloConstruction::Construction(G4LogicalVolume* Log, G4double calo_dist, G4double calo_angle)
{
  int nbcol=13;
  int nbrow=16;
  
  //------------------------------------------------------ materials
  
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4int ncomponents;

  G4Material* Al = 
  new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  G4Material* LH2 = 
  new G4Material("LH2",  z= 1., a=1.00794*g/mole, density= 0.0723*g/cm3);

  G4Element* elPb = new G4Element("Lead","Pb",z=82.,a=207.2*g/mole);
  G4Element* elF = new G4Element("Fluor","F",z=9.,a=18.998*g/mole);

  G4Material* PbF2 = new G4Material("PbF2",density=7.77*g/cm3,ncomponents=2);
  PbF2->AddElement(elPb,1);
  PbF2->AddElement(elF, 2);

  G4Element* elCu= new G4Element("Copper","Cu",z=29.,a=63.546*g/mole);
  G4Element* elZn= new G4Element("Stain","Zn",z=30.,a=65.38*g/mole);

  G4Material* Laiton= new G4Material("Laiton", density=8.564*g/cm3, ncomponents=2);
  Laiton->AddElement(elCu,4);
  Laiton->AddElement(elZn,1);

  G4Element* elC= new G4Element("Carbon","C",z=6.,a=12.01*g/mole);
  G4Element* elH= new G4Element("Hydrogen","H",z=1.,a=1.008*g/mole);
  G4Element* elCl= new G4Element("Chlorine","Cl",z=17.,a=35.453*g/mole);

  G4Material* G10= new G4Material("G10", z=8.455, a=16.91*g/mole, density=1.90*g/cm3);

  G4Material* CH2= new G4Material("CH2", density=0.89*g/cm3, ncomponents=2);
  CH2->AddElement(elC,1);
  CH2->AddElement(elH,2);

  G4Element* elFe= new G4Element("Iron","Fe",z=26.,a=55.845*g/mole);
  G4Element* elCr= new G4Element("Chromium","Cr",z=24.,a=51.996*g/mole);

  G4Material* Steel= new G4Material("Steel", density=7.68*g/cm3, ncomponents=3);
  Steel->AddElement(elFe,87);
  Steel->AddElement(elCr,11);
  Steel->AddElement(elC,2);

  G4Material* Tyvek= new G4Material("Tyvek", density=0.935*g/cm3, ncomponents=2);
  Tyvek->AddElement(elC,1);
  Tyvek->AddElement(elH,2);
 
  //Tedlar is composed of polycarbonates and PVC. I have assumed that there is one CH2-CH2 for one CHCl-CH2
  G4Material* Tedlar= new G4Material("Tedlar", density=1.2*g/cm3, ncomponents=3);
  Tedlar->AddElement(elC,4);
  Tedlar->AddElement(elH,7);
  Tedlar->AddElement(elCl,1);

  G4Material* PVC= new G4Material("PVC", density=1.30/2.*g/cm3, ncomponents=3);
  PVC->AddElement(elC,2);
  PVC->AddElement(elH,3);
  PVC->AddElement(elCl,1);

  G4NistManager* man = G4NistManager::Instance();
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

  //------------------------------------------------------ volumes
 
  DVCSCaloGeom* CaloGeom;

  //------------------------------ CuZn PMT carrier

  G4double flas_length=0.2*cm;

  G4LogicalVolume* Laiton_frame_log=CaloGeom->PMCarrierlog(Laiton);
  G4VisAttributes *at3=new G4VisAttributes(G4Color(0,0,100,1));
  Laiton_frame_log->SetVisAttributes(at3);

  //------------------------------ G10

  G4LogicalVolume* G10_log= CaloGeom->G10log(G10);
  G4VisAttributes *at1=new G4VisAttributes(G4Color(0,100,0,1));
  G10_log->SetVisAttributes(at1); 

  //------------------------------ a screw

  G4LogicalVolume* Screw_log=CaloGeom->Screwlog(Steel);
  G4VisAttributes *at=new G4VisAttributes(G4Color(100,0,0,1));
  Screw_log->SetVisAttributes(at);

  //------------------------------ a calorimeter block

  G4double block_length=18.6*cm;
  G4double block_width =3.*cm;

  G4LogicalVolume* CaloBlock_log = CaloGeom->CaloBlocklog(PbF2);

  //------------------------------ Tyvek paper (white)
  
  G4LogicalVolume* Tyvek_log= CaloGeom->Tyveklog(Tyvek);
  G4VisAttributes *atTy=new G4VisAttributes(G4Color(100,100,0,1));
  Tyvek_log->SetVisAttributes(atTy);
  
  //------------------------------ Tedlar paper (black)

  G4LogicalVolume* Tedlar_log= CaloGeom->Tedlarlog(Tedlar);
  G4VisAttributes *atTed=new G4VisAttributes(G4Color(100,0,100,1));
  Tedlar_log->SetVisAttributes(atTed);

  //------------------------------ Small Slim 
  /*
  G4LogicalVolume* SShim_log= CaloGeom->SShimlog(Laiton);
  SShim_log->SetVisAttributes(at1);
  */

  //------------------------------ Big Shim
  /*
  G4LogicalVolume* BShim_log= CaloGeom->BShimlog(Laiton);
  BShim_log->SetVisAttributes(at1);
  */
  //------------------------------ the calorimeter

  G4double calo_x = 2.002*double(nbcol)*block_width/2.;
  G4double calo_y = 1.2*double(nbrow)*block_width/2.;
  G4double calo_z = 55*cm;
  G4Box* Calo_box
    = new G4Box("Calo_box",calo_x,calo_y,calo_z);
  G4LogicalVolume* Calo_log = new G4LogicalVolume(Calo_box,Air,"Calo_log",0,0,0);
  
  Calo_log-> SetVisAttributes (at1);
  //------------------------------- now we put the blocks in the calorimeter
 
  //G4int nSShim=0;
  //G4int nBShim=0;
  G4double z_block=block_length/2.-calo_z+3*cm;
  G4double x_global=calo_x/1.85;

  /*
  for(G4int col=0;col<nbcol;col++)
    {
    for(G4int row=0;row<nbrow;row++)
      {
	x_block[col*nbrow+row] = -230.41 + col*3.014*cm;
	y_block[col*nbrow+row] = -229.39 + row*3.014*cm;
      }
    }
  */

  
  int run_number = dvcsGlobals::run_number;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");
  Float_t *x_blockDB = db->GetEntry_f("CALO_geom_X",run_number);
  Float_t *y_blockDB = db->GetEntry_f("CALO_geom_Y",run_number);

  for(G4int col=0;col<nbcol;col++)
    {
    for(G4int row=0;row<nbrow;row++)
      {
	x_block[col*nbrow+row] = x_blockDB[col*nbrow+row]*10.0;
	y_block[col*nbrow+row] = y_blockDB[col*nbrow+row]*10.0;
      }
    }
  

  for(G4int col=0;col<nbcol;col++)
    {
    for(G4int row=0;row<nbrow;row++)
      {
	G4double x, y, z_laiton, z_G10;
	G4double adjust_y=0;
	G4double adjust_x=0;
	
	//adjust_x and adjust_y avoid some overlaps in the building of the calorimeter. 
	if(run_number>=10000 && run_number<=12199){
	  //Fall 2014
	  if(col*nbrow+row==12) adjust_x=-0.007*cm;
	  if(col*nbrow+row==28) {adjust_x=-0.002*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==29) adjust_x=0.0001*cm;
	  if(col*nbrow+row==31) adjust_y=0.0005*cm;
	  if(col*nbrow+row==44) adjust_x=-0.001*cm;
	  if(col*nbrow+row==45) {adjust_x=0.002*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==60) adjust_x=-0.0055*cm;
	  if(col*nbrow+row==71) adjust_x=0.005*cm;
	  if(col*nbrow+row==75) adjust_x=0.0045*cm;
	  if(col*nbrow+row==76) {adjust_x=-0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==91) adjust_x=0.0005*cm;
	  if(col*nbrow+row==92) adjust_x=-0.0005*cm;
	  if(col*nbrow+row==93) {adjust_x=0.0015*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==94) adjust_x=-0.001*cm;
	  if(col*nbrow+row==95) adjust_y=0.002*cm;
	  if(col*nbrow+row==109) adjust_x=0.005*cm;
	  if(col*nbrow+row==111) adjust_x=0.0025*cm;
	  if(col*nbrow+row==124) adjust_x=0.002*cm;
	  if(col*nbrow+row==127) adjust_x=0.001*cm;
	  if(col*nbrow+row==145) adjust_y=0.002*cm;
	  if(col*nbrow+row==159) adjust_y=0.005*cm;
	  if(col*nbrow+row==161) adjust_y=0.001*cm;
	  if(col*nbrow+row==171) adjust_x=0.001*cm;
	  if(col*nbrow+row==173) adjust_y=-0.011*cm;
	  if(col*nbrow+row==174) adjust_y=0.011*cm;
	  if(col*nbrow+row==175) adjust_y=0.0005*cm;
	  if(col*nbrow+row==189) adjust_x=0.002*cm;
	  if(col*nbrow+row==190) adjust_x=0.002*cm;
	  if(col*nbrow+row==207) adjust_x=0.001*cm;
	}

	else if(run_number>=12200 && run_number<=13799){
	  //Spring 2016
	  if(col*nbrow+row==26) adjust_y=-0.001*cm;
	  if(col*nbrow+row==28) adjust_y=-0.001*cm;
	  if(col*nbrow+row==31) adjust_y=0.001*cm;
	  if(col*nbrow+row==43) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==45) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==57) adjust_x=0.002*cm;
	  if(col*nbrow+row==76) adjust_y=-0.001*cm;
	  if(col*nbrow+row==83) adjust_y=0.001*cm;
	  if(col*nbrow+row==93) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==98) {adjust_x=0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==110) adjust_y=0.001*cm;
	  if(col*nbrow+row==115) adjust_y=0.001*cm;
	  if(col*nbrow+row==125) {adjust_x=0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==130) {adjust_x=0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==140) adjust_y=0.001*cm;
	  if(col*nbrow+row==142) {adjust_x=-0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==155) {adjust_x=0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==157) {adjust_x=0.002*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==159) adjust_y=0.002*cm;
	  if(col*nbrow+row==165) adjust_x=0.001*cm;
	}

	else if(run_number>=13800){
	  //Fall 2016
	  if(col*nbrow+row==10) adjust_x=-0.001*cm;
	  if(col*nbrow+row==25) adjust_x=0.001*cm;
	  if(col*nbrow+row==27) adjust_y=-0.002*cm;
	  if(col*nbrow+row==28) adjust_y=0.002*cm;
	  if(col*nbrow+row==44) {adjust_x=-0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==46) {adjust_x=-0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==57) adjust_x=0.002*cm;
	  if(col*nbrow+row==61) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==63) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==65) {adjust_x=-0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==68) adjust_y=-0.001*cm;
	  if(col*nbrow+row==69) {adjust_x=-0.0015*cm; adjust_y=0.0015*cm;}
	  if(col*nbrow+row==78) {adjust_x=-0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==82) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==83) {adjust_x=0.001*cm; adjust_y=0.002*cm;}
	  if(col*nbrow+row==85) {adjust_x=-0.0015*cm; adjust_y=0.0015*cm;}
	  if(col*nbrow+row==95) {adjust_x=0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==98) {adjust_x=0.001*cm; adjust_y=-0.002*cm;}
	  if(col*nbrow+row==100) {adjust_x=0.0015*cm; adjust_y=-0.0015*cm;}
	  if(col*nbrow+row==127) adjust_y=0.001*cm;
	  if(col*nbrow+row==131) {adjust_x=-0.001*cm; adjust_y=0.001*cm;}
	  if(col*nbrow+row==142) {adjust_x=-0.002*cm; adjust_y=0.002*cm;}
	  if(col*nbrow+row==146) {adjust_x=0.001*cm; adjust_y=-0.001*cm;}
	  if(col*nbrow+row==157) {adjust_x=0.002*cm; adjust_y=-0.002*cm;}
	  if(col*nbrow+row==163) {adjust_x=-0.002*cm; adjust_y=0.002*cm;}
	  if(col*nbrow+row==178) {adjust_x=0.002*cm; adjust_y=-0.002*cm;}
	  if(col*nbrow+row==195) adjust_x=0.001*cm;
	}
	
	/*
	//2010
	if (col*nbrow+row==6) adjust_y=-0.003*cm;
	if (col*nbrow+row==7) adjust_y=0.0025*cm;
	if (col*nbrow+row==29){ 
	  adjust_x=0.001*cm; adjust_y=0.001*cm;
	}
	if (col*nbrow+row==51) adjust_y=-0.0035*cm;
	if (col*nbrow+row==52) adjust_y=-0.0025*cm;
	if (col*nbrow+row==57){ 
	  adjust_x=0.001*cm; adjust_y=0.001*cm;
	}
	if (col*nbrow+row==59){ 
	  adjust_x=0.001*cm; adjust_y=0.001*cm;
	}
	if (col*nbrow+row==68){ 
	  adjust_x=0.001*cm; adjust_y=-0.001*cm;
	}
	if (col*nbrow+row==70){ 
	  adjust_x=0.001*cm; adjust_y=-0.001*cm;
	}
	if (col*nbrow+row==100){ 
	  adjust_x=0.001*cm; adjust_y=0.001*cm;
	}
	if (col*nbrow+row==109){ 
	  adjust_x=0.001*cm; adjust_y=0.001*cm;
	}
	if (col*nbrow+row==155) adjust_x=0.001*cm;
	if (col*nbrow+row==53) adjust_y=0.0005*cm;
	if (col*nbrow+row==124) adjust_y=0.003*cm;
	if (col*nbrow+row==80) adjust_x=0.0015*cm;
	if (col*nbrow+row==139) adjust_x=0.005*cm;
	if (col*nbrow+row==68) adjust_x=-0.001*cm;
	if (col*nbrow+row==189) adjust_x=0.003*cm;
	if (col*nbrow+row==191){ 
	  adjust_x=0.001*cm; adjust_y=0.001*cm;
	}
	if (col*nbrow+row==195) adjust_x=0.003*cm;
	*/

	x=(x_block[col*nbrow+row]/10)*cm+adjust_x*cm;
	y=(y_block[col*nbrow+row]/10)*cm+adjust_y*cm;
	
	
	//=============================*******-----------))))))))::::::((((((((--------*********================================
	// adjust_x = 0;
	// adjust_y = 0;
	
	// x_block[col*nbrow+row] = -230.41 + col*3.2*cm;
	// y_block[col*nbrow+row] = -229.39 + row*3.2*cm;

	// x = x_block[col*nbrow+row];
	// y = y_block[col*nbrow+row];
	//=============================*******-----------))))))))::::::((((((((--------*********================================
	

	z_laiton=-0.5*cm-1*flas_length-0.999*block_length/2.+z_block;
	z_G10=-flas_length/2-0.999*block_length/2.+z_block;
	
	/*
	if (col*nbrow+row==SShim[nSShim]){
	  nSShim=nSShim+1;
	  CaloGeom->SShim(SShim_log, Calo_log, G4ThreeVector(
	  (x_block[col*nbrow+row]/10)/2*cm+(x_block[col*nbrow+row+1]/10)/2*cm+x_global,
	  (y_block[col*nbrow+row]/10)/2*cm+(y_block[col*nbrow+row+1]/10)/2*cm,z_block), col*nbrow+row);
	}

	if (col*nbrow+row==BShim[nBShim]){
	  nBShim=nBShim+1;
	  CaloGeom->BShim(BShim_log, Calo_log, G4ThreeVector(
	  (x_block[col*nbrow+row]/10)/2*cm+(x_block[col*nbrow+row+1]/10)/2*cm+x_global,
	  (y_block[col*nbrow+row]/10)/2*cm+(y_block[col*nbrow+row+1]/10)/2*cm,z_block), col*nbrow+row);
	}
	*/

	CaloGeom->G10(G10_log, Calo_log, G4ThreeVector(x+x_global,y,z_G10), col*nbrow+row);
	
	CaloGeom->Screw(Screw_log, Calo_log, G4ThreeVector(x+x_global,y,z_laiton), col*nbrow+row);
	
	CaloGeom->PMCarrier(Laiton_frame_log, Calo_log, G4ThreeVector(x+x_global,y,z_laiton), col*nbrow+row);
	
	CaloGeom->Tedlar(Tedlar_log, Calo_log, G4ThreeVector(x+x_global,y,z_block), col*nbrow+row);

	CaloGeom->Tyvek(Tyvek_log, Calo_log, G4ThreeVector(x+x_global,y,z_block), col*nbrow+row);

	CaloGeom->CaloBlock(CaloBlock_log, Calo_log, G4ThreeVector(x+x_global,y,z_block), col*nbrow+row);
	
      }
    }

  //---------------------------- Aluminium Frame for the blocks
  
  G4double Side_width=2.*cm;
  G4double basis_width=39.9*cm;
  
  G4LogicalVolume* Frame_log = CaloGeom->Alframelog(Al);
  CaloGeom->Alframe(Frame_log, Calo_log,G4ThreeVector(-4.5*cm+(basis_width+Side_width)/2.+x_global,-0.3*cm,-0.1*cm+z_block));

  //----------------------------- CH2 cover for the Aluminium frame

  G4LogicalVolume* CHcover_log = CaloGeom->CHcoverlog(CH2);
  CaloGeom->CHcover(CHcover_log, Calo_log, G4ThreeVector(-4.5*cm+x_global,-0.3*cm,-0.1*cm+z_block));
 
 
  //------------------------------- Black Box
  //It is made with 3 boxes. Two are in Aluminium. Between these two Al boxes,  there is a PVC box. 

  G4LogicalVolume* OutAl_log = CaloGeom->BB_OutAllog(Al);
  CaloGeom->BB_OutAl(OutAl_log, Calo_log, G4ThreeVector(-4.5*cm+x_global,0,z_block));

  G4LogicalVolume* PVC_log = CaloGeom->BB_PVClog(PVC);
  CaloGeom->BB_PVC(PVC_log, Calo_log, G4ThreeVector(-4.5*cm+x_global,0,z_block));


  G4LogicalVolume* InAl_log = CaloGeom->BB_InAllog(PVC);
  CaloGeom->BB_InAl(InAl_log, Calo_log, G4ThreeVector(-4.5*cm+x_global,0,z_block));

   //------------------------------------- CHshield

   G4LogicalVolume* CHshield_log=CaloGeom->CHshieldlog(CH2);
   CaloGeom->CHshield(CHshield_log, Calo_log, G4ThreeVector(-4.5*cm+x_global,0,z_block));

   //////////////////// Now we place the Calo in the hall ////////////////////
   
   G4RotationMatrix rm;
   rm.rotateY(calo_angle);
   // y pos is not 0 it is -0.0673 which comes from survey and is placed in db
   G4VPhysicalVolume* Calo_phys = new G4PVPlacement(G4Transform3D(rm,G4ThreeVector((calo_dist+calo_z-3*cm)*sin(calo_angle)-(calo_x/1.85+0*4.5*cm)*cos(calo_angle), -0.0673*cm, (calo_dist+calo_z - 3*cm)*cos(calo_angle)+(calo_x/1.85+0*4.5*cm)*sin(calo_angle))),Calo_log,"Calo",Log,false,0); 

   

  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------

  // G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // DVCSCaloSD* aCaloSD = new DVCSCaloSD( "CaloSD" );
  // SDman->AddNewDetector( aCaloSD );
  // CaloBlock_log->SetSensitiveDetector( aCaloSD );

  return;
}

