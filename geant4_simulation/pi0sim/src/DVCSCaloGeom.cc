// $Id: DVCSDetectorConstruction.cc,v 1.9 2006/06/29 17:47:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-03-patch-01 $
//

#include "DVCSCaloGeom.hh"
//#include "DVCSField.hh"
#include <TString.h>
#include "G4Trap.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "globals.hh"
//#include "DVCSCaloSD.hh"
#include "G4SDManager.hh"
#include "G4Orb.hh"
#include "G4VisAttributes.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

DVCSCaloGeom::DVCSCaloGeom()
  :  SShim_log(0), BShim_log(0), Screw_log(0), Tedlar_log(0), Tyvek_log(0), G10_log(0), PMCarrier_log(0), CHcover_log(0), CaloBlock_log(0), Alframe_log(0), CHshield_log(0), CHshield_phys(0), Alframe_phys(0), CaloBlock_phys(0), CHcover_phys(0), PMCarrier_phys(0), G10_phys(0), Tyvek_phys(0), Tedlar_phys(0), Screw_phys(0), BShim_phys(0), SShim_phys(0)
{;}

DVCSCaloGeom::~DVCSCaloGeom()
{
}

// Methods to create SmallShim /////////////////////////////


G4LogicalVolume* DVCSCaloGeom::SShimlog(G4Material* mat){

  G4double SShim_length=14.6*cm;
  G4double SShim_width=2*cm;
  G4double SShim_height=0.0254*cm;


  G4Box* SShim_box = new G4Box("SShim_box",SShim_width/2.,SShim_height/2.,SShim_length/2); 
  
  G4LogicalVolume* SShim_log= new G4LogicalVolume(SShim_box, mat,"SShim_log",0,0,0);

  return SShim_log;
}

G4VPhysicalVolume* DVCSCaloGeom::SShim(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){
  
  G4VPhysicalVolume* SShim_phys = new G4PVPlacement(0, Vect,
						    Log1,Form("SmallShim_%d",n),Log2,false,0, surf_check);
  
  return SShim_phys;
}
				 
// Method to create Big Shims /////////////////////////////////////

G4LogicalVolume* DVCSCaloGeom::BShimlog(G4Material* mat){

  G4double BShim_length=14.6*cm;
  G4double BShim_width=2*cm;
  G4double BShim_height=0.0381*cm;


  G4Box* BShim_box = new G4Box("BShim_box",BShim_width/2.,BShim_height/2.,BShim_length/2);
  
  G4LogicalVolume* BShim_log= new G4LogicalVolume(BShim_box, mat,"SShim_log",0,0,0);

  return BShim_log;
}

G4VPhysicalVolume* DVCSCaloGeom::BShim(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){
  
  G4VPhysicalVolume* BShim_phys = new G4PVPlacement(0, Vect,
						    Log1,Form("BigShim_%d",n),Log2,false,0, surf_check );
  
  return BShim_phys;
}

//Methods to create a screw //////////////////////


G4LogicalVolume* DVCSCaloGeom::Screwlog(G4Material* mat){

   G4double screw_diameter=0.25*cm;
   G4double screw_length=1*cm;

   G4Tubs *Screw_tub=new G4Tubs("Screw_tub", 0, 0.999*screw_diameter/2., screw_length/2,0,360*degree);

  G4LogicalVolume* Screw_log= new G4LogicalVolume(Screw_tub, mat,"Screw_log",0,0,0);

  return Screw_log;
}

void DVCSCaloGeom::Screw(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){

  G4double x_screw, y_screw;
  
  for(G4int screw_1=0;screw_1<2;screw_1++)
	  {
	    for (G4int screw_2=0;screw_2<2;screw_2++)
	      {
		x_screw=(2*screw_1-1)*1.2*cm;
		y_screw=(2*screw_2-1)*1.2*cm;
		G4VPhysicalVolume* Screw_phys = new G4PVPlacement(0,G4ThreeVector(x_screw,y_screw,0)+Vect,Log1,Form("Screw_%d",n),Log2,false,0, surf_check );
	      }
	  }
  
  return;
}

// Methods to create Tedlar ////////////////////////////////////////

G4LogicalVolume* DVCSCaloGeom::Tedlarlog(G4Material* mat){

  
  G4double clin_length=0.999*18.6*cm+1.*0.2*cm+0.5*cm;
  G4double flas_length=0.2*cm;
  G4double flas_height=3.*cm;
  G4double flas_width=3.*cm;

  G4Box* Tedlarout_box = new G4Box("Tyvekout_box",1.0044*flas_width/2.,1.0044*flas_height/2.,(clin_length+flas_length/2)/2);

  G4Box* Tedlarin_box = new G4Box("Tyvekin_box",1.0034*flas_width/2.,1.0034*flas_height/2.,(clin_length+flas_length/2)/1.999);

  G4SubtractionSolid* Tedlar_box= new G4SubtractionSolid("Tedlar_box", Tedlarout_box, Tedlarin_box);
  

  G4LogicalVolume* Tedlar_log= new G4LogicalVolume(Tedlar_box, mat,"Tedlar_log",0,0,0);

  return Tedlar_log;
}

void DVCSCaloGeom::Tedlar(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){
  
  G4double flas_length=0.2*cm;

  G4VPhysicalVolume* Tedlar_phys = new G4PVPlacement(0, Vect+G4ThreeVector(0,0,-0.5*cm+flas_length/2), Log1,Form("Tedlar_%d",n),Log2,false,0, surf_check);
  
  return;
}
				 
// Methods to create Tyvek paper ////////////////////////////////

G4LogicalVolume* DVCSCaloGeom::Tyveklog(G4Material* mat){

  G4double block_width =3.*cm;
  G4double block_height=3.*cm;
  G4double block_length=18.6*cm;
  
  G4Box* Tyvekout_box = new G4Box("Tyvekout_box",0.99995*block_width/2,0.99995*block_height/2,block_length/2.);

  G4Box* Tyvekin_box = new G4Box("Tyvekin_box",0.99905*block_width/2,0.99905*block_height/2,block_length/2.*1.0001);

  G4SubtractionSolid* Tyvek_box= new G4SubtractionSolid("Tyvek_box", Tyvekout_box, Tyvekin_box);

  G4LogicalVolume* Tyvek_log= new G4LogicalVolume(Tyvek_box, mat,"Tyvek_log",0,0,0);

  return Tyvek_log;
}

void DVCSCaloGeom::Tyvek(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){

  G4VPhysicalVolume* Tyvek_phys = new G4PVPlacement(0, Vect, Log1,Form("Tyvek_%d",n),Log2,false,0, surf_check);
  
  return;
}

// Methods to create G10 //////////////////////////////////////////////

		
G4LogicalVolume* DVCSCaloGeom::G10log(G4Material* mat){
  
  G4double center_diameter=2.8*cm;
  G4double flas_length=0.2*cm;
  G4double flas_height=3.*cm;
  G4double flas_width=3.*cm;


  G4Box* preG10_box= new G4Box("Flasque_box", 0.999*flas_width/2, 0.999*flas_height/2, flas_length/2);

  G4Tubs *G10_entrance_tub=new G4Tubs("G10_entrance_tub", 0, center_diameter/2., flas_length/1.999,0,360*degree);

  G4SubtractionSolid* G10_box= new G4SubtractionSolid("G10final_box", preG10_box, G10_entrance_tub);

  G4LogicalVolume* G10_log= new G4LogicalVolume(G10_box, mat,"G10_log",0,0,0);

  return G10_log;
}

G4VPhysicalVolume* DVCSCaloGeom::G10(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){
  
  G4VPhysicalVolume* G10_phys = new G4PVPlacement(0, Vect,
						  Log1,Form("G10_%d",n),Log2,false,0, surf_check);
  
  return G10_phys;
}

// Methods to create PMTcarriers in laiton ////////////////////////////

G4LogicalVolume* DVCSCaloGeom::PMCarrierlog(G4Material* mat){

  
  G4double clin_length=0.999*18.6*cm+1.*0.2*cm+0.5*cm;
  G4double clin_height=2.6*cm;
  G4double clin_width=0.005*cm;
  
 G4Box* Clinquant_box= new G4Box("clinquant", clin_width/2, clin_height/2, clin_length/2);

 G4double flas_length=0.2*cm;
 G4double flas_height=3.*cm;
 G4double flas_width=3.*cm;

 G4Box* laitonsquare_box= new G4Box("Flasque_box", flas_width/2, flas_height/2, flas_length/2);

 G4double screw_diameter=0.25*cm;
 G4double center_diameter=2.8*cm;

 G4Tubs *Block_entrance_tub=new G4Tubs("Block_entrance_tub", 0, center_diameter/2., flas_length,0,360*degree);

 G4Tubs *Screwhole_tub=new G4Tubs("Screwhole_tub", 0, screw_diameter/2, flas_length/1.999,0,360*degree);

 G4SubtractionSolid* Flasquepart_box= new G4SubtractionSolid("Flasquepart_box", laitonsquare_box, Block_entrance_tub);

 G4SubtractionSolid* Flasquepart1_box= new G4SubtractionSolid("Flasquepart1_box", Flasquepart_box, Screwhole_tub, 0, G4ThreeVector(-1.2*cm,-1.2*cm,0));

 G4SubtractionSolid* Flasquepart2_box= new G4SubtractionSolid("Flasquepart2_box", Flasquepart1_box, Screwhole_tub, 0, G4ThreeVector(1.2*cm,-1.2*cm,0));

 G4SubtractionSolid* Flasquepart3_box= new G4SubtractionSolid("Flasquepart3_box", Flasquepart2_box, Screwhole_tub, 0, G4ThreeVector(-1.2*cm,1.2*cm,0));

 G4SubtractionSolid* Flasque_box= new G4SubtractionSolid("Flasque_box", Flasquepart3_box, Screwhole_tub, 0, G4ThreeVector(1.2*cm,1.2*cm,0));
 
 G4UnionSolid *Laiton_partial=new G4UnionSolid("Flasque+1Clinquant",Flasque_box,Clinquant_box,0,G4ThreeVector(flas_width/2+clin_width/2,0,clin_length/2)); 

G4UnionSolid *Laiton_frame=new G4UnionSolid("Flasque+2Clinquant",Laiton_partial,Clinquant_box,0,G4ThreeVector(-flas_width/2-clin_width/2,0,clin_length/2)); 

 G4LogicalVolume* PMCarrier_log= new G4LogicalVolume(Laiton_frame, mat, "PMCarrier_log",0,0,0);

  return PMCarrier_log;
}

G4VPhysicalVolume* DVCSCaloGeom::PMCarrier(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){
  
  G4VPhysicalVolume* PMCarrier_phys = new G4PVPlacement(0, Vect,
							Log1,Form("PMCarrier_%d",n),Log2,false,0, surf_check);
  
  return PMCarrier_phys;
}
	
// Methods to create CHcover /////////////////////////////////////////////////			 
				 
G4LogicalVolume* DVCSCaloGeom::CHcoverlog(G4Material* mat){

  G4double CHcover_length=1.*cm;
  G4double Frame_length=20.9*cm;
  G4double Side_width=2.*cm;
  G4double Side_height=54.*cm;

  G4double basis_width=39.9*cm;
  G4double lowpart_height=3.*cm;
  G4double highpart_height=2.*cm;

  G4Box* CHcover_box = new G4Box("CHcover_box",basis_width/2.,Side_height/2.,CHcover_length/2);
  
  G4LogicalVolume* CHcover_log= new G4LogicalVolume(CHcover_box, mat,"CHcover_log",0,0,0);

  return CHcover_log;
}

G4VPhysicalVolume* DVCSCaloGeom::CHcover(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect){

  G4double CHcover_length=1.*cm;
  G4double Frame_length=20.9*cm;
  
  G4VPhysicalVolume* CHcover_phys = new G4PVPlacement(0, Vect+G4ThreeVector(0,0,-(Frame_length+CHcover_length)/2.),
						      Log1,"CHcover",Log2,false,0, surf_check);
  
  return CHcover_phys;
}

// Methods to create CaloBlock in PbF2 ///////////////////////////////////////

G4LogicalVolume* DVCSCaloGeom::CaloBlocklog(G4Material* mat){

  //===========Remember that bl_width and bl_height should be 3.cm and 3.cm, instead of 3.014 cm
  G4double block_width =3.*cm;
  G4double block_height=3.*cm;
  G4double block_length=18.6*cm;

  G4double block_x = 0.999*block_width/2.;
  G4double block_y = 0.999*block_height/2.;
  G4double block_z = 0.999*block_length/2.;

  G4Box* CaloBlock_box = new G4Box("CaloBlock_box",block_x,block_y,block_z);
  
  G4LogicalVolume* CaloBlock_log= new G4LogicalVolume(CaloBlock_box, mat,"CaloBlock_log",0,0,0);

  return CaloBlock_log;
}

void DVCSCaloGeom::CaloBlock(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n){
  
  G4VPhysicalVolume* CaloBlock_phys = new G4PVPlacement(0, Vect,
							Log1,Form("CaloBlock_%d",n),Log2,false,100+n, surf_check);
  
  return;
}
				 
// Methods to create the frame of Aluminium which surrounds the blocks/////////

G4LogicalVolume* DVCSCaloGeom::Alframelog(G4Material* mat){
 
  G4double Frame_length=20.9*cm;
  G4double Side_width=2.*cm;
  G4double Side_height=54.*cm;

  G4double basis_width=39.9*cm;
  G4double lowpart_height=3.*cm;
  G4double highpart_height=2.*cm;

  G4Box* Side_box = new G4Box("Side_box",Side_width/2.,Side_height/2.,Frame_length/2); 

  G4Box* Lowpart_box = new G4Box("lowpart_box",basis_width/2.,lowpart_height/2.,Frame_length/2);

  G4Box* Highpart_box = new G4Box("highpart_box",basis_width/2.,highpart_height/2.,Frame_length/2);

  G4UnionSolid *Lowpart_Side1=new G4UnionSolid("Lowpart+1Side",Side_box,Lowpart_box,0,G4ThreeVector(-(basis_width+Side_width)/2.,(-Side_height+lowpart_height)/2.,0)); 

  G4UnionSolid *Lowpart_Side2=new G4UnionSolid("Lowpart+2Sides",Lowpart_Side1,Side_box,0,G4ThreeVector(-basis_width-Side_width,0,0));

  G4UnionSolid *Frame_box=new G4UnionSolid("Frame_box",Lowpart_Side2,Highpart_box,0,G4ThreeVector(-(basis_width+Side_width)/2.,(Side_height-highpart_height)/2.,0));

 G4LogicalVolume* Alframe_log= new G4LogicalVolume(Frame_box, mat, "Alframe_log",0,0,0);

  return Alframe_log;
}

G4VPhysicalVolume* DVCSCaloGeom::Alframe(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect){
  
  G4VPhysicalVolume* Alframe_phys = new G4PVPlacement(0, Vect,
						      Log1,"Alframe",Log2,false,0, surf_check);
  
  return Alframe_phys;
}
				 
// Methods to create Polyethylene shield ///////////////////////////////////

G4LogicalVolume* DVCSCaloGeom::CHshieldlog(G4Material* mat){

  G4double CHshield_height=60*cm;
  G4double CHshield_length=9.6*cm;
  G4double CHshieldtrap_width=CHshield_length*sin(15*degree);
  G4double CHshieldbox_width=45*cm-CHshieldtrap_width;

G4RotationMatrix *rot = new G4RotationMatrix();
rot->rotateX(90*degree);
 
G4Trap *CHshield_trap = new G4Trap("CHshield_trap", CHshield_height, CHshield_length, CHshieldtrap_width, 0.01*cm);

 G4Box *CHshield_box= new G4Box("CHshield_box", CHshieldbox_width/2., CHshield_height/2., CHshield_length/2.);

 G4UnionSolid* CHshield_final= new G4UnionSolid("CHshield_final",CHshield_box, CHshield_trap, rot ,G4ThreeVector(CHshieldbox_width/2+CHshieldtrap_width/3.5,0,0));
  
  
  G4LogicalVolume* CHshield_log= new G4LogicalVolume(CHshield_final, mat,"CHshield_log",0,0,0);

  return CHshield_log;
}

void DVCSCaloGeom::CHshield(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect){

  G4double CHshield_length=9.6*cm;
  G4double CHshieldtrap_width=CHshield_length*sin(15*degree);
  G4double block_length=18.6*cm;
  G4double CHcover_length=1.*cm;
  G4double flas_length=0.2*cm;
  G4double screw_length=1*cm;
  G4double BlackBox_length=0.6*cm;
  
  G4VPhysicalVolume* CHshield_phys = new G4PVPlacement(0, Vect+G4ThreeVector(-CHshieldtrap_width/2.,0,-block_length/2.-screw_length-flas_length-CHcover_length/2-BlackBox_length-CHshield_length/2.),
						       Log1,"CHshield_phys",Log2,false,0, surf_check);
  
  return;
}

// Methods to build the blackbox in which is the calorimeter, It is made of 3 boxes in the simulation/////////////////////////////////////////////////////////

// Outside Aluminium box of the Black box

G4LogicalVolume* DVCSCaloGeom::BB_OutAllog(G4Material* mat){

  G4double BB_length=91.6*cm;
  G4double BB_width=60.6*cm;
  G4double BB_height=89.21*cm;
  G4double Al_thickness=0.1*cm;

  G4Box* OutAl1_box = new G4Box("OutAl1_box",BB_width/2.,BB_height/2.,BB_length/2);
   G4Box* InnerAl1_box = new G4Box("InnerAl1_box",BB_width/2.-Al_thickness,BB_height/2.-Al_thickness,BB_length/2.-Al_thickness);
   G4SubtractionSolid* OutAl_box= new G4SubtractionSolid("OutAl_box", OutAl1_box, InnerAl1_box);

   G4LogicalVolume* OutAl_log = new G4LogicalVolume(OutAl_box,mat,"OutAl_log",0,0,0);

  return OutAl_log;
}

void DVCSCaloGeom::BB_OutAl(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect){

  G4double BB_length=91.6*cm;
  G4double BB_width=60.6*cm;
  G4double BB_height=89.21*cm;
  G4double Al_thickness=0.1*cm;
  G4double PVC_thickness=0.4*cm;
  G4double CHcover_length=1.*cm;
  G4double Side_width=2.*cm;
  G4double basis_width=39.9*cm;
  G4double Frame_length=20.9*cm;
  
  G4VPhysicalVolume* OutAl_phys = new G4PVPlacement(0,Vect+G4ThreeVector(
-BB_width/2.+2*Al_thickness+PVC_thickness+basis_width/2.+Side_width,
BB_height/2.-Al_thickness-PVC_thickness/2.-52.38*cm,
BB_length/2.-(Frame_length+CHcover_length)/2.-0.1*cm-CHcover_length/2.),Log1,"OutAl",Log2,false,0, surf_check);
  
  return;
}

// PVC box between the 2 Aluminium boxes 


G4LogicalVolume* DVCSCaloGeom::BB_PVClog(G4Material* mat){

  G4double BB_length=91.6*cm;
  G4double BB_width=60.6*cm;
  G4double BB_height=89.21*cm;
  G4double Al_thickness=0.1*cm;
  G4double PVC_thickness=0.4*cm;

   G4Box* InnerAl1_box = new G4Box("InnerAl1_box",BB_width/2.-Al_thickness,BB_height/2.-Al_thickness,BB_length/2.-Al_thickness);

   G4Box* InnerPVC_box = new G4Box("InnerPVC_box",BB_width/2.-Al_thickness-PVC_thickness,BB_height/2.-Al_thickness-PVC_thickness,BB_length/2.-Al_thickness-PVC_thickness);
   G4SubtractionSolid* PVC_box= new G4SubtractionSolid("PVC_box", InnerAl1_box, InnerPVC_box);

   G4LogicalVolume* PVC_log = new G4LogicalVolume(PVC_box,mat,"PVC_log",0,0,0);

  return PVC_log;
}

void DVCSCaloGeom::BB_PVC(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect){

  G4double BB_length=91.6*cm;
  G4double BB_width=60.6*cm;
  G4double BB_height=89.21*cm;
  G4double Al_thickness=0.1*cm;
  G4double PVC_thickness=0.4*cm;
  G4double CHcover_length=1.*cm;
  G4double Side_width=2.*cm;
  G4double basis_width=39.9*cm;
  G4double Frame_length=20.9*cm;
  
  G4VPhysicalVolume* PVC_phys = new G4PVPlacement(0,Vect+G4ThreeVector(
-BB_width/2.+2*Al_thickness+PVC_thickness+basis_width/2.+Side_width,
BB_height/2.-Al_thickness-PVC_thickness/2.-52.38*cm,
BB_length/2.-(Frame_length+CHcover_length)/2.-0.1*cm-CHcover_length/2.),Log1,"OutAl",Log2,false,0, surf_check);
  
  return;
}

// Inside Aluminium box of the Black box

G4LogicalVolume* DVCSCaloGeom::BB_InAllog(G4Material* mat){

  G4double BB_length=91.6*cm;
  G4double BB_width=60.6*cm;
  G4double BB_height=89.21*cm;
  G4double Al_thickness=0.1*cm;
  G4double PVC_thickness=0.4*cm;

  G4Box* InnerAl2_box = new G4Box("InnerAl1_box",BB_width/2.-2*Al_thickness-PVC_thickness,BB_height/2.-2*Al_thickness-PVC_thickness,BB_length/2.-2*Al_thickness-PVC_thickness);

   G4Box* InnerPVC_box = new G4Box("InnerPVC_box",BB_width/2.-Al_thickness-PVC_thickness,BB_height/2.-Al_thickness-PVC_thickness,BB_length/2.-Al_thickness-PVC_thickness);
   G4SubtractionSolid* InAl_box= new G4SubtractionSolid("InAl_box", InnerPVC_box, InnerAl2_box);

   G4LogicalVolume* InAl_log = new G4LogicalVolume(InAl_box,mat,"InAl_log",0,0,0);

  return InAl_log;
}

void DVCSCaloGeom::BB_InAl(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect){

  G4double BB_length=91.6*cm;
  G4double BB_width=60.6*cm;
  G4double BB_height=89.21*cm;
  G4double Al_thickness=0.1*cm;
  G4double PVC_thickness=0.4*cm;
  G4double CHcover_length=1.*cm;
  G4double Side_width=2.*cm;
  G4double basis_width=39.9*cm;
  G4double Frame_length=20.9*cm;
  
  G4VPhysicalVolume* InAl_phys = new G4PVPlacement(0,Vect+G4ThreeVector(
-BB_width/2.+2*Al_thickness+PVC_thickness+basis_width/2.+Side_width,
BB_height/2.-Al_thickness-PVC_thickness/2.-52.38*cm,
BB_length/2.-(Frame_length+CHcover_length)/2.-0.1*cm-CHcover_length/2.),Log1,"OutAl",Log2,false,0, surf_check);
  
  return;
}
