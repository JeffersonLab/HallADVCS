
#ifndef DVCSCaloGeom_H
#define DVCSCaloGeom_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

class DVCSCaloGeom : public G4VUserDetectorConstruction
{
  public:
  
    DVCSCaloGeom();
    ~DVCSCaloGeom();

  G4LogicalVolume* SShimlog(G4Material* mat);

  G4VPhysicalVolume* SShim(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);

   G4LogicalVolume* BShimlog(G4Material* mat);

  G4VPhysicalVolume* BShim(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);

  G4LogicalVolume* Screwlog(G4Material* mat);

  void Screw(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n); 

  G4LogicalVolume* G10log(G4Material* mat);

  G4VPhysicalVolume* G10(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);

  G4LogicalVolume* Tedlarlog(G4Material* mat);

  void Tedlar(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);


  G4LogicalVolume* Tyveklog(G4Material* mat);

  void Tyvek(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);

  G4LogicalVolume* CHcoverlog(G4Material* mat);

  G4VPhysicalVolume* CHcover(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect);

  G4LogicalVolume* Alframelog(G4Material* mat);

  G4VPhysicalVolume* Alframe(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect); 

  G4LogicalVolume* CaloBlocklog(G4Material* mat);

  void CaloBlock(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);

  G4LogicalVolume* PMCarrierlog(G4Material* mat);

  G4VPhysicalVolume* PMCarrier(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect, G4int n);

  G4LogicalVolume* CHshieldlog(G4Material* mat);

  void CHshield(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect);

  G4LogicalVolume* BB_OutAllog(G4Material* mat);

  void BB_OutAl(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect);

  G4LogicalVolume* BB_PVClog(G4Material* mat);

  void BB_PVC(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect);

  G4LogicalVolume* BB_InAllog(G4Material* mat);

  void BB_InAl(G4LogicalVolume* Log1, G4LogicalVolume* Log2, G4ThreeVector Vect);

  private:

  static const bool surf_check = 0;

    // Logical volumes
    //
    G4LogicalVolume* SShim_log;
    G4LogicalVolume* CHshield_log;
    G4LogicalVolume* CaloBlock_log;
    G4LogicalVolume* PMCarrier_log;
    G4LogicalVolume* G10_log;
    G4LogicalVolume* Screw_log;
    G4LogicalVolume* Tyvek_log;
    G4LogicalVolume* Tedlar_log;
    G4LogicalVolume* BShim_log;
    G4LogicalVolume* Alframe_log;
    G4LogicalVolume* CHcover_log;

    // Physical volumes
    //    
    G4VPhysicalVolume* SShim_phys;
    G4VPhysicalVolume* PMCarrier_phys;
    G4VPhysicalVolume* G10_phys; 
    G4VPhysicalVolume* CaloBlock_phys;
    G4VPhysicalVolume* Screw_phys;
    G4VPhysicalVolume* Tyvek_phys;
    G4VPhysicalVolume* Tedlar_phys;
    G4VPhysicalVolume* BShim_phys;
    G4VPhysicalVolume* Alframe_phys;
    G4VPhysicalVolume* CHcover_phys;
    G4VPhysicalVolume* CHshield_phys;
};

#endif

