//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DVCSDetectorConstruction.hh,v 1.6 2006/06/29 17:47:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-03-patch-01 $
//

#ifndef DVCSCaloConstruction_H
#define DVCSCaloConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "DVCSCaloGeom.hh"

class DVCSCaloConstruction : public G4VUserDetectorConstruction
{
  public:

    DVCSCaloConstruction();
    ~DVCSCaloConstruction();

  void Construction(G4LogicalVolume* Log, G4double calo_dist, G4double calo_angle);

  private:
    
    // Logical volumes
    //

    G4LogicalVolume* calorimeterBlock_log;
    G4LogicalVolume* Calo_log;
    G4LogicalVolume* CaloBlock_log;
    G4LogicalVolume* Laiton_frame_log;
    G4LogicalVolume* G10_log;
    G4LogicalVolume* Screw_log;
    G4LogicalVolume* Tyvek_log;
    G4LogicalVolume* Tedlar_log;
    G4LogicalVolume* SShim_log;
    G4LogicalVolume* BShim_log;
    G4LogicalVolume* Frame_log;
    G4LogicalVolume* CHcover_log;
    G4LogicalVolume* OutAl_log;
    G4LogicalVolume* PVC_log;
    G4LogicalVolume* InAl_log;
  

    // Physical volumes
    //
 
    G4VPhysicalVolume* CaloBlock_phys;
    G4VPhysicalVolume* Calo_phys;
    G4VPhysicalVolume* Laiton_frame_phys;
    G4VPhysicalVolume* G10_phys;
    G4VPhysicalVolume* Screw_phys;
    G4VPhysicalVolume* Tyvek_phys;
    G4VPhysicalVolume* Tedlar_phys;
    G4VPhysicalVolume* BShim_phys;
    G4VPhysicalVolume* Frame_phys;
    G4VPhysicalVolume* CHcover_phys;
    G4VPhysicalVolume* OutAl_phys;
    G4VPhysicalVolume* PVC_phys;
    G4VPhysicalVolume* InAl_phys;
};

#endif

