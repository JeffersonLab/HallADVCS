#include "dvcsHRS_windowSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "G4RunManager.hh"
#include "dvcsRunAction.hh"
#include "dvcsPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

HRSwindowSD::HRSwindowSD(G4String name):
  G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");

  /*
  file_out = new TFile("Energy_In_HRS_window.root", "Recreate");
  
  atree = new TTree("em_energy", "Energy of em at HRS window");
  atree->Branch("tot_em_energy", &tot_em_energy, "tot_em_energy/D");
  */
}

HRSwindowSD::~HRSwindowSD()
{  
  //atree->Write();
  //file_out->Close();
}


void HRSwindowSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new dvcsHitsCollection

    (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
  
}

G4bool HRSwindowSD::ProcessHits(G4Step *aStep, G4TouchableHistory*)
{
  G4StepPoint* point1 = aStep->GetPreStepPoint();
  G4StepPoint* point2 = aStep->GetPostStepPoint();
  G4Track* track      = aStep->GetTrack();
  G4double tot_energy;
  
  if (point1->GetStepStatus() == fGeomBoundary && track->GetParticleDefinition()->GetPDGCharge() < 0 && track->GetTrackID() == 1)
    {
      tot_em_energy = track->GetTotalEnergy()/GeV;
      em_momentum = track->GetMomentum()/GeV;
    }
}

void HRSwindowSD::EndOfEvent(G4HCofThisEvent*)
{
  
}

