#ifndef DVCS_HRS_WINDOW_SD
#define DVCS_HRS_WINDOW_SD 1

#include <G4VSensitiveDetector.hh>
#include "dvcsHit.hh"
#include "G4ThreeVector.hh"

#include <TFile.h>
#include <TTree.h>

class G4Step;
class G4HCofThisEvent;

class HRSwindowSD : public G4VSensitiveDetector
{
public:
  HRSwindowSD(G4String);
  ~HRSwindowSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  

public:
  //double Get_em_energy();
  //  G4ThreeVector Get_em_momentum();

private:
  //TFile *file_out;
  //TTree *atree;
  
  G4double tot_em_energy;
  G4ThreeVector em_momentum;

  dvcsHitsCollection *trackerCollection;
};

#endif
