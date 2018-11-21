#ifndef DVCSCaloSD_h
#define DVCSCaloSD_h 1

#include "G4VSensitiveDetector.hh"
#include "DVCSCaloHit.hh"

#include "TH1.h"
//#include "TCaloEvent.h"
#include "TNtuple.h"
#include "TFile.h"
//#include "TCaloCluster.h"
#include "TLorentzVector.h"
//#include "TDVCSEvent.h"
//#include "TDVCSEventMC.h"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DVCSCaloSD : public G4VSensitiveDetector
{
  public:
      DVCSCaloSD(G4String);
     ~DVCSCaloSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

  //TDVCSEvent *ev;
  //TCaloEvent *caloev;
  TVector3 *electrowrite;
  TVector3 *vertexwrite;
  TFile *tf;
  TTree *t;

  private:
  DVCSCaloHitsCollection *trackerCollection;
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

