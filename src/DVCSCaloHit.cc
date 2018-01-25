#include "DVCSCaloHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<DVCSCaloHit> DVCSCaloHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DVCSCaloHit::DVCSCaloHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DVCSCaloHit::~DVCSCaloHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DVCSCaloHit::DVCSCaloHit(const DVCSCaloHit& right)
  : G4VHit()
{
  trackID   = right.trackID;
  blockNb   = right.blockNb;
  edep      = right.edep;
  pos       = right.pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DVCSCaloHit& DVCSCaloHit::operator=(const DVCSCaloHit& right)
{
  trackID   = right.trackID;
  blockNb = right.blockNb;
  edep      = right.edep;
  pos       = right.pos;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DVCSCaloHit::operator==(const DVCSCaloHit& right) const
{
  return (this==&right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DVCSCaloHit::Print()
{
  // G4cout << "  trackID: " << trackID << "  blockNb: " << blockNb
  //        << "  energy deposit: " << G4BestUnit(edep,"Energy")
  // 	    << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
