#include "dvcsHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<dvcsHit> dvcsHitAllocator;

dvcsHit::dvcsHit() {}

dvcsHit::~dvcsHit() {}

dvcsHit::dvcsHit( const dvcsHit &right ) : G4VHit()
{
  trackID = right.trackID;
  blockNb = right.blockNb;
  edep    = right.edep;
  pos     = right.pos;
}

const dvcsHit & dvcsHit::operator=(const dvcsHit& right)
{
  trackID = right.trackID;
  blockNb = right.blockNb;
  edep    = right.edep;
  pos     = right.pos;
  return *this;
}

G4int dvcsHit::operator==( const dvcsHit& right )const
{
  return (this==&right) ? 1 : 0;
}

void dvcsHit::Print()
{
  //  G4cout << "  trackID: " << trackID << "  blockNb: " << blockNb
  //     << "  energy deposit: " << G4BestUnit(edep,"Energy")
  //     << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}
