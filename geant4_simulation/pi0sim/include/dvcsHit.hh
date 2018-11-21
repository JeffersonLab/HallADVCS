#ifndef DVCS_HIT
#define DVCS_HIT 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class dvcsHit : public G4VHit
{
public:
  dvcsHit();
  ~dvcsHit();
  dvcsHit(const dvcsHit&);
  const dvcsHit& operator=(const dvcsHit&);
  G4int operator==(const dvcsHit&) const;
  
  inline void* operator new(size_t);
  inline void operator delete(void*);
  
  void Print();

  void SetTrackID (G4int track) {trackID = track;}
  void SetBlockNb (G4int block) {blockNb = block;}
  void SetEdep    (G4double de) {edep = de;}
  void SetPos(G4ThreeVector xyz) { pos = xyz;}
  
  G4int GetTrackID () {return trackID;}
  G4int GetBlockNb () {return blockNb;}
  G4double GetEdep () {return edep;}
  G4ThreeVector GetPos() { return pos; }

private:
  
  G4int trackID;
  G4int blockNb;
  G4double edep;
  G4ThreeVector pos;
};

typedef G4THitsCollection<dvcsHit> dvcsHitsCollection;

extern G4Allocator<dvcsHit> dvcsHitAllocator;

inline void* dvcsHit::operator new(size_t)
{
  void *aHit;
  aHit = (void*)dvcsHitAllocator.MallocSingle();
  return aHit;
}

inline void dvcsHit::operator delete(void *aHit)
{
  dvcsHitAllocator.FreeSingle((dvcsHit*) aHit);
}

#endif
