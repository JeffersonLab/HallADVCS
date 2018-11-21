#ifndef DVCSCaloHit_h
#define DVCSCaloHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DVCSCaloHit : public G4VHit
{
  public:

      DVCSCaloHit();
     ~DVCSCaloHit();
      DVCSCaloHit(const DVCSCaloHit&);
      const DVCSCaloHit& operator=(const DVCSCaloHit&);
      G4int operator==(const DVCSCaloHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Print();

  public:
  
      void SetTrackID  (G4int track)      { trackID = track; };
      void SetBlockNb  (G4int block)      { blockNb = block; };  
      void SetEdep     (G4double de)      { edep = de; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      
      G4int GetTrackID()    { return trackID; };
      G4int GetBlockNb()    { return blockNb; };
      G4double GetEdep()    { return edep; };      
      G4ThreeVector GetPos(){ return pos; };
  
  private:
  
      G4int         trackID;
      G4int         blockNb;
      G4double      edep;
      G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<DVCSCaloHit> DVCSCaloHitsCollection;

extern G4Allocator<DVCSCaloHit> DVCSCaloHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* DVCSCaloHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) DVCSCaloHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DVCSCaloHit::operator delete(void *aHit)
{
  DVCSCaloHitAllocator.FreeSingle((DVCSCaloHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
