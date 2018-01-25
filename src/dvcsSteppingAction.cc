#include "dvcsSteppingAction.hh"

#include "dvcsDetectorConstruction.hh"
#include "dvcsEventAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"



dvcsSteppingAction::dvcsSteppingAction( dvcsDetectorConstruction *dd, dvcsEventAction *evv ):
  det(dd), ev_act(evv)
{;}


dvcsSteppingAction::~dvcsSteppingAction()
{;}

void dvcsSteppingAction::UserSteppingAction( const G4Step *aStep )
{
  G4VPhysicalVolume* volume 
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4StepPoint* point1 = aStep->GetPreStepPoint();
  G4Track *atrack = aStep->GetTrack();
  
  int copy_nmb = volume->GetCopyNo();

  double energy = 0;
  G4ThreeVector em_HRS_mom(0, 0, 0);

  if( volume == det->GetHRS_window() && point1->GetStepStatus() == fGeomBoundary &&
      atrack->GetParticleDefinition()->GetPDGCharge() < 0 && atrack->GetTrackID() == 1)
    {
      energy = atrack->GetTotalEnergy()/GeV;
      em_HRS_mom = atrack->GetMomentum()/GeV;
      ev_act->DefineHRS_em(em_HRS_mom, energy);
    }
  else if( copy_nmb >=100 && copy_nmb < 308 ) // copy_numbers from 100 to 307 correspond to calo
    {
      ev_act->AddEdep(copy_nmb - 100, aStep->GetTotalEnergyDeposit()/GeV);
    }
    
  
}
