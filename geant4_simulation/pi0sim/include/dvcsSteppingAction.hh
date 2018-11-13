#ifndef DVCS_STEPPING_ACTION
#define DVCS_STEPPING_ACTION 1

#include "G4UserSteppingAction.hh"
#include "dvcsDetectorConstruction.hh"
#include "dvcsEventAction.hh"

//class dvcsDetectorConstruction;
//class dvcsEventAction;

class dvcsSteppingAction : public G4UserSteppingAction
{
public:
  dvcsSteppingAction( dvcsDetectorConstruction *, dvcsEventAction* );
  ~dvcsSteppingAction();

  void UserSteppingAction( const G4Step* );

private:
  dvcsDetectorConstruction *det;
  dvcsEventAction *ev_act;
};

#endif
