#include "dvcsRunAction.hh"

#include "G4Run.hh"
#include "dvcsHist_Manager.hh"

dvcsRunAction::dvcsRunAction(dvcsHistManager *hh):  hist_man(hh)
{;}


dvcsRunAction::~dvcsRunAction()
{}


void dvcsRunAction::BeginOfRunAction(const G4Run* aRun)
{
  hist_man->Book();
}


void dvcsRunAction::EndOfRunAction(const G4Run*)
{ 
  hist_man->Save();
}
