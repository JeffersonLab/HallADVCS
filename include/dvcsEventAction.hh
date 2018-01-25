#ifndef DVCS_EVENT_ACTION
#define DVCS_EVENT_ACTION 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "TDVCSEvent.h"
#include "TCaloEvent.h"
#include "dvcsObj_Manager.hh"
#include "TLorentzVector.h"
#include "TRandom2.h"

class G4Event;
class dvcsHistManager;
class TLorentzVector;

class dvcsEventAction: public G4UserEventAction
{
public:
  dvcsEventAction( dvcsHistManager*, ObjManager* );
  //dvcsEventAction( dvcsHist_Manager * );
  ~dvcsEventAction();
  
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
 
  void DefinePrimaries(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, double);
  void DefineHRS_em(G4ThreeVector, double);
  void DefineWeights( double, double, double );
  void AddEdep(int, double);

private:
  
  //===================Primaries=============  
  
  bool hit_HRS;

  double Px[5];
  double Py[5];
  double Pz[5];
  double E[5];
  double vert_z;
  double smear_vertz;
  double calo_edep_[208];
  
  double PSF;
  double crs_sum;
  double crs_dif;

  TRandom2 rand2;

  int ev_number;
  dvcsHistManager *hist_man;
  TDVCSEvent *dvcs_event;
  TCaloEvent *calo_event;
  TLorentzVector *L_calo_phot;
  // static const double Ebeam;
  // static const double Mp;
};


#endif
