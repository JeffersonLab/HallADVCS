#ifndef DVCS_EVENT_ACTION
#define DVCS_EVENT_ACTION 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "/work/halla/dvcs/disk1/salina/soft-sfa/TDVCSEvent.h"
#include "/work/halla/dvcs/disk1/salina/soft-sfa/TCaloEvent.h"
#include "dvcsObj_Manager.hh"
#include "TLorentzVector.h"
#include "TRandom2.h"
#include <vector>
#include "/work/halla/dvcs/disk1/salina/soft-sfa/TARSWave2.h"
#include "TChain.h"
#include "TTree.h"
//#include "/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/RFunction.C"

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
 
  void DefinePrimaries(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, double);
  void DefineHRS_em(G4ThreeVector, double);
  void DefineWeights( double, double, double );
  void AddEdep(int, double);

private:
  
  //===================Primaries=============  
  
  bool hit_HRS;

  double Px_em_ini, Py_em_ini, Pz_em_ini, E_em_ini; // components of initial electron jst before interacting prevert radiated
  double Px_em_scat, Py_em_scat, Pz_em_scat, E_em_scat; // components of the scattered electron
  double Px_em_v, Py_em_v, Pz_em_v, E_em_v; // components of the scattered electron at vertex
  double Px_em_HRS, Py_em_HRS, Pz_em_HRS, E_em_HRS;  // components of the scettered electron when it reaches the HRS window
  double Px_phot1, Py_phot1, Pz_phot1, E_phot1; // components of the first photon from
  double Px_phot2, Py_phot2, Pz_phot2, E_phot2; // components of the first photon from

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
  // R_function *rfunc_gula;
};


#endif
