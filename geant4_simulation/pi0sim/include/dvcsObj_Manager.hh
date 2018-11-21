#ifndef DVCS_OBJ_MANAGER
#define DVCS_OBJ_MANAGER 1

#include "/work/halla/dvcs/disk1/salina/soft-2/TCaloEvent.h"
#include "/work/halla/dvcs/disk1/salina/soft-2/TDVCSEvent.h"
#include "TLorentzVector.h"

class ObjManager
{
public:
  ObjManager();
  
  TCaloEvent *GetCaloEvent_ptr();
  TLorentzVector *GetCaloPhot_ptr();
  TDVCSEvent *GetDVCSEvent_ptr();

private:
  TCaloEvent *calo_ev;
  TDVCSEvent *dvcs_ev;
  TLorentzVector *L_calo_phot;
};

#endif
