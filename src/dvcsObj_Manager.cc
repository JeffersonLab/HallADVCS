#include "dvcsObj_Manager.hh"

ObjManager::ObjManager()
{
  calo_ev = new TCaloEvent();
  dvcs_ev = new TDVCSEvent();
  L_calo_phot = new TLorentzVector();
}

TCaloEvent *ObjManager::GetCaloEvent_ptr()
{
  return calo_ev;
}

TLorentzVector *ObjManager::GetCaloPhot_ptr()
{
  return L_calo_phot;
}

TDVCSEvent *ObjManager::GetDVCSEvent_ptr()
{
  return dvcs_ev;
}
