#include "dvcsHist_Manager.hh"
#include "G4UnitsTable.hh"

#include "TFile.h"
#include "TTree.h"
#include "dvcsGlobals.hh"

#include <iostream>

using namespace std;

dvcsHistManager::dvcsHistManager( ObjManager *aa )
{
  calo_ev = aa->GetCaloEvent_ptr();
  cout<<"Calo_ev pointer in HistManager constructor is "<<calo_ev<<endl;
}

void dvcsHistManager::Book()
{
  int run_number = dvcsGlobals::run_number;
  
  file_out = new TFile(Form("./Rootfiles/dvcs_sim_no_smear_kin%d_1.root",run_number), "Recreate");
  tr1 = new TTree("tr1", "dvcs_tree");
  
  tr1->Branch("calo_event.", "TCaloEvent", &calo_ev, 3200, 12);
  tr1->Branch("ene1", &ene1, "ene1/D");
  tr1->Branch("ene2", &ene2, "ene2/D");
  tr1->Branch("xc1", &xc1, "xc1/D");
  tr1->Branch("xc2", &xc2, "xc2/D");
  tr1->Branch("yc1", &yc1, "yc1/D");
  tr1->Branch("yc2", &yc2, "yc2/D");
  tr1->Branch("m", &m, "m/D");
  tr1->Branch("mm2", &mm2, "mm2/D");
  tr1->Branch("size1", &size1, "size1/I");
  tr1->Branch("size2", &size2, "size2/I");
  tr1->Branch("psf", &psf, "psf/D");
  tr1->Branch("kx", &kx, "kx/D");
  tr1->Branch("ky", &ky, "ky/D");
  tr1->Branch("kz", &kz, "kz/D");
  tr1->Branch("kx_v", &kx_v, "kx_v/D");
  tr1->Branch("ky_v", &ky_v, "ky_v/D");
  tr1->Branch("kz_v", &kz_v, "kz_v/D");
  tr1->Branch("kpx", &kpx, "kpx/D");
  tr1->Branch("kpy", &kpy, "kpy/D");
  tr1->Branch("kpz", &kpz, "kpz/D");
  tr1->Branch("kpx_v", &kpx_v, "kpx_v/D");
  tr1->Branch("kpy_v", &kpy_v, "kpy_v/D");
  tr1->Branch("kpz_v", &kpz_v, "kpz_v/D");
  tr1->Branch("q1x", &q1x, "q1x/D");
  tr1->Branch("q2x", &q2x, "q2x/D");
  tr1->Branch("q1y", &q1y, "q1y/D");
  tr1->Branch("q2y", &q2y, "q2y/D");
  tr1->Branch("q1z", &q1z, "q1z/D");
  tr1->Branch("q2z", &q2z, "q2z/D");
  tr1->Branch("q1x_v", &q1x_v, "q1x_v/D");
  tr1->Branch("q2x_v", &q2x_v, "q2x_v/D");
  tr1->Branch("q1y_v", &q1y_v, "q1y_v/D");
  tr1->Branch("q2y_v", &q2y_v, "q2y_v/D");
  tr1->Branch("q1z_v", &q1z_v, "q1z_v/D");
  tr1->Branch("q2z_v", &q2z_v, "q2z_v/D");
  tr1->Branch("v", &v, "v/D");
  tr1->Branch("v_v", &v_v, "v_v/D");
  tr1->Branch("t", &t, "t/D");
  tr1->Branch("tmin", &tmin,"tmin/D");
  tr1->Branch("t_v", &t_v, "t_v/D");
  tr1->Branch("tmin_v", &tmin_v, "tmin_v/D");
  tr1->Branch("xB", &xB, "xB/D");
  tr1->Branch("xB_v", &xB_v, "xB_v/D");
  tr1->Branch("q2", &q2, "q2/D");
  tr1->Branch("q2_v", &q2_v, "q2_v/D");
  tr1->Branch("phi", &phi, "phi/D");
  tr1->Branch("phi_v", &phi_v, "phi_v/D");
  tr1->Branch("crs_sum", &crs_sum, "crs_sum/D");
  tr1->Branch("crs_dif", &crs_dif, "crs_dif/D");
  tr1->Branch("rval", &rval, "rval/D");
  tr1->Branch("rval_v", &rval_v, "rval_v/D");
}

void dvcsHistManager::Fill_Ntuple()
{
  tr1->Fill();
}

void dvcsHistManager::Save()
{
  file_out=tr1->GetCurrentFile();
  tr1->Write();
  file_out->Close();
}

void dvcsHistManager::SetCaloData(double a_ene1, double a_ene2, double a_xc1, double a_xc2, 
				  double a_yc1, double a_yc2, int a_size1, int a_size2)
{
  ene1 = a_ene1;
  ene2 = a_ene2;
  xc1 = a_xc1;
  xc2 = a_xc2;
  yc1 = a_yc1;
  yc2 = a_yc2;
  size1 = a_size1;
  size2 = a_size2;
  //cout<<"xc = "<<(*xc)<<"   yc = "<<(*yc)<<endl;
}

void dvcsHistManager::SetBeamGen( double a_kx_v, double a_ky_v, double a_kz_v )
{
  kx_v = a_kx_v;
  ky_v = a_ky_v;
  kz_v = a_kz_v;
}

void dvcsHistManager::SetBeamRec( double a_kx, double a_ky, double a_kz )
{
  kx = a_kx;
  ky = a_ky;
  kz = a_kz;
}

void dvcsHistManager::SetScatemGen( double a_kpx_v, double a_kpy_v, double a_kpz_v )
{
  kpx_v = a_kpx_v;
  kpy_v = a_kpy_v;
  kpz_v = a_kpz_v;
}

void dvcsHistManager::SetScatemRec( double a_kpx, double a_kpy, double a_kpz )
{
  kpx = a_kpx;
  kpy = a_kpy;
  kpz = a_kpz;
}

void dvcsHistManager::SetPhotGen( double a_q1x_v, double a_q2x_v, double a_q1y_v, double a_q2y_v, 
				  double a_q1z_v, double a_q2z_v )
{
  q1x_v = a_q1x_v;
  q2x_v = a_q2x_v;
  q1y_v = a_q1y_v;
  q2y_v = a_q2y_v;
  q1z_v = a_q1z_v;
  q2z_v = a_q2z_v;
}

void dvcsHistManager::SetPhotRec( double a_q1x, double a_q2x, double a_q1y, double a_q2y, 
				  double a_q1z, double a_q2z )
{
  q1x = a_q1x;
  q2x = a_q2x;
  q1y = a_q1y;
  q2y = a_q2y;
  q1z = a_q1z;
  q2z = a_q2z;
}

void dvcsHistManager::SetkineGen( double a_t_v, double a_xB_v, double a_q2_v, double a_phi_v, double a_tmin_v)
{
  t_v = a_t_v;
  xB_v = a_xB_v;
  q2_v = a_q2_v;
  phi_v = a_phi_v;
  tmin_v = a_tmin_v;
}

void dvcsHistManager::SetkineRec( double a_t, double a_xB, double a_q2, double a_phi, double a_tmin)
{
  t = a_t;
  xB = a_xB;
  q2 = a_q2;
  phi = a_phi;
  tmin = a_tmin;
}

void dvcsHistManager::SetMM2(double a_mm2)
{
  mm2 = a_mm2;
}

void dvcsHistManager::SetMinv( double a_m )
{
  m = a_m;
}

void dvcsHistManager::SetVertexes(double a_vz, double a_vz_v )
{
  v = a_vz;
  v_v = a_vz_v;
}

void dvcsHistManager::SetWeights(double a_PSF, double a_crs_sum, double a_crs_dif)
{
  psf = a_PSF;
  crs_sum = a_crs_sum;
  crs_dif = a_crs_dif;
}

void dvcsHistManager::SetRvalGen(double a_rval_v)
{
  rval_v = a_rval_v;
}

void dvcsHistManager::SetRvalRec(double a_rval)
{
  rval = a_rval;
}
