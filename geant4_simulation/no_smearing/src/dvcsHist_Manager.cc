#include "dvcsHist_Manager.hh"
#include "G4UnitsTable.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

dvcsHistManager::dvcsHistManager()
{
  
}

void dvcsHistManager::Book()
{
  file_out = new TFile("dvcs_sim_kine1.root", "Create");
  tr1 = new TTree("tr1", "dvcs_tree");
  
  tr1->Branch("ene", &ene, "ene/D");
  tr1->Branch("xc", &xc, "xc/D");
  tr1->Branch("yc", &yc, "yc/D");
  tr1->Branch("mm2", &mm2, "mm2/D");
  tr1->Branch("size", &size, "size/I");
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
  tr1->Branch("qx", &qx, "qx/D");
  tr1->Branch("qy", &qy, "qy/D");
  tr1->Branch("qz", &qz, "qz/D");
  tr1->Branch("qx_v", &qx_v, "qx_v/D");
  tr1->Branch("qy_v", &qy_v, "qy_v/D");
  tr1->Branch("qz_v", &qz_v, "qz_v/D");
  tr1->Branch("vz", &vz, "vz/D");
  tr1->Branch("vz_v", &vz_v, "vz_v/D");
  tr1->Branch("t", &t, "t/D");
  tr1->Branch("t_v", &t_v, "t_v/D");
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
  tr1->Write();
  file_out->Close();
}

void dvcsHistManager::SetCaloData(double a_ene, double a_xc, double a_yc, int a_size)
{
  ene = a_ene;
  xc = a_xc;
  yc = a_yc;
  size = a_size;
  
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

void dvcsHistManager::SetPhotGen( double a_qx_v, double a_qy_v, double a_qz_v )
{
  qx_v = a_qx_v;
  qy_v = a_qy_v;
  qz_v = a_qz_v;
}

void dvcsHistManager::SetPhotRec( double a_qx, double a_qy, double a_qz )
{
  qx = a_qx;
  qy = a_qy;
  qz = a_qz;
}

void dvcsHistManager::SetkineGen( double a_t_v, double a_xB_v, double a_q2_v, double a_phi_v)
{
  t_v = a_t_v;
  xB_v = a_xB_v;
  q2_v = a_q2_v;
  phi_v = a_phi_v;
}

void dvcsHistManager::SetkineRec( double a_t, double a_xB, double a_q2, double a_phi)
{
  t = a_t;
  xB = a_xB;
  q2 = a_q2;
  phi = a_phi;
}

void dvcsHistManager::SetMM2(double a_mm2)
{
  mm2 = a_mm2;
}

void dvcsHistManager::SetVertexes(double a_vz, double a_vz_v )
{
  vz = a_vz;
  vz_v = a_vz_v;
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
