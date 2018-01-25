#ifndef DVCS_HISTO_MANAGER
#define DVCS_HISTO_MANAGER 1

#include <TFile.h>
#include <TTree.h>

class dvcsHistManager
{
public:
  dvcsHistManager();
  
  void Book();
  void Fill_Ntuple();
  void Save();
  
  void SetCaloData(double, double, double, int); //ene, xc, yc, size
  void SetBeamGen(double, double, double); // kx_v, ky_v, kz_v
  void SetBeamRec(double, double, double); // kx, ky, kz
  void SetScatemGen(double, double, double); // kpx_v, kpy_v, kpz_v
  void SetScatemRec(double, double, double); // kpx, kpy, kpz
  void SetPhotGen(double, double, double); // kqx_v, kqy_v, kqz_v
  void SetPhotRec(double, double, double); // kqx, kqy, kqz
  void SetkineGen(double, double, double, double); // t, xB, Q2, phi(rad)
  void SetkineRec(double, double, double, double); // t, xB, Q2, phi(rad)
  void SetMM2( double ); // MM2
  void SetVertexes(double, double ); // vertex , smeared vertex
  void SetWeights(double, double, double); // PSF, crs_sum, crs_dif
  void SetRvalGen(double); // Generated R value
  void SetRvalRec(double); // Reconstructed R value
  
private:
  
  TFile *file_out;
  TTree *tr1;

  double ene; // energy of cluster in Calo
  double xc; // x coordinate of cluster in Calo
  double yc; // y coordinate of cluster in Calo
  double mm2; //missing mass squared
  int size; // number of blocks in a cluster
  double psf; // Phase Spase Factor
  double kx; // x component of beam electron's momentum
  double ky; // y component of beam electron's momentum
  double kz; // z component of beam electron's momentum
  double kx_v; // x component of reconstructed beam electron's momentum
  double ky_v; // y component of reconstructed beam electron's momentum
  double kz_v; // z component of reconstructed beam electron's momentum
  double kpx; // x component of scattered electron's momentum at HRS window
  double kpy; // y component of scattered electron's momentum at HRS window
  double kpz; // z component of scattered electron's momentum at HRS window
  double kpx_v; // x component of scattered electron's momentum at interaction vertex
  double kpy_v; // y component of scattered electron's momentum at interaction vertex
  double kpz_v; // z component of scattered electron's momentum at interaction vertex
  double qx; // x component of photon's momentum from calo
  double qy; // y component of photon's momentum from calo
  double qz; // z component of photon's momentum from calo
  double qx_v; // x component of photon's momentum at interaction vertex (from generatot)
  double qy_v; // y component of photon's momentum at interaction vertex (from generatot)
  double qz_v; // z component of photon's momentum at interaction vertex (from generatot)
  double vz; // z component of interaction coordinate This is smeareed from the vz_v(generated value) by sigma = 1.2(mm)/sin(theta_HRS)
  double vz_v; // z component of interaction coordinate taken from generetor
  double t; // t calculated using angles of reconstructed particles
  double t_v; // t calculated using information of generated particles
  double xB; // xB calculated using information of reconstructed particles
  double xB_v; // xB calculated using information of generated particles
  double q2; // Q2 calculated from reconstructed electron't kinematics
  double q2_v; // Q2 calculated at interaction vertex
  double phi; // phi calculated using reconstructed particles
  double phi_v; // phi calculated at interaction vertex
  double crs_sum; // Sigma plus + Sigma minus
  double crs_dif; // Sigma plus - Sigma minus
  double rval; // r_value of HRS acceptance
  double rval_v; // r_value of HRS acceptance at production vertex
};

#endif
