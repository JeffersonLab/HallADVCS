#include "dvcsEventAction.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "dvcsHist_Manager.hh"
#include "G4EventManager.hh"
#include "dvcsPrimaryGeneratorAction.hh"
#include "TVector3.h"
#include "TProcessID.h"
#include "TDVCSEvent.h"
#include "TCaloEvent.h"
#include "TDVCSGlobal2h"
#include "dvcsGlobals.hh"
#include "TRandom2.h"
//#include "/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/RFunction.C"
#include "/work/halla/dvcs/disk1/salina/RFunction_110717.C"
// const double dvcsEventAction::Ebeam = 5.552;
// const double dvcsEventAction::Mp = 0.938;

dvcsEventAction::dvcsEventAction( dvcsHistManager *aa, ObjManager *bb )
{
  hist_man = aa;
  int run_number = dvcsGlobals::run_number;
  gdvcs->SetWF(0);gdvcs->SetRun(run_number);gdvcs->ForceUpdate(); // This should organized by an automatic way !! The production is 8000
  
  rand2.SetSeed(0);
  
  cout<<"In Event Action Run number is "<<run_number<<endl;

  dvcs_event = bb->GetDVCSEvent_ptr();
  calo_event = bb->GetCaloEvent_ptr();
  //L_calo_phot = bb->GetCaloPhot_ptr();
  
  dvcs_event->SetCaloEvent(calo_event);
  
  cout<<"Caloriemeter angle (deg): "<<dvcs_event->GetGeometry()->GetCaloTheta()*TMath::RadToDeg()<<endl;
  cout<<"Calorimter distance (cm): "<<dvcs_event->GetGeometry()->GetCaloDist()<<endl;
  cout<<"Calorimter X offset (cm): "<<calo_event->GetGeometry()->GetCenterXPos()<<endl;
  cout<<"Calorimter Y offset (cm): "<<calo_event->GetGeometry()->GetCenterYPos()<<endl;

  
  for( int i = 0; i < 208; i++ )
    {
      cout<<"Calo block position x and y (cm): "<<calo_event->GetGeometry()->GetBlockXPos(i)<<"\t"
  	  <<calo_event->GetGeometry()->GetBlockYPos(i)<<endl;
    }
  
  hit_HRS = 0;
  cout<<"Event Action Constructor"<<endl;
}

dvcsEventAction::~dvcsEventAction()
{
}

void dvcsEventAction::BeginOfEventAction(const G4Event* evt)
{
  if( dvcsGlobals::hit_HRS_CALO_flag )
    {
      Px_em_HRS = 0;
      Py_em_HRS = 0;
      Pz_em_HRS = 0;
      E_em_HRS = 0;
            
      for( int i = 0; i < 208; i++ )
	{
	  calo_edep_[i] = 0;
	}
    }
}

void dvcsEventAction::EndOfEventAction(const G4Event* evt)
{
  if( dvcsGlobals::hit_HRS_CALO_flag && hit_HRS)
    {
      double Ebeam = dvcsGlobals::Ebeam;
      double Mp = 0.938;
      double mpi = 0.1349766;
      
      double clust_x1;
      double clust_y1;
      float clust_E1;
      int clust_size1;
      
      float phot_px1;
      float phot_py1;
      float phot_pz1;

      double clust_x2;
      double clust_y2;
      float clust_E2;
      int clust_size2;
      
      float phot_px2;
      float phot_py2;
      float phot_pz2;


      double fphi;
      double tM;
      double mm2;
      double minv;

      Double_t tjitter;
      
      int n_clust;

      double xB;
      double Q2;

      int run_number = dvcsGlobals::run_number;

      double trigsimthresh;
      if(run_number==361)trigsimthresh=1.1;
      if(run_number==362 ||run_number==363)trigsimthresh=1.6;
      if(run_number==481)trigsimthresh=0.5;
      if(run_number==482)trigsimthresh=0.9;
      if(run_number==483)trigsimthresh=1.1;
      if(run_number==484)trigsimthresh=1.5;
      if(run_number==601)trigsimthresh=0.8;
      if(run_number==603)trigsimthresh=1.0;
     
      Int_t ObjectNumber=TProcessID::GetObjectCount();

      for( int j = 0; j < 208; j++ )
	{
	  TCaloBlock *block=calo_event->AddBlock(j);// j is block number
	  double e_smear = rand2.Gaus(calo_edep_[j], sqrt(calo_edep_[j]/175.));
	  //block->SetBlockEnergy(e_smear);//smeared calo energy
	  block->SetBlockEnergy(calo_edep_[j]);//No smearing
	}
      
      dvcs_event->SetVertex(0, 0, smear_vertz);
      calo_event->TriggerSim(trigsimthresh);
      calo_event->DoClustering();
      n_clust = calo_event->GetNbClusters();

      if( n_clust >= 2 )
	{
	  //  1-st photon
	  calo_event->GetCluster(0)->Analyze();
	  TLorentzVector L_calo_phot1 = dvcs_event->GetPhoton(0, 7, 0);
	  L_calo_phot1 = L_calo_phot1;//1.04*L_calo_phot1;

	  clust_x1 = calo_event->GetCluster(0)->GetX();
	  clust_y1 = calo_event->GetCluster(0)->GetY();
	  clust_E1 = L_calo_phot1.E();//calo_event->GetCluster(0)->GetE();
	  clust_size1 = calo_event->GetCluster(0)->GetClusSize();
	  
	  phot_px1 = L_calo_phot1.Px();
	  phot_py1 = L_calo_phot1.Py();
	  phot_pz1 = L_calo_phot1.Pz();
	  // 2-nd photon
	  calo_event->GetCluster(1)->Analyze();
	  TLorentzVector L_calo_phot2 = dvcs_event->GetPhoton(1, 7, 0);
	  L_calo_phot2 = L_calo_phot2;//1.04*L_calo_phot2;

	  clust_x2 = calo_event->GetCluster(1)->GetX();
	  clust_y2 = calo_event->GetCluster(1)->GetY();
	  clust_E2 = L_calo_phot2.E();//calo_event->GetCluster(1)->GetE();
	  clust_size2 = calo_event->GetCluster(1)->GetClusSize();
	  
	  phot_px2 = L_calo_phot2.Px();
	  phot_py2 = L_calo_phot2.Py();
	  phot_pz2 = L_calo_phot2.Pz();
       
	  TLorentzVector k(Px_em_ini,Py_em_ini,Pz_em_ini,E_em_ini);
	  
	  //cout<<"Beam energy at vertex = "<<E_em_ini<<"  Incident beam energy "<<Ebeam<<endl;
	  
	  TLorentzVector kp(Px_em_HRS, Py_em_HRS, Pz_em_HRS,E_em_HRS);
	  TVector3 q3pion(phot_px1 + phot_px2, phot_py1 + phot_py2, phot_pz1 + phot_pz2); // 3 vector of 2 photons
	  TLorentzVector qpi0;//
	  TLorentzVector pion2=L_calo_phot1 + L_calo_phot2;
	  qpi0.SetPxPyPzE(phot_px1+phot_px2,phot_py1+phot_py2,phot_pz1+phot_pz2,clust_E1+clust_E2);
	  //cout<<"pi0 inv. mass = "<<qpi0.M()<<" pion2.M()  = "<<pion2.M()<<endl;
	  Q2 =-(kp-k).Mag2();//2*k.Mag()*kp.Mag()*(1 - cos(k.Angle(kp)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))
	  double nu=E_em_ini-(kp.Vect()).Mag();
	  xB = Q2/(2.*Mp*nu);
	  TVector3 q = k.Vect() - kp.Vect();//virtual photon momenta
	  double cosqqp=(q.Dot(qpi0.Vect()))/(q.Mag()*((qpi0).Vect()).Mag());//cosine of angle between vitrual photon and pion  
	  double sintheta_gg=sin(acos(cosqqp));
	  double q_exp=sqrt(nu*nu+Q2);
	  double s = Mp*Mp+Q2*((1.-xB)/xB);
	  double W = sqrt(s);
	  double F2 = (s-mpi*mpi-Mp*Mp)*(s-mpi*mpi-Mp*Mp)/4.-(mpi*mpi*Mp*Mp);
	  double qpmax =((s+mpi*mpi-pow(Mp,2))*q_exp +2*(nu+Mp)*sqrt(F2))/(2*pow(W,2));
	  double q0pmax = sqrt(pow(mpi,2)+pow(qpmax,2));
	  double qp=((pow(W,2)+pow(mpi,2)-pow(Mp,2))*q_exp*cosqqp+2*(nu+Mp)*sqrt(F2-pow(mpi*q_exp*sintheta_gg,2)))/(2*(pow(W,2)+pow(q_exp*sintheta_gg,2)));
	  double q0p = sqrt(pow(mpi,2)+pow(qp,2));
	  double qpmin =((pow(W,2)+pow(mpi,2)-pow(Mp,2))*q_exp-2*(nu+Mp)*sqrt(F2))/(2*pow(W,2));
	  double q0pmin = sqrt(pow(mpi,2)+pow(qpmin,2));
	  double t_min = -Q2+pow(mpi,2)-2*nu*q0pmax+2*q_exp*qpmax;
	  tM = -Q2+pow(mpi,2)-2*nu*q0p+2*q_exp*qp*cosqqp;
	  
	  //cout<<"tmin Eric = "<<t_min<<" t Eric = "<<tM<<endl; 
	  //TVector3 q = k - kp;
      
	  TVector3 v1=q.Cross(kp.Vect());
	  TVector3 v2=q.Cross(q3pion);
	  fphi = v1.Angle(v2);
	  if(q.Dot(v1.Cross(v2))<0) fphi = 2.*TMath::Pi()-fphi;
	  // double cos=(q.Dot(qp))/(q.Mag()*qp.Mag());
	  // double tM=(Q2*Mp+2.*nu*Mp*(nu-sqrt(nu*nu+Q2)*cos))/(sqrt(nu*nu+Q2)*cos-nu-Mp);
	  double t_min_2=(Q2*Mp+2.*nu*Mp*(nu-sqrt(nu*nu+Q2)))/(sqrt(nu*nu+Q2)-nu-Mp);    
	  TLorentzVector Lb, Lp, Lkp, L_mis;
	  Lb.SetPxPyPzE(Px_em_ini, Py_em_ini, Pz_em_ini, E_em_ini);
	  Lp.SetPxPyPzE(0, 0, 0, Mp);
	  Lkp.SetPxPyPzE(Px_em_HRS, Py_em_HRS, Pz_em_HRS, E_em_HRS);
	  
	  //tM = (Lb - Lkp - L_calo_phot2 - L_calo_phot1).M2();
	  //cout<<"tmin_me_Rafo = "<<t_min_2<<" t 4 vectors = "<<tM<<endl;

	  //      cout<<"Lgp.M() = "<<(Lgp.M())<<"   Lgp.P() = "<<(Lgp.P())<<endl;
	  //cout<<"Lkp.M() = "<<(Lkp.M())<<"   Lkp.P() = "<<(Lkp.P())<<endl;
	  
	  mm2 = (Lb + Lp - L_calo_phot2 - L_calo_phot1 - Lkp).M2();
	  minv = (L_calo_phot1 + L_calo_phot2).M();
	  
	  // cout<<"n_rec = "<<n_rec<<endl;
	  // cout<<" clust_E1[0] =  "<<clust_E1[0]<<" clust_E2[0] =  "<<clust_E2[0]<<endl;
      
	  hist_man->SetCaloData(clust_E1, clust_E2, clust_x1, clust_x2, clust_y1, clust_y2, clust_size1, clust_size2);
	  hist_man->SetPhotRec(phot_px1, phot_px2, phot_py1, phot_py2, phot_pz1, phot_pz2); //Set Reconstructed photon kinematics: Px, Py, Pz
	  hist_man->SetkineRec(tM, xB, Q2, fphi,t_min);
	  hist_man->SetMM2(mm2);
	  hist_man->SetMinv(minv);
	  hist_man->Fill_Ntuple();
	}
      calo_event->Reset();

      TProcessID::SetObjectCount(ObjectNumber);
    }
  
  hit_HRS = false;

}

void dvcsEventAction::DefinePrimaries(TLorentzVector *L_em_init, TLorentzVector *L_em_scat, TLorentzVector *L_em_scat_v, 
				      TLorentzVector *L_phot1, TLorentzVector *L_phot2, TLorentzVector *L_prot, double z)
{
  Px_em_ini = L_em_init->Px(); Px_em_scat = L_em_scat->Px();  Px_phot1 = L_phot1->Px(); Px_phot2 = L_phot2->Px();
  Py_em_ini = L_em_init->Py(); Py_em_scat = L_em_scat->Py();  Py_phot1 = L_phot1->Py(); Py_phot2 = L_phot2->Py();
  Pz_em_ini = L_em_init->Pz(); Pz_em_scat = L_em_scat->Pz();  Pz_phot1 = L_phot1->Pz(); Pz_phot2 = L_phot2->Pz();
  E_em_ini = L_em_init->E();  E_em_scat = L_em_scat->E();     E_phot1  = L_phot1->E();  E_phot2  = L_phot2->E();

  double HRS_angle = dvcsGlobals::HRS_angle;
  double Ebeam = dvcsGlobals::Ebeam;
  vert_z = z;
  smear_vertz = vert_z + rand2.Gaus(0, 1.2/sin(HRS_angle))*0.1; // 0.1 is for mm->cm
  
  hist_man->SetVertexes(smear_vertz, vert_z );
  
  //cout<<"vertz = "<<vert_z<<"\t vert_z_v = "<<smear_vertz<<endl;

  hist_man->SetScatemGen(Px_em_scat, Py_em_scat, Pz_em_scat);
  hist_man->SetPhotGen(Px_phot1, Px_phot2, Py_phot1, Py_phot2, Pz_phot1, Pz_phot2);
  double beam_x = 0;
  double beam_y = 0;
  double beam_z = E_em_ini;
  hist_man->SetBeamGen(Px_em_ini, Py_em_ini, Pz_em_ini);
  hist_man->SetBeamRec(beam_x, beam_y, beam_z);

  //===========Calculate Q2, t, phi=========
  
  TLorentzVector L_pi0 = *L_phot1 + *L_phot2;

  const double Mp = 0.938;
  const double mpi = 0.1349766;
  
  TVector3 k_v(0, 0, E_em_ini);
  
  TVector3 kp_v = L_em_scat_v->Vect();//L_em_scat->Vect();----changed.....Mongi
  TVector3 qp_v = L_pi0.Vect();
  
  double Q2_v = 2*k_v.Mag()*kp_v.Mag()*(1 - cos(k_v.Angle(kp_v)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))
  
  TVector3 q_v = k_v - kp_v;
  
  TVector3 v1_v=q_v.Cross(kp_v);
  TVector3 v2_v=q_v.Cross(qp_v);
  double fphi_v=v1_v.Angle(v2_v);
  if(q_v.Dot(v1_v.Cross(v2_v))<0) fphi_v=2.*TMath::Pi()-fphi_v;
  
  double cos_v=(q_v.Dot(qp_v))/(q_v.Mag()*qp_v.Mag());
  double nu_v=E_em_ini-kp_v.Mag();
  double xB_v = Q2_v/(2*Mp*nu_v);
  double cosqqp=cos_v;//(q_v.Dot(L_pi0.Vect()))/(q_v.Mag()*((L_pi0).Vect()).Mag());//cosine of angle between vitrual photon and pion                        
  double sintheta_gg=sin(acos(cosqqp));
  double q_exp=sqrt(nu_v*nu_v+Q2_v);
  double s = Mp*Mp+Q2_v*((1.-xB_v)/xB_v);
  double W = sqrt(s);
  double F2 = (s-mpi*mpi-Mp*Mp)*(s-mpi*mpi-Mp*Mp)/4.-(mpi*mpi*Mp*Mp);
  double qpmax =((s+mpi*mpi-pow(Mp,2))*q_exp +2*(nu_v+Mp)*sqrt(F2))/(2*pow(W,2));
  double q0pmax = sqrt(pow(mpi,2)+pow(qpmax,2));
  double qp=((pow(W,2)+pow(mpi,2)-pow(Mp,2))*q_exp*cosqqp+2*(nu_v+Mp)*sqrt(F2-pow(mpi*q_exp*sintheta_gg,2)))/(2*(pow(W,2)+pow(q_exp*sintheta_gg,2)));
  double q0p = sqrt(pow(mpi,2)+pow(qp,2));
  double qpmin =((pow(W,2)+pow(mpi,2)-pow(Mp,2))*q_exp-2*(nu_v+Mp)*sqrt(F2))/(2*pow(W,2));
  double q0pmin = sqrt(pow(mpi,2)+pow(qpmin,2));
  //double t_min = -Q2_v+pow(mpi,2)-2*nu_v*q0pmax+2*q_exp*qpmax;
  //tM = -Q2_v+pow(mpi,2)-2*nu_v*q0p+2*q_exp*qp*cosqqp;
  double tMin_v=-Q2_v+pow(mpi,2)-2*nu_v*q0pmax+2*q_exp*qpmax;//(Q2_v*Mp+2.*nu_v*Mp*(nu_v-sqrt(nu_v*nu_v+Q2_v)))/(sqrt(nu_v*nu_v+Q2_v)-nu_v-Mp);
  double tM_v = -Q2_v+pow(mpi,2)-2*nu_v*q0p+2*q_exp*qp*cosqqp;//2*Mp*(Mp - L_prot->E());

  hist_man->SetkineGen(tM_v, xB_v, Q2_v, fphi_v,tMin_v);

  // TLorentzVector L_Eb_tmp;  L_Eb_tmp.SetPxPyPzE(0, 0, Ebeam, Ebeam);
  // TLorentzVector L_init_phot = L_Eb_tmp - *L_em_scat;

  // cout<<" t (using only angles) = "<<tM_v<<endl;
  // cout<<" t (using momenta) = "<<(*L_phot - L_init_phot).M2()<<endl;
  // cout<<" Their Differrence is  "<<(tM_v - (*L_phot - L_init_phot).M2())<<endl;
  

  //============Calculate R_value===========
  
  //  Double_t HRS_angle= dvcsGlobals::HRS_angle;  //23.91*TMath::Pi()/180.; // KINEMATICS 3
  Double_t pcentral= dvcsGlobals::HRS_momentum; // KINEMATICS 3
  int run_number = dvcsGlobals::run_number;
//commenting out nov-8-2018
 /* int run=10555;//set default run to kin36_1 for Alexa's rfunction
  if(run_number==362){run=14152;}
  else if(run_number==363){run=14480;}
  //else{cout<<"ERROR : run number out of range, rfunction will be set to the default kinematic, i.e kin36_1 :"<<endl;}
  
  //HRS_angle
  Double_t ry = -0.01*vert_z*TMath::Sin(TMath::ATan2(Px_em_scat,Pz_em_scat)); // 0.01 to convert cm to meters
  Double_t rdp= (sqrt(Px_em_scat*Px_em_scat+Py_em_scat*Py_em_scat+Pz_em_scat*Pz_em_scat)-pcentral)/pcentral;
  Double_t rtheta = TMath::ATan2(-Py_em_scat,TMath::Sqrt(Px_em_scat*Px_em_scat+Pz_em_scat*Pz_em_scat));
  Double_t rphi   = TMath::ATan2(Px_em_scat,Pz_em_scat) - HRS_angle;
  //R_function *rfunc_gula= new R_function(Form("/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/dis_cut_%d.root",run_number));
  double r_val = RFunction(run,rtheta,rdp,rphi,ry);//Alexa's Rfunction//dvcs_event->rfunction(ry,rdp,rtheta,rphi);//old Rfunction
    //rfunc_gula->Global_R_function(rphi,rdp,rtheta,ry);//dvcs_event->rfunction(ry,rdp,rtheta,rphi);
  hist_man->SetRvalGen(r_val);*/
//end of commented out r val nov-8-2018
}

void dvcsEventAction::DefineWeights(double a_PSF, double a_crs_sum, double a_crs_dif)
{
  PSF = a_PSF;
  crs_sum = a_crs_sum;
  crs_dif = a_crs_dif;

  hist_man->SetWeights(PSF, crs_sum, crs_dif);

  //cout<<"In Event Action psf = "<<PSF<<"  crs_sum = "<<crs_sum<<"   crs_dif = "<<crs_dif<<endl;
}

void dvcsEventAction::DefineHRS_em(G4ThreeVector mom, double e)
{
  Px_em_HRS = mom.getX();
  Py_em_HRS = mom.getY();
  Pz_em_HRS = mom.getZ();
  E_em_HRS = e;

  hist_man->SetScatemRec(Px_em_HRS, Py_em_HRS, Pz_em_HRS);
  hit_HRS = true;

  //Calculate R_value

  double HRS_angle = dvcsGlobals::HRS_angle;
  Double_t pcentral= dvcsGlobals::HRS_momentum; // KINEMATICS 3

  int run_number = dvcsGlobals::run_number;
  int run=10555;//set default run to kin36_1 for Alexa's rfunction                                                                                          
  if(run_number==362){run=14152;}
  else if(run_number==363){run=14480;}
  else if(run_number==481){run=12518;}
  else if(run_number==482){run=13009;}//check with alexa, changed from 13013
  else if(run_number==483){run=12843;}
  else if(run_number==484){run=13136;}
  else if(run_number==601){run=15017;}
  else if(run_number==602){run=14628;}

//  cout<<"The run number is "<<run_number<<" and run is "<<run<<endl;


  //else{cout<<"ERROR : run number out of range, rfunction will be set to the default kinematic, i.e kin36_1 :"<<endl;}

  //HRS_angle
  Double_t ry = -0.01*smear_vertz*TMath::Sin(TMath::ATan2(Px_em_HRS,Pz_em_HRS));  // 0.01 to convert cm to meters
  Double_t rdp= (sqrt(Px_em_HRS*Px_em_HRS+Py_em_HRS*Py_em_HRS+Pz_em_HRS*Pz_em_HRS)-pcentral)/pcentral;
  Double_t rtheta = TMath::ATan2(-Py_em_HRS,TMath::Sqrt(Px_em_HRS*Px_em_HRS+Pz_em_HRS*Pz_em_HRS));
  Double_t rphi   = TMath::ATan2(Px_em_HRS,Pz_em_HRS) - HRS_angle;
  //R_function *rfunc_gula= new R_function(Form("/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/dis_cut_%d.root",run_number));
  //double r_val = rfunc_gula->Global_R_function(rphi,rdp,rtheta,ry);//dvcs_event->rfunction(ry,rdp,rtheta,rphi); 
  double r_val = RFunction(run,rtheta,rdp,rphi,ry);//Alexa's Rfunction/
  hist_man->SetRvalRec(r_val);

}


void dvcsEventAction::AddEdep( int block_No, double energy )
{
  if( block_No >=0 && block_No < 208 )
    {
      calo_edep_[block_No] = calo_edep_[block_No] + energy;
    }
  else
    {
      G4cout<<" Atention!! the calo index is out of range"<<G4endl;
    }
}
