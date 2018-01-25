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
#include "/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/RFunction.C"
#include "/lustre/expphy/volatile/halla/dvcs/mongi/alexarfunction/RFunction_110717.C"

dvcsEventAction::dvcsEventAction( dvcsHistManager *aa, ObjManager *bb )
{
  hist_man = aa;
  int run_number = dvcsGlobals::run_number;
  gdvcs->SetWF(0);gdvcs->SetRun(run_number);gdvcs->ForceUpdate(); // This should organized by an automatic way !! The production is 8000
  
  rand2.SetSeed(0);
  
  cout<<"In Event Action Run number is "<<run_number<<endl;

  dvcs_event = bb->GetDVCSEvent_ptr();
  calo_event = bb->GetCaloEvent_ptr();
  L_calo_phot = bb->GetCaloPhot_ptr();//comment this line out for pion sim
  
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
      ev_number = evt->GetEventID();//these 5 lines different for pi0sim
      Px[4] = 0;
      Py[4] = 0;
      Pz[4] = 0;
      E[4] = 0;
      
      for( int i = 0; i < 208; i++ )
	{
	  calo_edep_[i] = 0;
	}
    }
}

void dvcsEventAction::EndOfEventAction(const G4Event* evt)
{
  if( dvcsGlobals::hit_HRS_CALO_flag )
    {
      double Ebeam = dvcsGlobals::Ebeam;
      double Mp = 0.938;
      
      Int_t ObjectNumber=TProcessID::GetObjectCount();
      for( int j = 0; j < 208; j++ )
	{
	  TCaloBlock *block=calo_event->AddBlock(j);// j is block number
	  //	  block->SetBlockEnergy(rand2.Gaus(calo_edep_[j], sqrt(calo_edep_[j]/175.)));
	  block->SetBlockEnergy(calo_edep_[j]);
	}
      //      dvcs_event->SetVertex(0, 0, smear_vertz); // Maybe I should use smeared vertex here instead of the generated one
      dvcs_event->SetVertex(0, 0, vert_z); // Maybe I should use smeared vertex here instead of the generated one
      calo_event->TriggerSim(0.2);
      calo_event->DoClustering();
      int n_clust = calo_event->GetNbClusters();
  
      // Here somehow I need to make an array of TLorentzVectors, instead of just one L_calo_phot
      *L_calo_phot = TLorentzVector(0, 0, 0, 0);
      for( int k = 0; k < n_clust; k++ )
	{
	  calo_event->GetCluster(k)->Analyze();
	  *L_calo_phot = dvcs_event->GetPhoton(k, 7, 0);
	  //  *L_calo_phot = 1.04*(*L_calo_phot);//Removed 4% photon energy correction
		  //G4cout<<"L_calo_phot.E="<<L_calo_phot->E()<<endl;
	}
      //  cout<<"n_clust = "<<n_clust<<endl;

      if( n_clust >= 1 && hit_HRS )
	{
	  double clust_x = calo_event->GetCluster(0)->GetX();
	  double clust_y = calo_event->GetCluster(0)->GetY();
	  double clust_E = calo_event->GetCluster(0)->GetE();
	  int clust_size = calo_event->GetCluster(0)->GetClusSize();
      
	  hist_man->SetCaloData(clust_E, clust_x, clust_y, clust_size);
      
	  // cout<<"clust_E = "<<clust_E<<endl;
	  // cout<<"Calo_phot_E = "<<(L_calo_phot->E())<<endl;

	  double phot_px = L_calo_phot->Px();
	  double phot_py = L_calo_phot->Py();
	  double phot_pz = L_calo_phot->Pz();
      
	  hist_man->SetPhotRec(phot_px, phot_py, phot_pz); //Set Reconstructed photon kinematics: Px, Py, Pz
      
	  //hist_man->Fill_Ntuple(vert_z, &Px[0], &Py[0], &Pz[0], &E[0], &calo_edep_[0], calo_event);
	  //hist_man->Fill_Ntuple(vert_z, &Px[0], &Py[0], &Pz[0], &E[0], &calo_edep_[0]);
      
	  TVector3 k(0, 0, Ebeam);

	  TVector3 kp(Px[4], Py[4], Pz[4]);
	  TVector3 qp(phot_px, phot_py, phot_pz);
      
	  //cout<<"Pz[4] = "<<Pz[4]<<endl;

	  double Q2 = 2*k.Mag()*kp.Mag()*(1 - cos(k.Angle(kp)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))

	  TVector3 q = k - kp;
      
	  TVector3 v1=q.Cross(kp);
	  TVector3 v2=q.Cross(qp);
	  double fphi=v1.Angle(v2);
	  if(q.Dot(v1.Cross(v2))<0) fphi=2.*TMath::Pi()-fphi;

	  double cos=(q.Dot(qp))/(q.Mag()*qp.Mag());
	  double nu=Ebeam-kp.Mag();
	  double tM=(Q2*Mp+2.*nu*Mp*(nu-sqrt(nu*nu+Q2)*cos))/(sqrt(nu*nu+Q2)*cos-nu-Mp);
	  double xB = Q2/(2*Mp*nu);

	  hist_man->SetkineRec(tM, xB, Q2, fphi);

	  TLorentzVector Lb, Lp, Lgp, Lkp, L_mis;
	  Lb.SetPxPyPzE(0, 0, Ebeam, Ebeam);
	  Lp.SetPxPyPzE(0, 0, 0, Mp);
	  Lkp.SetPxPyPzE(Px[4], Py[4], Pz[4], E[4]);
	  Lgp = *L_calo_phot;
      
	  //      cout<<"Lgp.M() = "<<(Lgp.M())<<"   Lgp.P() = "<<(Lgp.P())<<endl;
	  //cout<<"Lkp.M() = "<<(Lkp.M())<<"   Lkp.P() = "<<(Lkp.P())<<endl;

	  double mm2 = (Lb + Lp - Lgp - Lkp).M2();

	  hist_man->SetMM2(mm2);

	  hist_man->Fill_Ntuple();
	}

      hit_HRS = false;
      calo_event->Reset();
      TProcessID::SetObjectCount(ObjectNumber);
    }
}

void dvcsEventAction::DefinePrimaries(TLorentzVector *L_em_init, TLorentzVector *L_em_scat, 
				      TLorentzVector *L_phot, TLorentzVector *L_prot, double z)
{
  Px[0] = L_em_init->Px(); Px[1] = L_em_scat->Px();  Px[2] = L_phot->Px(); Px[3] = L_prot->Px();
  Py[0] = L_em_init->Py(); Py[1] = L_em_scat->Py();  Py[2] = L_phot->Py(); Py[3] = L_prot->Py();
  Pz[0] = L_em_init->Pz(); Pz[1] = L_em_scat->Pz();  Pz[2] = L_phot->Pz(); Pz[3] = L_prot->Pz();
  E[0]  = L_em_init->E();  E[1] = L_em_scat->E();    E[2]  = L_phot->E();  E[3]  = L_prot->E();

  double HRS_angle = dvcsGlobals::HRS_angle;
  double Ebeam = dvcsGlobals::Ebeam;
  vert_z = z;
  smear_vertz = vert_z + rand2.Gaus(0, 1.2/sin(HRS_angle))*0.1; // 0.1 is for mm->cm
  
  hist_man->SetVertexes(smear_vertz, vert_z );
  
  //cout<<"vertz = "<<vert_z<<"\t vert_z_v = "<<smear_vertz<<endl;

  hist_man->SetScatemGen(Px[1], Py[1], Pz[1]);
  hist_man->SetPhotGen(Px[2], Py[2], Pz[2]);
  double beam_x = 0;
  double beam_y = 0;
  double beam_z = Ebeam;
  hist_man->SetBeamGen(Px[0], Py[0], Pz[0]);
  hist_man->SetBeamRec(beam_x, beam_y, beam_z);

  //===========Calculate Q2, t, phi=========
  
  const double Mp = 0.938;
  
  TVector3 k_v(0, 0, E[0]);//Changed from Ebeam elog 12+GeV/473
  
  TVector3 kp_v(Px[1], Py[1], Pz[1]);
  TVector3 qp_v(Px[2], Py[2], Pz[2]);
  
  double Q2_v = 2*k_v.Mag()*kp_v.Mag()*(1 - cos(k_v.Angle(kp_v)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))
  
  TVector3 q_v = k_v - kp_v;
  
  TVector3 v1_v=q_v.Cross(kp_v);
  TVector3 v2_v=q_v.Cross(qp_v);
  double fphi_v=v1_v.Angle(v2_v);
  if(q_v.Dot(v1_v.Cross(v2_v))<0) fphi_v=2.*TMath::Pi()-fphi_v;
  
  double cos_v=(q_v.Dot(qp_v))/(q_v.Mag()*qp_v.Mag());
  double nu_v=E[0]-kp_v.Mag();//Changed from Ebeam elog 12+GeV/473
  double xB_v = Q2_v/(2*Mp*nu_v);
  //double tM_v=(Q2_v*Mp+2.*nu_v*Mp*(nu_v-sqrt(nu_v*nu_v+Q2_v)*cos_v))/(sqrt(nu_v*nu_v+Q2_v)*cos_v-nu_v-Mp);
  double tM_v = 2*Mp*(Mp - L_prot->E());  

  hist_man->SetkineGen(tM_v, xB_v, Q2_v, fphi_v);

  // TLorentzVector L_Eb_tmp;  L_Eb_tmp.SetPxPyPzE(0, 0, Ebeam, Ebeam);
  // TLorentzVector L_init_phot = L_Eb_tmp - *L_em_scat;

  // cout<<" t (using only angles) = "<<tM_v<<endl;
  // cout<<" t (using momenta) = "<<(*L_phot - L_init_phot).M2()<<endl;
  // cout<<" Their Differrence is  "<<(tM_v - (*L_phot - L_init_phot).M2())<<endl;
  



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
  Px[4] = mom.getX();
  Py[4] = mom.getY();
  Pz[4] = mom.getZ();
  E[4] = e;

  hist_man->SetScatemRec(Px[4], Py[4], Pz[4]);
  hit_HRS = true;

 //Calculate R_value

  double HRS_angle = dvcsGlobals::HRS_angle;
  Double_t pcentral= dvcsGlobals::HRS_momentum; // KINEMATICS 3
  int run_number = dvcsGlobals::run_number;
  int run=10555;//set default run to kin36_1 for Alexa's rfunction                                                                                          
  if(run_number==362){run=14150;}//Mongi had 14152
  else if(run_number==363){run=14480;}
  else if(run_number==481){run=12518;}
  else if(run_number==482){run=13009;}//check with alexa, changed from 13013
  else if(run_number==483){run=12843;}
  else if(run_number==484){run=13136;}
  else if(run_number==601){run=15017;}
  else if(run_number==602){run=14628;}

  cout<<"The run number is "<<run_number<<" and run is "<<run<<endl;

  //else{cout<<"ERROR : run number out of range, rfunction will be set to the default kinematic, i.e kin36_1 :"<<endl;}

  //Calculate R-value using electron momentum at HRS entrance

  //HRS_angle
  Double_t ry = -0.01*smear_vertz*TMath::Sin(TMath::ATan2(Px[4],Pz[4]));  // 0.01 to convert cm to meters
  Double_t rdp= (sqrt(Px[4]*Px[4]+Py[4]*Py[4]+Pz[4]*Pz[4])-pcentral)/pcentral;
  Double_t rtheta = TMath::ATan2(-Py[4],TMath::Sqrt(Px[4]*Px[4]+Pz[4]*Pz[4]));
  Double_t rphi   = TMath::ATan2(Px[4],Pz[4]) - HRS_angle;
  //R_function *rfunc_gula= new R_function(Form("/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/dis_cut_%d.root",run_number));
  //double r_val = rfunc_gula->Global_R_function(rphi,rdp,rtheta,ry);//dvcs_event->rfunction(ry,rdp,rtheta,rphi); 
  double r_val = RFunction(run,rtheta,rdp,rphi,ry);//Alexa's Rfunction/
  hist_man->SetRvalRec(r_val);

  //Calculate R-value using electron momentum production vertex
  //HRS_angle
  Double_t ry_v = -0.01*smear_vertz*TMath::Sin(TMath::ATan2(Px[1],Pz[1]));  // 0.01 to convert cm to meters
  Double_t rdp_v= (sqrt(Px[1]*Px[1]+Py[1]*Py[1]+Pz[1]*Pz[1])-pcentral)/pcentral;
  Double_t rtheta_v = TMath::ATan2(-Py[1],TMath::Sqrt(Px[1]*Px[1]+Pz[1]*Pz[1]));
  Double_t rphi_v   = TMath::ATan2(Px[1],Pz[1]) - HRS_angle;
  //R_function *rfunc_gula= new R_function(Form("/lustre/expphy/volatile/halla/dvcs/mongi/gularfunction/dis_cut_%d.root",run_number));
  //double r_val = rfunc_gula->Global_R_function(rphi,rdp,rtheta,ry);//dvcs_event->rfunction(ry,rdp,rtheta,rphi); 
  double r_val_v = RFunction(run,rtheta_v,rdp_v,rphi_v,ry_v);//Alexa's Rfunction/
  hist_man->SetRvalGen(r_val_v);

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
