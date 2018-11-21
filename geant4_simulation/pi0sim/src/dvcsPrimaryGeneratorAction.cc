#include "dvcsPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "cmath"
#include <stdlib.h>
#include <time.h>
#include "dvcsEventAction.hh"
#include "dvcsGlobals.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

dvcsPrimaryGeneratorAction::dvcsPrimaryGeneratorAction( dvcsEventAction *event_action, G4int s1, G4int s2 ):
  vert_pos(0), L_em_init(0), L_em_scat(0), L_em_scat_v(0), L_final_phot1(0), L_final_phot2(0), L_final_prot(0)
{
  ev_numb = 0; // Initialize it to 0;
  G4int n_particles = 1;
  
  test_count1 = 0;
  test_count2 = 0;
  srand( time(NULL) );

  G4int seed1 = rand()%1000;
  G4int seed2 = rand()%1000;
  ev_act = event_action;
  
  //=================Generator Things==================================
  G4cout<<"********************************************** Kuku *************************************************"<<G4endl;
  G4cout<<" As a seeds for dvcs generator were taken values "<<s1<<" and "<<s2<<" This is in the GenAvtion Constructor"<<G4endl;
  G4cout<<"********************************************** Kuku *************************************************"<<G4endl;
  
  double Ebeam = dvcsGlobals::Ebeam;
  double HRS_angle = dvcsGlobals::HRS_angle;
  double HRS_mom = dvcsGlobals::HRS_momentum;
  double calo_angle = abs(dvcsGlobals::Calo_angle);
  double calo_dist = dvcsGlobals::Calo_distance;
  int targ_gen_proc_type = dvcsGlobals::target_gen_proc_type;
  double targ_dens = dvcsGlobals::target_density;
  double targ_offs = dvcsGlobals::target_offset;
  double targ_length = dvcsGlobals::target_length;

  L_em_scat_v = new TLorentzVector();

  //dvcsEventAction::Ebeam = Ebeam;
  
  cout<<"==========================================="<<endl;
  cout<<"Eb = "<<Ebeam<<"\t HRS_angle = "<<HRS_angle<<"  s1 = "<<s1<<"   s2 = "<<s2<<endl;
  cout<<"==========================================="<<endl;

  gEv = new TGenPi0(Ebeam, targ_gen_proc_type, s1, s2); // P target
  gEv->SetTargetParam(targ_length, targ_offs, targ_dens); // P target

  cout<<"===============TT EE MM PP OO RR AA RR RR YY======================="<<endl;
  cout<<"Ebeam = "<<Ebeam<<"  gen proc type = "<<targ_gen_proc_type<<"   targ_length = "<<targ_length<<"  targ_offs = "<<targ_offs<<"   targ_dens = "<<targ_dens<<endl;

  gEv->SetGeometry(HRS_angle, HRS_mom, calo_angle, calo_dist); // kin1  // HRS central angle, e- momentum, Calo angle, Calo dist in mm

  cout<<"HRS angle in radians = "<<HRS_angle<<"   HRS_mom in GeV/c = "<<HRS_mom<<"calo_angle in radians = "<<calo_angle<<"calo dist in mm = "<<calo_dist<<endl;

  particleGun = new G4ParticleGun(n_particles);
}

dvcsPrimaryGeneratorAction::~dvcsPrimaryGeneratorAction()
{
  delete particleGun;
  delete L_em_scat_v;
}

void dvcsPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  dvcsGlobals::hit_HRS_CALO_flag = false;
  // test_count1 = test_count1 + 1;
  // cout<<"test_count1 = "<<test_count1<<endl;
  //================Generator Things==============================

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

  G4int electronID = 11;
  G4int photonID = 22;

  bool hit_HRS = false;
  TVector3 beam_dir;
  
  gEv->GenerateVertex();//Generates the vertex
  gEv->ExtBrem();//Make external pre-vertex radiative corrections
  gEv->GenKin();//Computes scattered electron kinematics
  gEv->IntRCBef();//Internal radiative corrections (before vertex)
  
  //cout<<"Init electron E = "<<gEv->GetInitialElectron()->E()<<endl;

  //cout<<"Kuku"<<endl;
  if(gEv->ComputeElectron())//Compute the scattered electron 4-vector
    {
      //  cout<<"after Compute Electron++++++++++++++"<<endl;
      
      ev_numb = ev_numb + 1;
      
      *L_em_scat_v = *gEv->GetScatteredElectron();
      
      gEv->IntRCAft();//Internal radiative corrections (after vertex)
      
      TLorentzVector L_b(0, 0, dvcsGlobals::Ebeam, dvcsGlobals::Ebeam);
      TLorentzVector L_scat_el = *gEv->GetScatteredElectron();
      TLorentzVector L_init_phot = L_b - L_scat_el;
      double nu = L_init_phot.E();
      
      TVector3 k_v = TVector3(L_b.Vect());
      TVector3 kp_v = TVector3(L_scat_el.Vect());  
      
      double Q2 = 2*k_v.Mag()*kp_v.Mag()*(1 - cos(k_v.Angle(kp_v)) ); // Q2 = 2p1p2*(1 - cos(tehta_12))
      
      //cout<<"Q2 diff = "<<(L_init_phot.M2() + Q2)<<endl;
      double Mp = 0.938;
      double xB = Q2/(2*Mp*nu);
      // Double_t q3=TMath::Sqrt(Q2+TMath::Power(nu,2.));
      // Double_t q0primemax=0.5*Q2*(1.-xB)/(xB*(Mp+nu-q3));
      // Double_t q0primemin=0.5*Q2*(1.-xB)/(xB*(Mp+nu+q3));
      
      // Double_t tmax=-Q2-2.*q0primemax*(nu-q3);
      // Double_t tmin=-Q2-2.*q0primemin*(nu+q3);
	  
      gEv->Settmin( -3.); //Use this method to change tmin (default -2 GeV)
      gEv->Settmax(0.);  //Use this method to change tmax (default 0 GeV)
	  
      gEv->ComputePi0(); // \gamma* p -> pi0 p
      gEv->TwoBodyDecay(0.1349766,0.,0.); // pi0 mass
      
      gEv->ApplySpecVerAcc();//Rotates all final vectors around the beam axis
	  
      hit_HRS = false;
      
      //Check for spectro and Calo accep.
      if(gEv->HitsSpectro(gEv->GetScatteredElectron()) && gEv->HitsCalo(gEv->GetFinalPhoton1()) && gEv->HitsCalo(gEv->GetFinalPhoton2()) )
	{
	  //cout<<"Hited HRS and Calo ---------------------"<<endl;
	  
	  hit_HRS = true;
	  
	  L_final_phot1 = gEv->GetFinalPhoton1();
	  L_final_phot2 = gEv->GetFinalPhoton2();
	  
	  L_final_prot = gEv->GetFinalProton();
	  L_em_init = gEv->GetInitialElectron();
	  L_em_scat = gEv->GetScatteredElectron();
	  beam_dir = L_em_scat->Vect();
	  vert_pos = gEv->GetVertex();
	  
	  //G4cout<<"rad effecton on electron:"<<(L_em_scat_v->P() - L_em_scat->P())<<G4endl;
	  //G4cout<<"rad effecton on electron:"<<L_em_scat_v<<"     "<<L_em_scat<<G4endl;

	  //gEv->Print();
	      
	  if( hit_HRS )
	    {
	      dvcsGlobals::hit_HRS_CALO_flag = true;
	      //ev_act->DefineWeights(gEv->GetPSF(), gEv->XSecSum(), gEv->XSecDif()); // I think we donn't have this for pi0 yet
	      ev_act->DefineWeights(gEv->GetPSF(), 1., 1.);
	      //		  cout<<gEv->GetPSF()<<"\t"<<gEv->XSecSum()<<"\t"<<gEv->XSecDif()<<endl;
	      
	      ev_act->DefinePrimaries(L_em_init, L_em_scat, L_em_scat_v, L_final_phot1, L_final_phot2, L_final_prot, vert_pos->Z());
	      
	      particleGun->SetParticleDefinition(particleTable->FindParticle(electronID));
	      
	      G4ThreeVector ee(beam_dir.Px(),beam_dir.Py(),beam_dir.Pz());
	      particleGun->SetParticleEnergy(L_em_scat->E()*GeV);
	      particleGun->SetParticlePosition(G4ThreeVector(vert_pos->X()*cm, vert_pos->Y()*cm, vert_pos->Z()*cm));
	      particleGun->SetParticleMomentumDirection(ee);

	      particleGun->GeneratePrimaryVertex(anEvent);
	      
	      // test_count2 = test_count2 + 1;
	      // cout<<"test_count2------------------------------- = "<<test_count2<<endl;

	      particleGun->SetParticleDefinition(particleTable->FindParticle(photonID));
	      beam_dir = L_final_phot1->Vect();
	      G4ThreeVector ph1(beam_dir.Px(),beam_dir.Py(),beam_dir.Pz());
	      particleGun->SetParticleEnergy(L_final_phot1->E()*GeV);
	      particleGun->SetParticlePosition(G4ThreeVector(vert_pos->X()*cm, vert_pos->Y()*cm, vert_pos->Z()*cm));
	      particleGun->SetParticleMomentumDirection(ph1);
	      particleGun->GeneratePrimaryVertex(anEvent);

	      particleGun->SetParticleDefinition(particleTable->FindParticle(photonID));
	      beam_dir = L_final_phot2->Vect();
	      G4ThreeVector ph2(beam_dir.Px(),beam_dir.Py(),beam_dir.Pz());

	      particleGun->SetParticleEnergy(L_final_phot2->E()*GeV);
	      particleGun->SetParticlePosition(G4ThreeVector(vert_pos->X()*cm, vert_pos->Y()*cm, vert_pos->Z()*cm));
	      particleGun->SetParticleMomentumDirection(ph2);
	      particleGun->GeneratePrimaryVertex(anEvent);		  
	    }
	}
    }
  gEv->Clear();
  
}
