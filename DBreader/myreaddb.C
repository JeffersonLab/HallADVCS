{
  gSystem->Load("libDVCS.so");

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=361; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=362; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=363; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=481; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=482; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=483; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=484; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=601; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;
  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  ////////////////////////////////
  Int_t run=603; // Run number
  ////////////////////////////////

  //Reading the database for this run number
  Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run);
  Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run);
  Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run);
  Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run);
  Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run);

  Int_t *targettype=db->GetEntry_i("SIMU_param_TargetType",run);
  Int_t *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run);

  cout<<"Run number : "<<run<<endl;

  cout<<"Beam energy : "<<beamenergy[0]<<" GeV"<<endl; 
  cout<<"Calorimeter distance : "<<calodistance[0]<<" cm"<<endl;
  cout<<"Calorimeter angle : "<<caloangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS angle : "<<hrsangle[0]*TMath::RadToDeg()<<" deg"<<endl; 
  cout<<"HRS momentum : "<<hrsmomentum[0]<<" GeV"<<endl; 
  cout<<"Target type : "<<targettype[0]<<endl;
  cout<<"Simulation process : "<<genprocess[0]<<endl;
  
  delete beamenergy;
  delete calodistance;
  delete caloangle;
  delete hrsangle;
  delete hrsmomentum;
  delete targettype;
  delete genprocess;

  delete db;
}
