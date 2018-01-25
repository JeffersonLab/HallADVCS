{
  gSystem->Load("libDVCS.so");

  TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");

  Double_t *beamenergy=new Double_t[1]; // Beam energy
  Double_t *calodistance=new Double_t[1]; // Calorimeter distance
  Double_t *caloangle=new Double_t[1]; // Calorimeter angle
  Double_t *hrsangle=new Double_t[1]; // HRS angle
  Double_t *hrsmomentum=new Double_t[1]; // HRS momentum

  Int_t *targettype=new Int_t[1]; // Target type
  Int_t *genprocess=new Int_t[1]; // Simulation process

  /////////////////////////////////////////
  double vars[6] = 
    {
      4822,
      8.843,
      2,
      15.184,
      20.244,
      3.996
    };
  Int_t run=int(vars[0]); // Run number: IMPORTANT

  // Simulation parameter for this run number
  beamenergy[0]=vars[1];//GeV
  calodistance[0]=vars[2]*100.;//cm
  caloangle[0]=vars[3]*TMath::DegToRad();//rad (notice positive value)
  hrsangle[0]=vars[4]*TMath::DegToRad();//rad
  hrsmomentum[0]=vars[5];//GeV
  targettype[0]=0;//LH2 target
  genprocess[0]=0;//proton DVCS
  //////////////////////////////////////
  //////////////////////////////////////
  

  //Filling the database for this run number with these parameters
  db->AddEntry_d("BEAM_param_Energy",run,run,beamenergy,"Sample kinematics");
  db->AddEntry_d("CALO_geom_Dist",run,run,calodistance,"Sample kinematics");
  db->AddEntry_d("CALO_geom_Yaw",run,run,caloangle,"Sample kinematics");
  db->AddEntry_d("SIMU_param_HRSangle",run,run,hrsangle,"Sample kinematics");
  db->AddEntry_d("SIMU_param_HRSmomentum",run,run,hrsmomentum,"Sample kinematics");
  db->AddEntry_i("SIMU_param_TargetType",run,run,targettype,"Sample kinematics");
  db->AddEntry_i("SIMU_param_GenProcess",run,run,genprocess,"Sample kinematics");

  delete db;
}
