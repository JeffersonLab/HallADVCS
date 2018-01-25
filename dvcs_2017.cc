#include "G4UImanager.hh"
#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//#include "G4UIterminal.hh"
//#include "G4UItcsh.hh"

#include "dvcsSteppingAction.hh"
#include "dvcsEventAction.hh"
#include "dvcsRunAction.hh"
#include "dvcsPhysicsList.hh"
#include "dvcsDetectorConstruction.hh"
#include "dvcsPrimaryGeneratorAction.hh"
#include "dvcsHist_Manager.hh"
#include "TDVCSEvent.h"
#include "TCaloEvent.h"
#include "dvcsObj_Manager.hh"
#include "dvcsGlobals.hh"
#include "TDVCSDB.h"
#include <fstream>

using namespace std;

int    dvcsGlobals::run_number = 0;
double dvcsGlobals::Ebeam = 0;
double dvcsGlobals::HRS_angle = 0; 
double dvcsGlobals::HRS_momentum = 0; 
double dvcsGlobals::Calo_distance = 0;
double dvcsGlobals::Calo_angle = 0; 
int    dvcsGlobals::target_type = 0;
int    dvcsGlobals::target_gen_proc_type = 0;
double dvcsGlobals::target_density = 0;
double dvcsGlobals::target_offset = 0;
double dvcsGlobals::target_length = 0;
bool dvcsGlobals::hit_HRS_CALO_flag = false;

int main( int argc, char** argv)
{
  //==============================Initialization of members of the dvcsGlobals classs==================================
  
  if( argc == 3 || argc == 4 )
    {
      TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"munoz","");
      
      int run_number = atoi(argv[2]);

      TCaloEvent *calo_test_event = new TCaloEvent(run_number); // This is just for initializing calo with this run number
      
      //Reading the database for this run number
      Double_t *beamenergy=db->GetEntry_d("BEAM_param_Energy",run_number);
      Double_t *calodistance=db->GetEntry_d("CALO_geom_Dist",run_number);
      Double_t *caloangle=db->GetEntry_d("CALO_geom_Yaw",run_number);
      Double_t *hrsangle=db->GetEntry_d("SIMU_param_HRSangle",run_number);
      Double_t *hrsmomentum=db->GetEntry_d("SIMU_param_HRSmomentum",run_number);
      
      int *targettype=db->GetEntry_i("SIMU_param_TargetType",run_number);
      int *genprocess=db->GetEntry_i("SIMU_param_GenProcess",run_number);
      
      dvcsGlobals::run_number = run_number;
      dvcsGlobals::Ebeam = beamenergy[0];
      dvcsGlobals::Calo_distance = calodistance[0]*10; // in DB calo_distance is in cm therefre the factor 10 is arised
      dvcsGlobals::Calo_angle = -caloangle[0];
      dvcsGlobals::HRS_angle = hrsangle[0];
      dvcsGlobals::HRS_momentum = hrsmomentum[0];
      dvcsGlobals::target_type = targettype[0];
      dvcsGlobals::target_gen_proc_type = genprocess[0];

      int targ_type = dvcsGlobals::target_type;
      cout<<"targ_type = "<<targ_type<<endl;
      if( targ_type == 0 ) { dvcsGlobals::target_density = 0.0723; } // Hydrogen
      else if( targ_type == 1 ) { dvcsGlobals::target_density = 0.167; }  // Deuteron
      else
      {
	cout<<"Target type from DB is "<<targ_type<<endl;
	cout<<"Please check the target type in the DB. It should be 0 (proton) or 1 (deutron)"<<endl;
	cout<<"pragram is exiting..."<<endl;
	exit(11);
      }
      int target_gen_proc_type = dvcsGlobals::target_gen_proc_type;
      if( target_gen_proc_type == 0 ) {cout<<"Process to be genarated is DVCS on proton"<<endl;}
      else if ( target_gen_proc_type == 1 ) {cout<<"Process to be genarated is DVCS on neutron"<<endl;}
      else if ( target_gen_proc_type == 2 ) {cout<<"Process to be genarated is DVCS on deutron"<<endl;}
      else
	{
	  cout<<"process type in DB is "<<target_gen_proc_type<<endl;
	  cout<<"Please check it, it should be 0, 1 or 2, for DVCS on p, n or d"<<endl;
	}
      
      delete beamenergy;
      delete calodistance;
      delete caloangle;
      delete hrsangle;
      delete hrsmomentum;
      delete targettype;
      delete genprocess;
      delete db;
      delete calo_test_event;
      
    }
  else
    {
      cout<<"Please check number of arguments. Seems it is not correct"<<endl;
      cout<<"you should run the script with the following command:"<<endl;
      cout<<"dvcs dvcs_run.mac run_number"<<endl;
      cout<<"where:"<<endl;
      cout<<"dvcs - executable"<<endl;
      cout<<"dvcs_run_mac - a script which specifes number of events and few other things needed by GEANT4"<<endl;
      cout<<"run number - the run number in DB, where Calo angle, HRS_angle and these kind of variables are specified"<<endl;
      cout<<"------------------------------------------------"<<endl;
      cout<<"Program is exiting..."<<endl;
      // exit(1);
    }

  dvcsGlobals::target_offset = 0.;
  dvcsGlobals::target_length = 15;
  
  ofstream out_dat("dvcs_run.log");
  
  out_dat<<"-------------------------Run Settings-------------------------"<<endl;
  out_dat<<"Run number = "<<dvcsGlobals::run_number<<endl;
  out_dat<<"Beam energy in GeV = "<<dvcsGlobals::Ebeam<<endl;
  out_dat<<"Calo distance in mm = "<<dvcsGlobals::Calo_distance<<endl;
  out_dat<<"Calo angle in deg = "<<dvcsGlobals::Calo_angle<<endl;
  out_dat<<"HRS angle in deg = "<<dvcsGlobals::HRS_angle<<endl;
  out_dat<<"HRS Central momentum in GeV/c = "<<dvcsGlobals::dvcsGlobals::HRS_momentum<<endl;
  out_dat<<"Target type = "<<dvcsGlobals::target_type<<endl;
  out_dat<<"Target offset in cm = "<<dvcsGlobals::target_offset<<endl;
  out_dat<<"Target length in cm = "<<dvcsGlobals::target_length<<endl;
  out_dat<<"-----------------------End ofRun Settings---------------------"<<endl;
  //=========================End of Initialization of members of the dvcsGlobals classs==============================


  ObjManager *obj_man = new ObjManager();
  
  //TDVCSEvent *dvcs_event = new TDVCSEvent();
  //TCaloEvent *Calo_event = new TCaloEvent();

  G4RunManager *runManager = new G4RunManager;
  G4VisManager* visManager = new G4VisExecutive;

  dvcsDetectorConstruction *detector = new dvcsDetectorConstruction;
  runManager->SetUserInitialization(detector);
  
  G4VUserPhysicsList *phys_list = new dvcsPhysicsList;
  runManager->SetUserInitialization(phys_list);
  
  //dvcsHist_Manager *hist = new dvcsHist_Manager( obj_man );
  dvcsHistManager *hist = new dvcsHistManager();
  G4cout<<"dvcs_HistMan="<<hist<<G4endl;

  dvcsEventAction *event_action = new dvcsEventAction( hist, obj_man );
  runManager->SetUserAction(event_action);
  
  G4UserRunAction* run_action = new dvcsRunAction( hist );
  runManager->SetUserAction(run_action);

  G4int seed1 = 11;
  G4int seed2 = 12;
  
  if( argc == 4 )
    {
      G4cout<<"********************************************** Kuku *************************************************"<<G4endl;
      seed1 = atoi(argv[2]);
      seed2 = atoi(argv[3]);
      G4cout<<" As a seeds for dvcs generator were taken values "<<seed1<<" and "<<seed2<<G4endl;
      G4cout<<"********************************************** Kuku *************************************************"<<G4endl;
    }
  G4VUserPrimaryGeneratorAction *gen_action = new dvcsPrimaryGeneratorAction(event_action, seed1, seed2);
  runManager->SetUserAction(gen_action);
  
  dvcsSteppingAction *step_action = new dvcsSteppingAction(detector, event_action);
  runManager->SetUserAction(step_action);
  
  runManager->Initialize();

  /*  
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  */
  visManager->Initialize();
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc == 3)   // batch mode  
    {

      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);

// #ifdef G4UI_USE
//       G4UIExecutive * ui = new G4UIExecutive(argc,argv);
// #ifdef G4VIS_USE
      
// #endif
//       //      ui->SessionStart();
//       delete ui;
// #endif

    }
  else           // interactive mode : define UI session
    { 
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute dvcs_vis.mac");
#endif
      ui->SessionStart();
      delete ui;
#endif
      
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }
  
  
  
  //G4int numberOfEvents = 1000;
  //runManager->BeamOn(numberOfEvents);
  delete visManager;

  delete hist;
  delete runManager;
  return 0;
}
