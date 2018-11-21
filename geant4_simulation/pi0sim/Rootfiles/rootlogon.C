//////////////////////////////////////////////////////////////////////////
//                                                                          
// rootlogon.C                                                            
//                                                                          
// Load Lib, paths and possible scripts to the analyzer upon start
//                                                                          
//////////////////////////////////////////////////////////////////////////
//	
// Author : Jin Huang <mailto:jinhuang@jlab.org>    Mar 2008
//          Chao Gu         Update for g2p use      Mar 2012
//          Kalyan Allada   Update for GMp expt.    Dec 2013           
//
//////////////////////////////////////////////////////////////////////////

{
  gSystem->Load("libDVCS.so");
  // gSystem->Load("/work/halla/dvcs/disk1/DVCS3/soft-2/libDVCS.so");
   //gSystem->Load("./F1/libF1.so");
  TString Arch(gSystem->GetBuildArch());
  TString Arch32("linux");
  TString Arch64("linuxx8664gcc");
  /*  if(Arch==Arch32){
    printf("\nrootlogon.C: Loading Replay Core Library..."); 
    gSystem->Load("./ReplayCore_C.so");
    gSystem->Load("./Fpp64/libFpp.so");
    gSystem->Load("./Gmp_Rastered_Beam/libGmp_Rastered_Beam.so");
    //    gSystem->Load("libFpp.so");
  }

  else if(Arch==Arch64){
    printf("\nrootlogon.C: Loading Replay Core Library..."); 
    gSystem->Load("./ReplayCore_C.so");
    gSystem->Load("./Fpp64/libFpp.so");
    gSystem->Load("./Gmp_Rastered_Beam/libGmp_Rastered_Beam.so");
    //gSystem->Load("libFpp64.so");
  }
  */
    //Load more libs here, if necessary. 
    //Make sure it's in path of $LD_LIBRARY_PATH

    printf("\nrootlogon.C: Adding include directories...");

    gSystem->AddIncludePath("-I$ANALYZER");
    gInterpreter->AddIncludePath("$ANALYZER/");

    gSystem->AddIncludePath("-I$ANALYZER/src");
    gInterpreter->AddIncludePath("$ANALYZER/src/");

    gSystem->AddIncludePath("-I$ANALYZER/hana_decode");
    gInterpreter->AddIncludePath("$ANALYZER/hana_decode/");

    gSystem->AddIncludePath("-I$ANALYZER/hana_scaler");
    gInterpreter->AddIncludePath("$ANALYZER/hana_scaler/");

    //gSystem->AddIncludePath("-I$LD_LIBRARY_PATH");
    //gInterpreter->AddIncludePath("$LD_LIBRARY_PATH");
    printf("\nrootlogon.C: Done!\n\n");
    //   gSystem->Load("libDVCS.so");
}
