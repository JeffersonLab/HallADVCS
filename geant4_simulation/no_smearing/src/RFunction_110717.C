#include "TMath.h"
#include "TF1.h"

//------------ Alexa Johnson --------------
//------------ alexa@jlab.org -------------
//----------------for DVCS ----------------
//---- Last updated: November 7, 2017   ------


Double_t RFunction(Int_t RunNumber, Double_t th, Double_t dp, Double_t ph, Double_t y){


  // Run number is the the run number assigned by Coda to the corresponding root file.
  // variable: th     is THETA in TARGET COORDINATES, root branch name = "L.tr.tg_th"
  // variable: dp     is DELTA in TARGET COORDINATES, root branch name = "L.tr.tg_dp"
  // variable: ph     is PHI in TARGET COORDINATES,   root branch name = "L.tr.tg_ph"
  // variable: y      is Y in TARGET COORDINATES,     root branch name = "L.tr.tg_y"

  // THIS DOES NOT CONTAIN ANY CUTS ON Z-VERTEX!


  Int_t kinparam_one;
  Int_t kinparam_two;

  //****************************************
  //****************************************
  //**  All runs prior to Spring 2016   ****
  //****************************************
  //****************************************

  if(RunNumber <=11217){
    kinparam_one = 0;
    kinparam_two = 0;
    // ideal r-cut = 0.003
  }

  //****************************************
  //****************************************
  //****        Spring 2016 runs        ****
  //****************************************
  //****************************************


  else if(RunNumber >= 12508 && RunNumber <=12647){
    kinparam_one = 48;
    kinparam_two = 1;
    // ideal r-cut = 0.003
  }
  else if((RunNumber >= 13000 && RunNumber <= 13009) || (RunNumber >= 13183 && RunNumber <= 13237)){
    kinparam_one = 48;
    kinparam_two = 2;
    //ideal r-cut = 0.003
  }
  else if(RunNumber >= 12838 && RunNumber <=12992){
    kinparam_one = 48;
    kinparam_two = 3;
    //ideal r-cut = 0.003
  }
  else if((RunNumber >= 13100 && RunNumber <=13162) || (RunNumber >= 13279 && RunNumber <= 13418)){
    kinparam_one = 48;
    kinparam_two = 4;
    //ideal r-cut = 0.004
  }


  //****************************************
  //****************************************
  //****        Fall 2016 runs        ******
  //****************************************
  //****************************************

  else if(RunNumber>=14150 && RunNumber<=14260){
    kinparam_one = 36;
    kinparam_two = 2;
    //ideal r-cut = 0.005
  }
  else if(RunNumber >=14476 && RunNumber<=14525){
    kinparam_one = 36;
    kinparam_two = 3;
    //ideal r-cut = 0.005;
  }
  else if(RunNumber >=14528 && RunNumber <= 14924){
    kinparam_one = 60;
    kinparam_two = 3;
    //ideal r-cut = 0.005;
  }
  else if(RunNumber >=15003 && RunNumber <= 15100){
    kinparam_one = 60;
    kinparam_two = 1;
    //ideal r-cut = 0.005;
  }



  //**************************************************

  else {
    printf("**********************************\n");
    printf("The run number you are using was not identified with a given kinematic setting.  Check that the run you are using is valid, and modify this script to include it with the approriate kinematic setting. \n");
    printf("**********************************\n");

 
  } 
  //***************************************************
  //***************************************************
  //***************************************************
  //            begin R-Function here
  //***************************************************
  //***************************************************
  //***************************************************

 
  Double_t r_thdp_1, r_thdp_2, r_thdp_3, r_phdp_1, r_phdp_2,r_phdp_3,r_phdp_4, r_phy_1, r_phy_2, r_phy_3, r_phy_4, r_thph_1, r_thph_2;
  // these are the names of functions defining initial cuts.  if a cut was deemed redundant, it is assigned a large value of 100000 to not affect the resulting R-Value.
  const int n_2D_cuts = 13;
  Double_t R;

  //**********************************************
  //**********************************************
  //------   Rfunction for x = 36 settings   -----
  //**********************************************
  //**********************************************
  if(kinparam_one == 36){
    if(kinparam_two==2){   
      //-----------Cuts from the theta delta plane ----------------                  
      //r_thdp_1 =(13.33*dp + 0.593 + th)/TMath::Sqrt(TMath::Power(13.33,2)+1);
      //r_thdp_2 = (-1.47*dp +0.097 - th)/TMath::Sqrt(TMath::Power(1.47,2)+1);
      //r_thdp_3 = (-4.25*dp +0.185 - th)/TMath::Sqrt(TMath::Power(4.25,2)+1);
      r_thdp_1 = -(-1.7*th*th - 0.06*th - 0.04 - dp)/TMath::Sqrt(TMath::Power((-3.4*th -0.06),2) + 1);
      r_thdp_2 = (-0.574*dp +0.097 - th)/TMath::Sqrt(TMath::Power(0.574,2)+1);
      r_thdp_3 = (-4.76*dp +0.219 - th)/TMath::Sqrt(TMath::Power(4.76,2)+1);
      //-----------Cuts from the phi delta plane ----------------                  
      //r_phdp_1 = (-0.075*dp + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.075,2)+1);
      r_phdp_1 = (-0.1*dp + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_2 = 100000;
      r_phdp_3 = 100000;
      r_phdp_4 = 100000;
      //-----------Cuts from the phi y plane ----------------                  
      //r_phy_1 = (0.175*y + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.175,2)+1);
      r_phy_1 = (0.22*y + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.22,2)+1);
      r_phy_2 = (-0.12*y + 0.028 - ph)/TMath::Sqrt(TMath::Power(0.12,2)+1);
      r_phy_3 = 100000;
      r_phy_4 = 100000;
      //-----------Cuts from the theta phi plane ----------------                  
      //r_thph_1 = th + 0.04;
      //r_thph_2 = 0.055 - th;
      r_thph_1 = th + 0.053;
      r_thph_2 = 0.06 - th;  
    }
    if(kinparam_two==3){  
      //-----------Cuts from the theta delta plane ----------------                  
      //r_thdp_1 =(13.33*dp + 0.593 + th)/TMath::Sqrt(TMath::Power(13.33,2)+1);
      //r_thdp_2 = (-1.47*dp +0.097 - th)/TMath::Sqrt(TMath::Power(1.47,2)+1);
      //r_thdp_3 = (-4.25*dp +0.185 - th)/TMath::Sqrt(TMath::Power(4.25,2)+1);
      r_thdp_1 = -(-1.7*th*th - 0.06*th - 0.04 - dp)/TMath::Sqrt(TMath::Power((-3.4*th -0.06),2) + 1);
      r_thdp_2 = (-0.574*dp +0.097 - th)/TMath::Sqrt(TMath::Power(0.574,2)+1);
      r_thdp_3 = (-4.76*dp +0.219 - th)/TMath::Sqrt(TMath::Power(4.76,2)+1);
      //-----------Cuts from the phi delta plane ----------------                  
      //r_phdp_1 = (-0.075*dp + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.075,2)+1);
      r_phdp_1 = (-0.1*dp + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_2 = 100000;
      r_phdp_3 = 100000;
      r_phdp_4 = 100000;
      //-----------Cuts from the phi y plane ----------------                  
      //r_phy_1 = (0.175*y + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.175,2)+1);
      r_phy_1 = (0.22*y + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.22,2)+1);
      r_phy_2 = (-0.12*y + 0.028 - ph)/TMath::Sqrt(TMath::Power(0.12,2)+1);
      r_phy_3 = 100000;
      r_phy_4 = 100000;
      //-----------Cuts from the theta phi plane ----------------                  
      //r_thph_1 = th + 0.04;
      //r_thph_2 = 0.055 - th;
      r_thph_1 = th + 0.053;
      r_thph_2 = 0.06 - th;  
    }
  }



  //**********************************************
  //**********************************************
  //------   Rfunction for x = 60 settings   -----
  //**********************************************
  //**********************************************


  if(kinparam_one == 60){
    if(kinparam_two==1){   
      //-----------Cuts from the theta delta plane ----------------                  
      //r_thdp_1 =(13.33*dp + 0.593 + th)/TMath::Sqrt(TMath::Power(13.33,2)+1);
      //r_thdp_2 = (-1.47*dp +0.097 - th)/TMath::Sqrt(TMath::Power(1.47,2)+1);
      //r_thdp_3 = (-4.25*dp +0.185 - th)/TMath::Sqrt(TMath::Power(4.25,2)+1);
      r_thdp_1 = -(-1.7*th*th - 0.06*th - 0.04 - dp)/TMath::Sqrt(TMath::Power((-3.4*th -0.06),2) + 1);
      r_thdp_2 = (-0.574*dp +0.097 - th)/TMath::Sqrt(TMath::Power(0.574,2)+1);
      r_thdp_3 = (-4.76*dp +0.219 - th)/TMath::Sqrt(TMath::Power(4.76,2)+1);
      //-----------Cuts from the phi delta plane ----------------                  
      //r_phdp_1 = (-0.075*dp + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.075,2)+1);
      r_phdp_1 = (-0.1*dp + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_2 = 100000;
      r_phdp_3 = 100000;
      r_phdp_4 = 100000;
      //-----------Cuts from the phi y plane ----------------                  
      //r_phy_1 = (0.175*y + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.175,2)+1);
      r_phy_1 = (0.22*y + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.22,2)+1);
      r_phy_2 = (-0.12*y + 0.028 - ph)/TMath::Sqrt(TMath::Power(0.12,2)+1);
      r_phy_3 = 100000;
      r_phy_4 = 100000;
      //-----------Cuts from the theta phi plane ----------------                  
      //r_thph_1 = th + 0.04;
      //r_thph_2 = 0.055 - th;
      r_thph_1 = th + 0.053;
      r_thph_2 = 0.06 - th;  
  
    }
    if(kinparam_two==3){  
      //-----------Cuts from the theta delta plane ----------------                  
      //r_thdp_1 =(13.33*dp + 0.593 + th)/TMath::Sqrt(TMath::Power(13.33,2)+1);
      //r_thdp_2 = (-1.47*dp +0.097 - th)/TMath::Sqrt(TMath::Power(1.47,2)+1);
      //r_thdp_3 = (-4.25*dp +0.185 - th)/TMath::Sqrt(TMath::Power(4.25,2)+1);
      r_thdp_1 = -(-1.7*th*th - 0.06*th - 0.04 - dp)/TMath::Sqrt(TMath::Power((-3.4*th -0.06),2) + 1);
      r_thdp_2 = (-0.574*dp +0.097 - th)/TMath::Sqrt(TMath::Power(0.574,2)+1);
      r_thdp_3 = (-4.76*dp +0.219 - th)/TMath::Sqrt(TMath::Power(4.76,2)+1);
      //-----------Cuts from the phi delta plane ----------------                  
      //r_phdp_1 = (-0.075*dp + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.075,2)+1);
      r_phdp_1 = (-0.1*dp + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_2 = 100000;
      r_phdp_3 = 100000;
      r_phdp_4 = 100000;
      //-----------Cuts from the phi y plane ----------------                  
      //r_phy_1 = (0.175*y + 0.0275 + ph)/TMath::Sqrt(TMath::Power(0.175,2)+1);
      r_phy_1 = (0.22*y + 0.029 + ph)/TMath::Sqrt(TMath::Power(0.22,2)+1);
      r_phy_2 = (-0.12*y + 0.028 - ph)/TMath::Sqrt(TMath::Power(0.12,2)+1);
      r_phy_3 = 100000;
      r_phy_4 = 100000;
      //-----------Cuts from the theta phi plane ----------------                  
      //r_thph_1 = th + 0.04;
      //r_thph_2 = 0.055 - th;
      r_thph_1 = th + 0.053;
      r_thph_2 = 0.06 - th;  
    }
  }


  //**********************************************
  //**********************************************
  //------   Rfunction for x = 48 settings   -----
  //**********************************************
  //**********************************************
 
  if(kinparam_one == 48){
    if(kinparam_two==1){
      //-----------Cuts from the theta delta plane ----------------           
      r_thdp_1 =-(-13.3*dp -0.55 - th)/TMath::Sqrt(TMath::Power(13.3,2)+1);
      r_thdp_2 = (-4*dp +0.18 - th)/TMath::Sqrt(TMath::Power(4,2)+1);
      //r_thdp_3 = (-0.85*dp +0.07 - th)/TMath::Sqrt(TMath::Power(0.85,2)+1);
      r_thdp_3 = (-0.6*dp +0.07 - th)/TMath::Sqrt(TMath::Power(0.6,2)+1);
     
      //-----------------------------------------------
      //----------- Cuts from the phi delta plane ---------------
      r_phdp_1 = (-0.15*dp + 0.03 - ph)/TMath::Sqrt(TMath::Power(0.15,2)+1);
      r_phdp_2 =- (0.125*dp - 0.03 - ph)/TMath::Sqrt(TMath::Power(0.125,2)+1);
      r_phdp_3 = -dp + 0.05;
      r_phdp_4 = dp + 0.045;
      //-----------------------------
      //---------- Cuts from the phi y plane -----------
      //r_phy_1 = (-0.25*y + 0.025 - ph)/TMath::Sqrt(TMath::Power(0.25,2)+1);
      r_phy_1 = (-0.25*y + 0.03 - ph)/TMath::Sqrt(TMath::Power(0.25,2)+1);
      r_phy_2 = (0.5775*y + 0.042 - ph)/TMath::Sqrt(TMath::Power(0.5775,2)+1);
      r_phy_3 =- (0.538*y - 0.048 - ph)/TMath::Sqrt(TMath::Power(0.538,2)+1);
      r_phy_4 =- (-0.225*y - 0.03 - ph)/TMath::Sqrt(TMath::Power(0.225,2)+1);

      //r_thph_1 = (0.063 - th);
      //r_thph_2 = -(-0.053 - th);
      r_thph_1 = (0.06 - th);
      r_thph_2 = -(-0.05 - th); 
            
    }
    if(kinparam_two==2){
      //-----------Cuts from the theta delta plane ----------------           
      //r_thdp_1 =-(-4*dp -0.18 - th)/TMath::Sqrt(TMath::Power(4,2)+1);
      //r_thdp_2 = (-1.82*dp +0.08 - th)/TMath::Sqrt(TMath::Power(1.82,2)+1);
      r_thdp_1 =-(-4*dp -0.17 - th)/TMath::Sqrt(TMath::Power(4,2)+1);
      r_thdp_2 = (-1.9*dp +0.09 - th)/TMath::Sqrt(TMath::Power(1.9,2)+1);
      r_thdp_3 = (-0.4*dp +0.04 - th)/TMath::Sqrt(TMath::Power(0.4,2)+1);
     //-----------------------------------------------
      //----------- Cuts from the phi delta plane ---------------
      r_phdp_1 = 100000; // redundant
      //r_phdp_2 =- (0.1*dp - 0.04 - ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_2 =- (0.1*dp - 0.05 - ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_3 = -dp + 0.05;
      r_phdp_4 = 1000000;
      //------------------------------------------------
      //---------- Cuts from the phi y plane -----------
      // r_phy_1 = (-0.275*y + 0.04 - ph)/TMath::Sqrt(TMath::Power(0.275,2)+1);
      r_phy_1 = (-0.923*y + 0.0372 - ph)/TMath::Sqrt(TMath::Power(0.923,2)+1);
      r_phy_2 = 100000;
      r_phy_3 = 100000; 
      //r_phy_4 =- (-0.325*y - 0.038 - ph)/TMath::Sqrt(TMath::Power(0.325,2)+1);
      //r_phy_4 =- (-0.975*y - 0.038 - ph)/TMath::Sqrt(TMath::Power(0.975,2)+1);
      r_phy_4 =- (-0.975*y - 0.04 - ph)/TMath::Sqrt(TMath::Power(0.975,2)+1);


      r_thph_1 = (0.0325 - th);
      r_thph_2 = -(-0.025 - th);     
    }
    if(kinparam_two==3){
      //-----------Cuts from the theta delta plane ----------------           
      //r_thdp_1 =-(-9*dp -0.3545 - th)/TMath::Sqrt(TMath::Power(9,2)+1);
      r_thdp_2 = (-2.96*dp +0.137 - th)/TMath::Sqrt(TMath::Power(2.96,2)+1);
      //r_thdp_3 = (-0.545*dp +0.065 - th)/TMath::Sqrt(TMath::Power(0.545,2)+1);
      r_thdp_3 = (-0.492*dp +0.054 - th)/TMath::Sqrt(TMath::Power(0.492,2)+1);
      r_thdp_1 =-(-7.14*dp -0.293 - th)/TMath::Sqrt(TMath::Power(7.14,2)+1);
     //-----------------------------------------------
      //----------- Cuts from the phi delta plane ---------------
      r_phdp_1 = 100000;
      //r_phdp_2 =- (0.1075*dp - 0.038 - ph)/TMath::Sqrt(TMath::Power(0.1075,2)+1);
      r_phdp_2 =- (0.1075*dp - 0.045 - ph)/TMath::Sqrt(TMath::Power(0.1075,2)+1);	    
      r_phdp_3 = -dp + 0.05;
      r_phdp_4 = 100000;
      //------------------------------------------------
      //---------- Cuts from the phi y plane -----------
      //r_phy_1 =(-0.225*y + 0.03 - ph)/TMath::Sqrt(TMath::Power(0.225,2)+1);
      //r_phy_2 = (0.972*y + 0.0572 - ph)/TMath::Sqrt(TMath::Power(0.972,2)+1);
      //r_phy_3 = -(0.025*y - 0.0375 - ph)/TMath::Sqrt(TMath::Power(0.025,2)+1);
      //r_phy_4 =- (-0.325*y - 0.038 - ph)/TMath::Sqrt(TMath::Power(0.325,2)+1);
      //r_phy_1 =(-1.11*y + 0.032 - ph)/TMath::Sqrt(TMath::Power(1.11,2)+1);
      //r_phy_2 = -(-1.23*y - 0.039 - ph)/TMath::Sqrt(TMath::Power(1.23,2)+1);
      r_phy_1 =(-0.69*y + 0.034 - ph)/TMath::Sqrt(TMath::Power(0.69,2)+1);
      r_phy_2 = -(-0.71*y - 0.036 - ph)/TMath::Sqrt(TMath::Power(0.71,2)+1);
      r_phy_3 = 10000;
      r_phy_4 = 10000;

      r_thph_1 = (0.05 - th);
      //r_thph_2 = -(-0.03 - th);  
      r_thph_2 = -(-0.036 - th);  
   
    }
    if(kinparam_two==4){    
      //-----------Cuts from the theta delta plane ----------------           
      // r_thdp_1 =-(-4*dp -0.18 - th)/TMath::Sqrt(TMath::Power(4,2)+1);
      r_thdp_2 = (-2.1*dp +0.1 - th)/TMath::Sqrt(TMath::Power(2.1,2)+1);
      // r_thdp_3 = (-0.4*dp +0.05 - th)/TMath::Sqrt(TMath::Power(0.4,2)+1);
      r_thdp_1 =-(-4.5*dp -0.19 - th)/TMath::Sqrt(TMath::Power(4.5,2)+1);
      r_thdp_3 = (-0.45*dp +0.045 - th)/TMath::Sqrt(TMath::Power(0.45,2)+1);
     //-----------------------------------------------
      //----------- Cuts from the phi delta plane ---------------
      r_phdp_1 = 100000;
      //r_phdp_2 =- (0.1*dp - 0.04 - ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_2 =- (0.1*dp - 0.045 - ph)/TMath::Sqrt(TMath::Power(0.1,2)+1);
      r_phdp_3 = -dp + 0.05;
      r_phdp_4 = 100000;
      //------------------------------------------------
      //---------- Cuts from the phi y plane -----------
      //r_phy_1 = (-0.225*y + 0.0325 - ph)/TMath::Sqrt(TMath::Power(0.225,2)+1);
      r_phy_1 =(-0.69*y + 0.034 - ph)/TMath::Sqrt(TMath::Power(0.69,2)+1);
      r_phy_2 = -(-0.71*y - 0.036 - ph)/TMath::Sqrt(TMath::Power(0.71,2)+1);
      r_phy_3 = 100000;
      r_phy_4 = 100000;
      //r_phy_4 =- (-0.2375*y - 0.036 - ph)/TMath::Sqrt(TMath::Power(0.2375,2)+1);

      //r_thph_1 = (0.045 - th);
      r_thph_2 = -(-0.03 - th);
      r_thph_1 = (0.038 - th);

    }
  }



  //**********************************************
  //**********************************************
  // Rfunction for runs taken prior to Spring 2016  
  //**********************************************
  //**********************************************

  if(kinparam_one == 0){
    if(kinparam_two==0){
      //-----------Cuts from the theta delta plane ----------------           
      r_thdp_1 =-(-13.3*dp -0.55 - th)/TMath::Sqrt(TMath::Power(13.3,2)+1);
      r_thdp_2 = (-4*dp +0.18 - th)/TMath::Sqrt(TMath::Power(4,2)+1);
      //r_thdp_3 = (-0.85*dp +0.07 - th)/TMath::Sqrt(TMath::Power(0.85,2)+1);
      r_thdp_3 = (-0.6*dp +0.07 - th)/TMath::Sqrt(TMath::Power(0.6,2)+1);
     
      //-----------------------------------------------
      //----------- Cuts from the phi delta plane ---------------
      r_phdp_1 = (-0.15*dp + 0.03 - ph)/TMath::Sqrt(TMath::Power(0.15,2)+1);
      r_phdp_2 =- (0.125*dp - 0.03 - ph)/TMath::Sqrt(TMath::Power(0.125,2)+1);
      r_phdp_3 = -dp + 0.05;
      r_phdp_4 = dp + 0.045;
      //-----------------------------
      //---------- Cuts from the phi y plane -----------
      //r_phy_1 = (-0.25*y + 0.025 - ph)/TMath::Sqrt(TMath::Power(0.25,2)+1);
      r_phy_1 = (-0.25*y + 0.03 - ph)/TMath::Sqrt(TMath::Power(0.25,2)+1);
      r_phy_2 = (0.5775*y + 0.042 - ph)/TMath::Sqrt(TMath::Power(0.5775,2)+1);
      r_phy_3 =- (0.538*y - 0.048 - ph)/TMath::Sqrt(TMath::Power(0.538,2)+1);
      r_phy_4 =- (-0.225*y - 0.03 - ph)/TMath::Sqrt(TMath::Power(0.225,2)+1);

      //r_thph_1 = (0.063 - th);
      //r_thph_2 = -(-0.053 - th);
      r_thph_1 = (0.06 - th);
      r_thph_2 = -(-0.05 - th); 

            
    }
  }

  //**********************************************
  //**********************************************
  //------   Find min value = R value        -----



  Double_t R_Array[13] = {r_thdp_1,r_thdp_2,r_thdp_3,r_phdp_1,r_phdp_2,r_phdp_3,r_phdp_4,r_phy_1,r_phy_2,r_phy_3,r_phy_4,r_thph_1,r_thph_2};

  R = TMath::MinElement(n_2D_cuts,R_Array);

  return R;
}





