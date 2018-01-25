#!/bin/tcsh
# source this script to set up the G4 DVCS12
setenv JLAB_VERSION 1.3

#default installation on ifarm and farm nodes, change it to use different installation
#setenv JLAB_ROOT /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/jlab_software
setenv JLAB_ROOT /site/12gev_phys
setenv ROOTSYS /apps/root/5.34.36/root
#setenv DVCS_GEN /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/dvcs_gen_rafopar
#setenv DVCS_GEN /u/home/wmhenry/geant4/kijunpark/workingOS7/dvcs_gen_CentOS7
setenv DVCS_GEN /w/halla-scifs17exp/dvcs/disk1/bill/geant4/workingOS7/dvcs_gen_CentOS7
#setenv DVCS_GEN /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/dvcs_gen_CentOS7
#no need to modify variables below generally
#source $JLAB_ROOT/$JLAB_VERSION/ce/jlab.csh       #setup basic enviroment variables 
#source /u/home/wmhenry/geant4/kijunpark/workingOS7/jlabDVCS12.csh       #setup basic enviroment variables
source /w/halla-scifs17exp/dvcs/disk1/bill/geant4/workingOS7/jlabDVCS12.csh       #setup basic enviroment variables
#setenv ANALYZER /work/halla/dvcs/disk1/carlos/64bits/analyzer
#setenv ANALYZER  /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/my64bits/analyzer
#setenv ANALYZER /u/home/wmhenry/geant4/kijunpark/workingOS7/analyzer
setenv ANALYZER /w/halla-scifs17exp/dvcs/disk1/bill/geant4/workingOS7/analyzer 
#setenv DVCS_SOFT /work/halla/dvcs/disk1/carlos/64bits/soft-2
#setenv DVCS_SOFT  /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/my64bits/soft-2
#setenv DVCS_SOFT /u/home/wmhenry/geant4/kijunpark/workingOS7/soft-2
setenv DVCS_SOFT /w/halla-scifs17exp/dvcs/disk1/bill/geant4/workingOS7/soft-2
setenv DB_DIR /work/halla/dvcs/disk1/dvcs/onlana/DB
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ANALYZER}:${ROOTSYS}/lib:${DVCS_SOFT}:/work/halla/dvcs/disk4/rafopar/dvcs_gen
setenv PATH ${ROOTSYS}/bin:${ANALYZER}:${PATH}

