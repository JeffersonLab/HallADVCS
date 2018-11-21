#!/bin/tcsh
# source this script to set up the G4 DVCS12
setenv JLAB_VERSION 1.3

#default installation on ifarm and farm nodes, change it to use different installation
#setenv JLAB_ROOT /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/jlab_software
setenv JLAB_ROOT /site/12gev_phys
setenv ROOTSYS /apps/root/5.34.36/root
setenv DVCS_GEN /work/halla/dvcs/disk1/salina/dvcs_gen
#/volatile/halla/dvcs/wmhenry/geant4/workingOS7/dvcs_gen_CentOS7
#no need to modify variables below generally
#setup basic enviroment variables 
source jlabDVCS12.csh      
setenv ANALYZER /w/halla-scifs17exp/dvcs/disk1/bill/geant4/workingOS7/analyzer 
setenv DVCS_SOFT /work/halla/dvcs/disk1/salina/soft-sfa

#setenv DB_DIR /work/halla/dvcs/disk1/dvcs/onlana/DB
setenv DB_DIR /work/halla/dvcs/disk1/salina/DB
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ANALYZER}:${ROOTSYS}/lib:${DVCS_SOFT}:/work/halla/dvcs/disk1/salina/soft-2
setenv LD_LIBRARY_PATH ${ANALYZER}:${ROOTSYS}/lib:${DVCS_SOFT}:${LD_LIBRARY_PATH}
setenv PATH ${ROOTSYS}/bin:${ANALYZER}:${PATH}

