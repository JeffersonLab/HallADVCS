#!/bin/csh -f

setenv overwrite "yes"
if($1 == "keepmine") then
	setenv overwrite "no"
endif

# Get date from:
# git log -1
set CE_DATE = "(Thu March 31 2016)"

# Looking for user defined JLAB_VERSION
if( ! $?JLAB_VERSION || $overwrite == "yes") then
	setenv JLAB_VERSION 1.3
endif



# Software packages
#set packages = (clhep xercesc mysql qt geant4 root jlab12)
#### testing mysql version mismatch keep out old version
set packages = (clhep xercesc qt mysql geant4 root jlab12)                                                                                                                                                        

# Only print out if there's a prompt
alias echo 'if($?prompt) echo \!*  '


# Do not edit below this line
##################################################################################################
echo


# Looking for custom defined OSRELEASE
set DEFAULT_OSRELEASE = `osrelease.pl`
if($?OSRELEASE) then
	if($OSRELEASE != $DEFAULT_OSRELEASE) then
		echo " >> User defined OSRELEASE set to:"  $OSRELEASE
	endif
	else
        setenv OSRELEASE $DEFAULT_OSRELEASE
endif

# JLAB_SOFTWARE is where all the architecture software will be
#setenv JLAB_SOFTWARE /site/12gev_phys/2.1/Linux_CentOS7.2.1511-x86_64-gcc4.8.5
setenv JLAB_SOFTWARE /site/12gev_phys/1.3/Linux_CentOS7.2.1511-x86_64-gcc4.8.5
#setenv JLAB_SOFTWARE /lustre/expphy/volatile/clas/claseg3/parkkj/DVCS/G4/dvcs/jlab_software/1.2

echo
echo " > Common Environment Version: <"$JLAB_VERSION">  "$CE_DATE
echo " > Running as "`whoami` on `hostname`
echo " > OS Release:    "$OSRELEASE
echo " > JLAB_ROOT set to:     "$JLAB_ROOT

source $JLAB_ROOT/$JLAB_VERSION/ce/versions.env
if( -d $JLAB_SOFTWARE) then
	echo " > JLAB_SOFTWARE set to: "$JLAB_SOFTWARE
else
	mkdir -p $JLAB_SOFTWARE
	echo " > '$JLAB_SOFTWARE' is not a directory. Creating it."
endif
echo

# Sourcing packages. This will set the LD_LIBRARY_PATH. 
# Powerpcs and Darwins will be set right below
if( ! $?LD_LIBRARY_PATH) then
	setenv LD_LIBRARY_PATH ""
endif
if( ! $?PYTHONPATH) then
	setenv PYTHONPATH "."
endif

foreach p ($packages)
	source $JLAB_ROOT/$JLAB_VERSION/ce/$p".env"
#	source ./ce/$p".env"
end


# for powerpcs: LIBPATH
if ( $?LIBPATH ) then
  setenv LIBPATH ${LD_LIBRARY_PATH}:${LIBPATH}
else
  setenv LIBPATH ${LD_LIBRARY_PATH}
endif
# for Darwins systems: DYLD_LIBRARY_PATH
if ( $?DYLD_LIBRARY_PATH ) then
  setenv DYLD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${DYLD_LIBRARY_PATH}
else
  setenv DYLD_LIBRARY_PATH ${LD_LIBRARY_PATH}
endif

echo

unalias echo
