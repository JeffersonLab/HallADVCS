#!/usr/bin/perl
#use strict;
#use warnings;
#use Cwd;
#use File::Basename;
#to run, enter kin followed by the script name, e.g. perl sim_pi0_command.pl 484 
my $kinnumber = $ARGV[0];
my $pi0dir = '/work/halla/dvcs/disk1/salina/pi0sim/build';

#checking if the exec file for simulation (pi0_2017) exists
#while ($filename = readdir(DIR)) {
#if (-e "$pi0dir/pi0_2017"){
 #  print "$filename\n";
#   system("./pi0_2017 run.mac $kinnumber\n");
  # print $init."\n";
  #  }
if (-e "$pi0dir"){
   system(chdir($pi0dir)); 
   system("./pi0_2017 pi0_run.mac $kinnumber\n");
   }
#done


