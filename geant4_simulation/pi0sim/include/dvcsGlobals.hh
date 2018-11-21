#include "/work/halla/dvcs/disk1/salina/soft-2/TDVCSDB.h"

struct dvcsGlobals
{
  static int run_number;
  static double Ebeam; // Beam energy in GeV's
  static double HRS_angle; // central angle of HRS window
  static double HRS_momentum; // HRS momentum central value
  static double Calo_distance; // Calorimeter distance
  static double Calo_angle;    // Calorimeter angle
  static int    target_type;   // target type p or deuteron
  static int    target_gen_proc_type; // target type for determining process (like DVCS on p, n, or d)
  static double target_density; // density of target
  static double target_offset; // shift of the target center from the SC center
  static double target_length; // target lengt
  static bool hit_HRS_CALO_flag;
  static TDVCSDB *db; // pointer to DB
};
