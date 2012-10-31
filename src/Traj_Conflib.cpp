// Traj_Conflib
#include "Traj_Conflib.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_Conflib::Traj_Conflib() : 
  energy_(0),
  radGyr_(0),
  timesFound_(0),
  conflibAtom_(0)
{}

//------------------------------------------------------------------------
bool Traj_Conflib::ID_TrajFormat(CpptrajFile& fileIn) {
  // If the file name is conflib.dat, assume this is a conflib.dat file 
  // from LMOD. Cant think of a better way to detect this since there is no 
  // magic number but the file is binary.
  if ( FileName_ == "conflib.dat" )
  {
    mprintf("  LMOD CONFLIB file\n");
    return true;
  }
  return false;
}

// Traj_Conflib::closeTraj()
void Traj_Conflib::closeTraj() {
  CloseFile();
}

// Traj_Conflib::openTraj()
int Traj_Conflib::openTraj() {
  if (OpenFile()) return 1;
  return 0;
}

// Traj_Conflib::setupTrajin()
int Traj_Conflib::setupTrajin(std::string const& fname, Topology* trajParm,
                    TrajInfo& tinfo)
{
  long unsigned int confFrame;
  int Frames = 0;

  // Conflib is double,double,int,natom*3*double
  confFrame = (((trajParm->Natom() * 3) + 2) * sizeof(double)) + sizeof(int);
  Frames = (int) (file_size_ / confFrame);

  if ( (file_size_ % confFrame) != 0 ) {
    mprintf("Warning: %s: Could not accurately predict # frames. This usually \n",
            BaseFileStr());
    mprintf("         indicates a corrupted trajectory. Will attempt to read %i frames.\n",
            Frames);
  }
  conflibAtom_ = trajParm->Natom();
  return Frames;
}

// Traj_Conflib::readFrame()
int Traj_Conflib::readFrame(int set, double *X, double *V,double *box, double *T) {

  if (IO->Read(&energy_,sizeof(double)) < 0) return 1;
  IO->Read(&radGyr_,sizeof(double));
  IO->Read(&timesFound_,sizeof(int));
  IO->Read(X,sizeof(double)*conflibAtom_*3); 

  if (debug_>0) mprinterr("CONFLIB %10i: %10.4lf %10.4lf %6i %10.4lf %10.4lf %10.4lf\n",
                         set, energy_, radGyr_, timesFound_, X[0], X[1], X[2]);
  return 0;
}

// Traj_Conflib::setupTrajout()
int Traj_Conflib::setupTrajout(std::string const& fname, Topology* trajParm,
                     int NframesToWrite, TrajInfo const& tinfo)
{
  mprintf("Error: conflib writes not yet implemented.\n");
  return 1;
}

// Traj_Conflib::writeFrame()
int Traj_Conflib::writeFrame(int set, double *X, double *V,double *box, double T) {
  mprintf("Error: conflib writes not yet implemented.\n");
  return 1;
}

// Traj_Conflib::info()
void Traj_Conflib::info() {
  mprintf("is an LMOD conflib file");
}

