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
  if ( fileIn.Filename().Base() == "conflib.dat" )
  {
    mprintf("  LMOD CONFLIB file\n");
    return true;
  }
  return false;
}

// Traj_Conflib::closeTraj()
void Traj_Conflib::closeTraj() {
  file_.CloseFile();
}

// Traj_Conflib::openTraj()
int Traj_Conflib::openTrajin() {
  return (file_.OpenFile());
}

// Traj_Conflib::setupTrajin()
int Traj_Conflib::setupTrajin(std::string const& fname, Topology* trajParm)
{
  int Frames;
  if (file_.OpenRead(fname)) return TRAJIN_ERR;
  size_t file_size = (size_t)file_.UncompressedSize();
  if (file_size > 0) {
    // Conflib frame is double,double,int,natom*3*double
    long unsigned int confFrame = (((trajParm->Natom() * 3) + 2) * sizeof(double)) + sizeof(int);
    Frames = (int) (file_size / confFrame);
    if ( (file_size % confFrame) != 0 ) {
      mprintf("Warning: %s: Could not accurately predict # frames. This can indicate either\n",
              file_.Filename().base());
      mprintf("Warning: the wrong topology is associated with this CONFLIB file or that the\n");
      mprintf("Warning: trajectory is corrupted. Will attempt to read %i frames.\n", Frames);
    }
  } else {
    Frames = TRAJIN_UNK;
  }
  conflibAtom_ = trajParm->Natom();
  return Frames;
}

// Traj_Conflib::readFrame()
int Traj_Conflib::readFrame(int set, double *X, double *V,double *box, double *T) {

  if (file_.Read(&energy_,sizeof(double)) < 0) return 1;
  file_.Read(&radGyr_,sizeof(double));
  file_.Read(&timesFound_,sizeof(int));
  file_.Read(X,sizeof(double)*conflibAtom_*3); 

  if (debug_>0) mprinterr("CONFLIB %10i: E=%10.4f RoG=%10.4f Found=%6i %12.4f %12.4f %12.4f\n",
                         set, energy_, radGyr_, timesFound_, X[0], X[1], X[2]);
  return 0;
}

// Traj_Conflib::setupTrajout()
int Traj_Conflib::setupTrajout(std::string const& fname, Topology* trajParm,
                               int NframesToWrite, bool append)
{
  mprinterr("Error: conflib writes not yet implemented.\n");
  return 1;
}

// Traj_Conflib::writeFrame()
int Traj_Conflib::writeFrame(int set, double *X, double *V,double *box, double T) {
  mprinterr("Error: conflib writes not yet implemented.\n");
  return 1;
}

// Traj_Conflib::info()
void Traj_Conflib::Info() {
  mprintf("is an LMOD conflib file");
}

