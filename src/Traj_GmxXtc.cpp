#include "Traj_GmxXtc.h"
#include "CpptrajStdio.h"

#ifndef NO_XDRFILE
/// CONSTRUCTOR
Traj_GmxXtc::Traj_GmxXtc() : xd_(0), vec_(0), natoms_(0), prec_(1000) {}

/// DESTRUCTOR
Traj_GmxXtc::~Traj_GmxXtc() {
  closeTraj();
  if (vec_ != 0) delete[] vec_;
}

// Traj_GmxXtc::ID_TrajFormat()
bool Traj_GmxXtc::ID_TrajFormat(CpptrajFile& infile) {
  // See if we can read # atoms from xtc file.
  // NOTE: read_xtc_natoms *should* take const char* - boo
  if ( read_xtc_natoms( (char*)infile.Filename().full(), &natoms_ ) != exdrOK )
    return false;
  return (natoms_ > 0);
}

// Traj_GmxXtc::setupTrajin()
int Traj_GmxXtc::setupTrajin(FileName const& fnameIn, Topology* trajParm)
{
  if (fnameIn.empty()) return 1;
  fname_ = fnameIn;
  // Read number of atoms
  if ( read_xtc_natoms( (char*)fname_.full(), &natoms_ ) != exdrOK ) {
    mprinterr("Error: Could not get number of atoms from XTC file.\n");
    return TRAJIN_ERR;
  }
  if (natoms_ != trajParm->Natom()) {
    mprinterr("Error: # atoms in XTC file (%i) does not match # atoms in parm %s (%i)\n",
              natoms_, trajParm->c_str(), trajParm->Natom());
    return 1;
  }
  // Allocate arrays for reading coords
  vec_ = new rvec[ natoms_ ];
  if (vec_ == 0) return 1;
  // FIXME need to upgrade to libxdrfile2 to get # frames...
  if (openTrajin()) return TRAJIN_ERR;
  Frame tmp( natoms_ );
  if (readFrame(0, tmp)) return TRAJIN_ERR;
  closeTraj();
  SetCoordInfo( CoordinateInfo(ReplicaDimArray(), tmp.BoxCrd(), false, false, false, false) );
  return TRAJIN_UNK;
}

int Traj_GmxXtc::setupTrajout(FileName const& fname, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  return 1;
}

// Traj_GmxXtc::openTrajin()
int Traj_GmxXtc::openTrajin() {
  xd_ = xdrfile_open(fname_.full(), "r");
  if (xd_ == 0) {
    mprinterr("Error: Could not open XTC file for read.\n");
    return 1;
  }
  return 0;
}

// Traj_GmxXtc::closeTraj()
void Traj_GmxXtc::closeTraj() {
  if (xd_ != 0) xdrfile_close(xd_);
  xd_ = 0;
}

// Traj_GmxXtc::readFrame()
int Traj_GmxXtc::readFrame(int set, Frame& frameIn) {
  float time;
  int step;
  int result = read_xtc(xd_, natoms_, &step, &time, box_, vec_, &prec_);
  mprintf("DEBUG: set %i step %i time %f\n", set, step, time);
  int idx = 0;
  for (int ix = 0; ix < natoms_; ix++)
    for (int kx = 0; kx < DIM; kx++)
      frameIn[idx++] = (double)vec_[ix][kx];
  idx = 0;
  Matrix_3x3 ucell;
  for (int ii = 0; ii < DIM; ii++)
    for (int ij = 0; ij < DIM; ij++)
      ucell[idx++] = (double)box_[ii][ij];
  frameIn.SetBox( Box(ucell) );
  if (result != exdrOK) return 1;
  return 0;
}

int Traj_GmxXtc::writeFrame(int set, Frame const& frameOut) {
  return 1;
}

void Traj_GmxXtc::Info() {

}

int Traj_GmxXtc::readVelocity(int set, Frame& frameIn) {
  return 1;
}

int Traj_GmxXtc::readForce(int set, Frame& frameIn) {
  return 1;
}


#else
// =============================================================================
Traj_GmxXtc::Traj_GmxXtc() {}

Traj_GmxXtc::~Traj_GmxXtc() {}

bool Traj_GmxXtc::ID_TrajFormat(CpptrajFile&) { return false; }

int Traj_GmxXtc::setupTrajin(FileName const&, Topology*) { return 1; }

int Traj_GmxXtc::setupTrajout(FileName const&,Topology*,CoordinateInfo const&,int,bool) {return 1;}

int Traj_GmxXtc::openTrajin() { return 1; }

void Traj_GmxXtc::closeTraj() { }

int Traj_GmxXtc::readFrame(int, Frame&) { return 1; }

int Traj_GmxXtc::writeFrame(int, Frame const&) { return 1; }

void Traj_GmxXtc::Info() { }

int Traj_GmxXtc::readVelocity(int, Frame&) { return 1; }

int Traj_GmxXtc::readForce(int, Frame&) { return 1; }
#endif
