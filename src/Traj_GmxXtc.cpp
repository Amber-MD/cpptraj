#include "Traj_GmxXtc.h"
#include "CpptrajStdio.h"
#include "Constants.h"

#ifndef NO_XDRFILE
/// CONSTRUCTOR
Traj_GmxXtc::Traj_GmxXtc() : xd_(0), vec_(0), dt_(1.0), natoms_(0), prec_(1000) {}

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
  if (fnameIn.empty()) return TRAJIN_ERR;
  fname_ = fnameIn;
  // Read number of atoms
  if ( read_xtc_natoms( (char*)fname_.full(), &natoms_ ) != exdrOK ) {
    mprinterr("Error: Could not get number of atoms from XTC file.\n");
    return TRAJIN_ERR;
  }
  if (natoms_ != trajParm->Natom()) {
    mprinterr("Error: # atoms in XTC file (%i) does not match # atoms in parm %s (%i)\n",
              natoms_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Allocate array for reading coords
  if (vec_ != 0) delete[] vec_;
  vec_ = new rvec[ natoms_ ];
  if (vec_ == 0) return TRAJIN_ERR;
  // FIXME need to upgrade to libxdrfile2 to get # frames...
  if (openTrajin()) return TRAJIN_ERR;
  Frame tmp( natoms_ );
  if (readFrame(0, tmp)) return TRAJIN_ERR;
  closeTraj();
  // No velocity, no temperature, yes time, no force
  SetCoordInfo( CoordinateInfo(ReplicaDimArray(), tmp.BoxCrd(), false, false, true, false) );
  return TRAJIN_UNK;
}

// Traj_GmxXtc::setupTrajout()
int Traj_GmxXtc::setupTrajout(FileName const& fnameIn, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  if (fnameIn.empty()) return 1;
  fname_ = fnameIn;
  if (!append) {
    SetCoordInfo( cInfoIn );
    natoms_ = trajParm->Natom();
    // Allocate array for writing coords
    if (vec_ != 0) delete[] vec_;
    vec_ = new rvec[ natoms_ ];
    if (vec_ == 0) return 1;
    // Open write
    xd_ = xdrfile_open(fname_.full(), "w");
    if (xd_ == 0) {
      mprinterr("Error: Could not open XTC file for write.\n");
      return 1;
    }
  } else {
    int nframes = setupTrajin( fname_, trajParm );
    if ( nframes == TRAJIN_ERR ) return 1;
    if ( nframes != TRAJIN_UNK )
      mprintf("\tAppending to XTC file starting at frame %i\n", nframes);
    // Re-open for append
    xd_ = xdrfile_open(fname_.full(), "a");
    if (xd_ == 0) {
      mprinterr("Error: Could not open XTC file for append.\n");
      return 1;
    }
  }
  return 0;
}

// Traj_GmxXtc::WriteHelp()
void Traj_GmxXtc::WriteHelp() {
  mprintf("\tdt : Time step to multiply set #s by (default 1.0). Ignored if time already present.\n");
}

// Traj_GmxXtc::processWriteArgs()
int Traj_GmxXtc::processWriteArgs(ArgList& argIn) {
  dt_ = argIn.getKeyDouble( "dt", 1.0 );
  return 0;
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
  if (result != exdrOK) return 1;
  //mprintf("DEBUG: set %i step %i time %f\n", set, step, time);
  frameIn.SetTime( time );
  int idx = 0;
  for (int ix = 0; ix < natoms_; ix++)
    for (int kx = 0; kx < DIM; kx++)
      frameIn[idx++] = (double)vec_[ix][kx] * Constants::NM_TO_ANG;
  idx = 0;
  Matrix_3x3 ucell;
  for (int ii = 0; ii < DIM; ii++)
    for (int ij = 0; ij < DIM; ij++)
      ucell[idx++] = (double)box_[ii][ij] * Constants::NM_TO_ANG;
  frameIn.SetBox( Box(ucell) );
  return 0;
}

// Traj_GmxXtc::writeFrame()
int Traj_GmxXtc::writeFrame(int set, Frame const& frameOut) {
  float time;
  if (CoordInfo().HasTime())
    time = (float)frameOut.Time();
  else
    time = (float)(dt_ * (double)set);
  Matrix_3x3 Ucell = frameOut.BoxCrd().UnitCell( Constants::ANG_TO_NM );
  int idx = 0;
  for (int ii = 0; ii < DIM; ii++)
    for (int ij = 0; ij < DIM; ij++)
      box_[ii][ij] = (float)Ucell[idx++];
  idx = 0;
  for (int ix = 0; ix < natoms_; ix++)
    for (int kx = 0; kx < DIM; kx++)
      vec_[ix][kx] = (float)frameOut[idx++] * Constants::ANG_TO_NM;
  int result = write_xtc(xd_, natoms_, set, time, box_, vec_, prec_);
  if (result != 0) return 1;
  return 0;
}

void Traj_GmxXtc::Info() {
  mprintf("is a GROMACS XTC file,");
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
