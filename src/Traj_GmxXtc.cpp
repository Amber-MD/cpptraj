#include "Traj_GmxXtc.h"
#include "CpptrajStdio.h"

#ifndef NO_XDRFILE
/// CONSTRUCTOR
Traj_GmxXtc::Traj_GmxXtc() : xd_(0), vec_(0), natoms_(0) {}

/// DESTRUCTOR
Traj_GmxXtc::~Traj_GmxXtc() {
  if (xd_ != 0) xdrfile_close(xd_);
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

int Traj_GmxXtc::setupTrajin(FileName const& fname, Topology* trajParm)
{
  return 1;
}

int Traj_GmxXtc::setupTrajout(FileName const& fname, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  return 1;
}

int Traj_GmxXtc::openTrajin() {
  return 1;
}

void Traj_GmxXtc::closeTraj() {

}

int Traj_GmxXtc::readFrame(int set, Frame& frameIn) {
  return 1;
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
