#include "Traj_GmxXtc.h"
#include "CpptrajStdio.h"

#ifndef NO_XDRFILE
/// CONSTRUCTOR
Traj_GmxXtc::Traj_GmxXtc() : xd_(0), vec_(0) {}

/// DESTRUCTOR
Traj_GmxXtc::~Traj_GmxXtc() {
  if (xd_ != 0) xdrfile_close(xd_);
  if (vec_ != 0) delete[] vec_;
}

#else
// =============================================================================
Traj_GmxXtc::Traj_GmxXtc() {}

Traj_GmxXtc::~Traj_GmxXtc() {}
#endif
