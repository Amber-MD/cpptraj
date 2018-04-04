#include "DataSet_Coords_Trajin.h"
#include "TrajectoryFile.h"
#include "CpptrajStdio.h"

DataSet_Coords_Trajin::DataSet_Coords_Trajin() :
  Traj_(0),
  nframes_(0)
{}

DataSet_Coords_Trajin::~DataSet_Coords_Trajin() {
  if (Traj_ != 0) delete Traj_;
}

// DataSet_Coords_Trajin::AddSingleTrajin()
int DataSet_Coords_Trajin::AddSingleTrajin(std::string const& fnameIn, ArgList& argIn,
                                           Topology* topIn)
{
  if (topIn==0) {
    mprinterr("Internal Error: DataSet_Coords_Trajin setup called with null topology.\n");
    return 1;
  }
  if (Traj_ != 0) delete Traj_;
  FileName fname( fnameIn );
  TrajectoryFile::TrajFormatType tformat;
  if ( (Traj_ = TrajectoryFile::DetectFormat( fname, tformat )) == 0 ) {
    mprinterr("Error: Could not determine trajectory %s format.\n", fname.full());
    return 1;
  }
  mprintf("\tReading '%s' as %s\n", fname.base(), TrajectoryFile::FormatString(tformat));

  if (Traj_->processReadArgs( argIn )) return 1;
  // Set up the format for reading and get the number of frames.
  nframes_ = Traj_->setupTrajin(fname, topIn);
  if (nframes_ == TrajectoryIO::TRAJIN_ERR) {
    mprinterr("Error: Could not set up '%s' for reading.\n", fname.full());
    return 1;
  }
  if (CoordsSetup(*topIn, Traj_->CoordInfo())) return 1;
  readFrame_.SetupFrameV( topIn->Atoms(), Traj_->CoordInfo() );

  return 0;
}

void DataSet_Coords_Trajin::Info() const { return ; }

// DataSet_Coords_Trajin::CoordsSetup()
int DataSet_Coords_Trajin::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn) {
  top_ = topIn;
  cInfo_ = cInfoIn;
  return 0;
}


