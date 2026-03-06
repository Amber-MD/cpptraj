#include "DataIO_Coords.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"
#include "TrajectoryFile.h"
#include "Trajin_Single.h"

/// CONSTRUCTOR
DataIO_Coords::DataIO_Coords() //:
  //is_parm_fmt_(false),
  //is_traj_fmt_(false)
{

}

// DataIO_Coords::ID_DataFormat()
/** NOTE: This is disabled intentionally. The problem is that some
  *       data formats looks very similar to trajectory formats,
  *       so that e.g. a cpptraj vector data file can look like
  *       an Amber ASCII trajectory.
  *       Users can use commands like 'loadcrd', or force the read
  *       with the 'as' keyword.
  */
bool DataIO_Coords::ID_DataFormat(CpptrajFile& infile)
{
  return false;
}

// DataIO_Coords::ReadHelp()
void DataIO_Coords::ReadHelp()
{

}

// DataIO_Coords::processReadArgs()
int DataIO_Coords::processReadArgs(ArgList& argIn)
{

  return 0;
}

/// \return True if this is a COORDS set we can append to
static inline bool can_append(DataSet::DataType typeIn) {
  return (typeIn == DataSet::COORDS ||
          typeIn == DataSet::FRAMES);
}

// DataIO_Coords::ReadData()
int DataIO_Coords::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  ClearAddedByMe();
  DataSet::DataType setType = DataSet::COORDS; // FIXME make user option
  //if (!is_parm_fmt_ && !is_traj_fmt_) {
    bool is_parm_fmt_ = false;
    bool is_traj_fmt_ = false;
    // Assume that ID_DataFormat() has not been called.
    //CpptrajFile tmpfile;
    //if (tmpfile.SetupWrite( fname, debug_ )) {
    //  mprinterr("Error: Could not setup check for parm/coords info in '%s'.\n",
    //            fname.full());
    //  return 1;
    //}
    // Needs to be either a topology format or a coords format
    ParmFile::ParmFormatType parm_format = ParmFile::DetectFormat( fname );
    TrajectoryFile::TrajFormatType traj_format = TrajectoryFile::DetectFormat( fname );
    is_parm_fmt_ = (parm_format != ParmFile::UNKNOWN_PARM);
    is_traj_fmt_ = (traj_format != TrajectoryFile::UNKNOWN_TRAJ);
    if (!is_parm_fmt_ && !is_traj_fmt_) {
      mprinterr("Error: '%s' does not have parm/coords info.\n", fname.full());
      return 1;
    }
  //}

  DataSet* dset = 0;
  if (!dsname.empty()) {
    // Is this set already present?
    DataSetList selectedDS = dsl.SelectSets( dsname );
    if (!selectedDS.empty()) {
      dset = selectedDS[0];
      if (selectedDS.size() > 1)
        mprintf("Warning: %s selects multiple data sets, only using the first (%s)\n", dsname.c_str(), dset->legend());
    }
    if (dset != 0) {
      if (!can_append(dset->Type())) {
        mprinterr("Error: Cannot append coordinates to existing set '%s'\n", dset->legend());
        return 1;
      } else
        mprintf("\tAppending to set '%s'\n", dset->legend());
    }
  }

  // Topology read/setup
  Topology topIn;
  Topology* topPtr = 0;
  if (dset == 0) {
    if (!is_parm_fmt_) {
      mprinterr("Error: '%s' does not contain any topology information.\n", fname.full());
      return 1;
    }
    // No data set yet; read topology info
    ParmFile pfile;
    ArgList topargs;
    if (pfile.ReadTopology( topIn, fname, topargs, debug_ )) {
      mprinterr("Error: Could not read topology information from '%s'\n", fname.full());
      return 1;
    }
    topPtr = &topIn;
  } else {
    topPtr = ((DataSet_Coords*)dset)->TopPtr();
  }

  // Trajectory setup
  Trajin_Single trajin;
  if (is_traj_fmt_) {
    trajin.SetDebug( debug_ );
    ArgList trajargs;
    if (trajin.SetupTrajRead( fname, trajargs, topPtr )) {
      mprinterr("Error: Could not set up trajectory info for '%s'\n", fname.full());
      return 1;
    }
  } 

  // If no data set yet, set it up
  if (dset == 0) {
    MetaData md( fname, dsname, -1 );
    dset = dsl.AddSet(setType, md);
    if (dset == 0) return 1;
    DataSet_Coords* coords = static_cast<DataSet_Coords*>( dset );
    if (coords->CoordsSetup( *topPtr, trajin.TrajCoordInfo() )) { // FIXME is this ok for no traj info?
      mprinterr("Error: Could not set up COORDS set %s\n", coords->legend());
      return 1;
    }
  }

  // Trajectory read
  if (is_traj_fmt_) {
    Frame frameIn;
    frameIn.SetupFrameV(topPtr->Atoms(), trajin.TrajCoordInfo());
    trajin.BeginTraj();
    trajin.Traj().PrintInfoLine();
    while (trajin.GetNextFrame( frameIn ))
      ((DataSet_Coords*)dset)->AddFrame( frameIn );
    trajin.EndTraj();
  }
  AddedByMe( dset );

  return 0;
}

// DataIO_Coords::WriteHelp()
void DataIO_Coords::WriteHelp()
{

}

// DataIO_Coords::processWriteArgs()
int DataIO_Coords::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Coords::WriteData()
int DataIO_Coords::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
