#include "DataIO_Coords.h"
#include "BufferedLine.h"
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
/** Special case: read in coordinates as CSV format. Assumes no topology,
  * creates a fake one. Intended to test against MDANCE.
  */
int DataIO_Coords::readAsCSV(DataSet* dsetIn, DataSet::DataType setType,
                             FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  static const char* SEP = ",\r";
  DataSet* dset = dsetIn;
  // First open the file and make sure there are commas
  BufferedLine infile;
  if (infile.OpenFileRead(fname)) {
    return 1;
  }
  std::string firstLine = infile.GetLine();
  if (firstLine.empty()) {
    mprinterr("Error: No lines in CSV file '%s'\n", fname.full());
    return 1;
  }
  ArgList line(firstLine, SEP);
  // Number of arguments is number of coords
  int ncoords = line.Nargs();
  if (ncoords < 3) {
    mprinterr("Error: Less than 3 coordinates in CSV file (%i)\n", ncoords);
    return 1;
  }
  if ( (ncoords%3) != 0 ) {
    mprinterr("Error: Number of coords (%i) is not a multiple of 3\n", ncoords);
    return 1;
  }
  int natoms = ncoords / 3;
  mprintf("\t%i atoms, %i coords.\n", natoms, ncoords);

  Topology top;
  Topology* topPtr = 0;
  if (dset == 0) {
    // Develop the pseudo-topology
    for (int iat = 0; iat != natoms; iat++)
      top.addTopAtom( Atom("C", "C"),
                      Residue("MOL", iat, ' ', ""), iat, false );
    top.CommonSetup(false, false);
    top.Summary();
    topPtr = &top;
  } else {
    topPtr = ((DataSet_Coords*)dset)->TopPtr();
    if (topPtr->Natom() != natoms) {
      mprinterr("Error: Atom mismatch between CSV (%i) and '%s' (%i)\n",
                natoms, dset->legend(), topPtr->Natom());
      return 1;
    }
  }

  // If no data set yet, set it up
  if (dset == 0) {
    MetaData md( fname, dsname, -1 );
    dset = dsl.AddSet(setType, md);
    if (dset == 0) return 1;
    DataSet_Coords* coords = static_cast<DataSet_Coords*>( dset );
    // Blank CoordinateInfo(), only COORDS
    if (coords->CoordsSetup( *topPtr, CoordinateInfo() )) {
      mprinterr("Error: Could not set up COORDS set %s\n", coords->legend());
      return 1;
    }
  }

  // Read coords
  int ifrm = 0;
  Frame frameIn( natoms );
  while (!firstLine.empty()) {
    if (ifrm > 0) {
      line.SetList( firstLine, SEP );
      if (line.Nargs() != ncoords) {
        mprinterr("Error: # of coordinates changes from %i to %i at line %i\n",
                  ncoords, line.Nargs(), infile.LineNumber());
        break;
      }
    }
    frameIn.ClearAtoms();
    //int icrd = 0;
    for (int iat = 0; iat != natoms; iat++) {
      double XYZ[3];
      XYZ[0] = line.getNextDouble(0);
      XYZ[1] = line.getNextDouble(0);
      XYZ[2] = line.getNextDouble(0);
      frameIn.AddXYZ( XYZ );
      //icrd += 3
    }
    //line.PrintDebug();
    // Sanity check
    if (line.CheckForMoreArgs()) {
      mprinterr("Error: Not enough double values read for line %i\n", infile.LineNumber());
      break;
    }
    ((DataSet_Coords*)dset)->AddFrame( frameIn );
    ifrm++;
    firstLine = infile.GetLine();
  }
  mprintf("\tRead in %i frames.\n", ifrm);
  AddedByMe( dset );

  return 0;
}

// DataIO_Coords::ReadData()
int DataIO_Coords::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  ClearAddedByMe();
  bool read_as_csv = false;
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
      // Special cases.
      // Check for .csv extension
      if (fname.Ext() == ".csv") {
        mprintf("\tAssuming coordinates stored in CSV format.\n");
        read_as_csv = true;
      } else {
        mprinterr("Error: '%s' does not have parm/coords info.\n", fname.full());
        return 1;
      }
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

  if (read_as_csv) {
    // Special case: read as CSV file
    return readAsCSV(dset, setType, fname, dsl, dsname);
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
