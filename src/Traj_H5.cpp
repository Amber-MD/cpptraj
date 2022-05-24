#include "Traj_H5.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#ifdef BINTRAJ
# include <netcdf.h>
# include "NC_Routines.h"
#endif

/// CONSTRUCTOR
Traj_H5::Traj_H5() {}

#ifdef HAS_HDF5
bool Traj_H5::HasConventions(int ncid) {
  std::string attrText = NC::GetAttrText(ncid, "conventions");
  if (attrText.empty())
    return false;
  if (attrText != "Pande") {
    mprinterr("Error: H5 file has conventions string that is not 'Pande'.\n");
    return false;
  }
  attrText = NC::GetAttrText(ncid, "conventionVersion");
  if ( attrText != "1.1")
    mprintf("Warning: H5 file has conventionVersion that is not 1.1 (%s)\n", attrText.c_str());
  return true;
}
#endif

/** Identify trajectory format. File should be setup for READ */
bool Traj_H5::ID_TrajFormat(CpptrajFile& fileIn) {
# ifdef HAS_HDF5
  int myNcid;
  if ( nc_open( fileIn.Filename().full(), NC_NOWRITE, &myNcid ) != NC_NOERR )
    return false;
  if (HasConventions(myNcid)) return true;
  nc_close( myNcid );
# else
  unsigned char buf[8];
  unsigned int nread = fileIn.Read(buf, 8);
  if (nread > 7 && buf[0] == 0x89 && buf[1] == 0x48 && buf[2] == 0x44 && buf[3] == 0x46 &&
                   buf[4] == 0x0d && buf[5] == 0x0a && buf[6] == 0x1a && buf[7] == 0x0a)
  {
    mprintf("Warning: File '%s' appears to be HDF5 but cpptraj was compiled without HDF5 support.\n", fname);
  }
# endif
  return false;
}

/** Print trajectory info to stdout. */
void Traj_H5::Info() {
  mprintf("is a <type>");
}

/** Close file. */
void Traj_H5::closeTraj() {

}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_H5::openTrajin() {

  return 0;
}

/** Read help */
void Traj_H5::ReadHelp() {

}

/** Process read arguments. */
int Traj_H5::processReadArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_H5::setupTrajin(FileName const& fname, Topology* trajParm)
{

  return TRAJIN_ERR;
}

/** Read specified trajectory frame. */
int Traj_H5::readFrame(int set, Frame& frameIn) {

  return 0;
}

/** Read velocities from specified frame. */
int Traj_H5::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_H5::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_H5::WriteHelp() {

}

/** Process write arguments. */
int Traj_H5::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_H5::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int Traj_H5::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_H5::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_H5::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_H5::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_H5::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_H5::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_H5::parallelCloseTraj() {

}
#endif
