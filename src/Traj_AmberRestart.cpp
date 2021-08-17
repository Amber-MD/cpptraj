// Traj_AmberRestart
#include <cstdio> // sscanf
#include <algorithm> // std::copy
#include "Traj_AmberRestart.h"
#include "Topology.h"
#include "ArgList.h"
#include "Frame.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NoTrailingWhitespace
#include "TextBlockBuffer.h"

// CONSTRUCTOR
Traj_AmberRestart::Traj_AmberRestart() :
  natom3_(0),
  numBoxCoords_(0),
  restartTime_(-1.0),
  restartTemp_(-1.0),
  time0_(1.0),
  dt_(1.0),
  singleWrite_( false),
  readAccess_(false),
  useVelAsCoords_(false),
  outputTemp_(false),
  outputTime_(true), // For backwards compat.
  prependExt_(false)
{}

/** Check for an integer (I5) followed by 0-2 scientific floats (E15.7) */
bool Traj_AmberRestart::ID_TrajFormat(CpptrajFile& fileIn) {
  // Assume file set up for read
  if (fileIn.OpenFile()) return false;
  bool isRestart = false;
  if ( fileIn.NextLine() !=0 ) { // Title
    const char* ptr = fileIn.NextLine(); // Natom [time [temp]]
    if (ptr != 0) {
      int i0;
      double D[3];
      int nread = sscanf(ptr, "%5i%15lf%15lf%lf", &i0, D, D+1, D+2);
      if (nread > 0 && nread < 4) {
        // Read at least 3 12.7 coordinates from next line.
        ptr = fileIn.NextLine();
        if (ptr != 0) {
          nread = sscanf(ptr, "%12lf%12lf%12lf", D, D+1, D+2);
          if (nread == 3)
            isRestart = true;
        }
      }
    }
  }
  fileIn.CloseFile();
  return isRestart;
}

// Traj_AmberRestart::closeTraj()
void Traj_AmberRestart::closeTraj() {
  file_.CloseFile();
}

// TODO This needs to be changed if ever used for trajout
// Traj_AmberRestart::openTrajin()
int Traj_AmberRestart::openTrajin() { return 0; }

void Traj_AmberRestart::WriteHelp() {
  mprintf("\tremdtraj: Write temperature to restart file (will also write time).\n"
          "\ttime0   : Time for first frame (if not specified time is not written).\n"
          "\tdt      : Time step for subsequent frames, t=(time0+frame)*dt; (default 1.0)\n"
          "\tkeepext : Keep filename extension; write '<name>.<num>.<ext>' instead.\n");
}

// Traj_AmberRestart::processWriteArgs()
int Traj_AmberRestart::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {
  outputTemp_ = argIn.hasKey("remdtraj");
  time0_ = argIn.getKeyDouble("time0", -1.0);
  dt_ = argIn.getKeyDouble("dt",1.0);
  singleWrite_ = argIn.hasKey("single");
  prependExt_ = argIn.hasKey("keepext");
  return 0;
}

// Traj_AmberRestart::setupTrajout()
/** Allocate a character buffer based on number of coords and whether 
  * velocities/box info is present.
  */
int Traj_AmberRestart::setupTrajout(FileName const& fname, Topology* trajParm,
                                    CoordinateInfo const& cInfoIn, 
                                    int NframesToWrite, bool append)
{
  if (append) {
    mprinterr("Error: Append not supported for Amber Restart.\n");
    return 1;
  }
  CoordinateInfo cInfo = cInfoIn;
  if (!cInfo.HasTemp() && outputTemp_) cInfo.SetTemperature(true);
  // If temperature requested write time as well or format will break.
  if (cInfo.HasTemp()) {
    outputTime_ = true;
    if (!cInfo.HasTime() && time0_ < 0.0) time0_ = 1.0;
  }
  if (outputTime_) {
    if (!cInfo.HasTime() && time0_ >= 0) cInfo.SetTime(true);
  } else
    cInfo.SetTime(false);
  SetCoordInfo( cInfo );
  if (file_.SetupWrite( fname, debug_ )) return 1;
  readAccess_ = false;
  // Set trajectory info
  natom3_ = trajParm->Natom() * 3;
  // Calculate the length of coordinate frame in bytes
  file_.SetupFrameBuffer( natom3_, 12, 6 ); 
  // Dont know ahead of time if velocities will be used, allocate space
  // just in case. Velocity will not be written if V input is null.
  file_.ResizeBuffer( natom3_ );
  // If box coords are present, allocate extra space for them
  if (CoordInfo().HasBox()) {
    if (!CoordInfo().TrajBox().Is_X_Aligned())
      mprintf("Warning: Unit cell is not X-aligned. Box cannot be properly stored as Amber ASCII restart.\n");
    numBoxCoords_ = 6;
    file_.ResizeBuffer( numBoxCoords_ );
  }
  // If number of frames to write == 1 set singleWrite so we dont append
  // frame # to filename.
  if (NframesToWrite==1) singleWrite_ = true;
  // Set up title
  std::string outTitle = Title();
  if (outTitle.empty()) {
    outTitle.assign("Cpptraj Generated Restart");
    outTitle.resize(80, ' ');
  } else {
    if ( outTitle.size() > 80) {
      mprintf("Warning: Amber restart title for %s too long: truncating.\n[%s]\n",
              file_.Filename().base(), outTitle.c_str());
      outTitle.resize(80);
    }
  }
  SetTitle( outTitle );

  return 0;
}

// Traj_AmberRestart::ReadHelp()
void Traj_AmberRestart::ReadHelp() {
  mprintf("\tusevelascoords: Use velocities in place of coordinates.\n");
}

// Traj_AmberRestart::processReadArgs(ArgList& argIn)
int Traj_AmberRestart::processReadArgs(ArgList& argIn) {
  useVelAsCoords_ = argIn.hasKey("usevelascoords");
  return 0;
}

// Traj_AmberRestart::setupTrajin()
/** Set up and read Amber restart file. Coordinate/velocities will be saved
  * here to avoid having to open the file again. Check that number of atoms 
  * matches number of atoms in associated parmtop. Check for box/velocity info.
  */
int Traj_AmberRestart::setupTrajin(FileName const& fname, Topology* trajParm)
{
  TextBlockBuffer tfile;
  tfile.SetDebug( debug_ );
  // In addition to the space needed to buffer the block of coordinates,
  // set up for 3 additional lines of 80 bytes: title, natoms etc, and
  // potential velocities/box.
  if (tfile.OpenFileRead( fname, trajParm->Natom()*3, 12, 6, 240 )) return TRAJIN_ERR;
  readAccess_ = true;

  // Read in title
  const char* ptr = tfile.Line();
  if (ptr == 0) {
    mprinterr("Error: Null encountered when reading Amber restart title.\n");
    return TRAJIN_ERR;
  }
  std::string title(ptr);
  SetTitle( NoTrailingWhitespace(title) );

  // Read in natoms, time, and Replica Temp if present
  ptr = tfile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read restart atoms/time.\n");
    return TRAJIN_ERR;
  }
  int restartAtoms = 0;
  bool hasTemp = false;
  bool hasTime = false;
  int nread = sscanf(ptr, "%i %lE %lE", &restartAtoms, &restartTime_, &restartTemp_);
  if (nread < 1) {
    mprinterr("Error: Unable to get restart atoms/time from line:\n"
              "Error:   [%s]\n", ptr);
    return TRAJIN_ERR;
  } else if (nread == 1) { // # atoms only
    restartTime_ = 0.0;
    restartTemp_ = -1.0;
  } else if (nread == 2) { // # atoms and time
    hasTime = true;
    restartTemp_ = -1.0;
  } else {                 // # atoms, time, and temperature
    hasTime = true;
    hasTemp = true; 
  }
  if (debug_ > 0) 
    mprintf("\tAmber restart: Atoms=%i Time=%lf Temp=%lf\n",restartAtoms,
            restartTime_, restartTemp_);

  // Check that natoms matches parm natoms
  if (restartAtoms != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in Amber Restart %s (%i) does not\n",
              tfile.Filename().base(), restartAtoms);
    mprinterr("       match number in associated parmtop (%i)\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  natom3_ = restartAtoms * 3;

  // Attemp to read coordinates from restart
  CRD_.resize( natom3_ );
  int eltsRead = tfile.BlockToDoubles( &CRD_[0] );
  if (debug_ > 0)
    mprintf("DEBUG: CRD: %i elts read\n", eltsRead);
  if (eltsRead != natom3_) {
    mprinterr("Error: Could not read coordinates from Amber restart.\n");
    return TRAJIN_ERR;
  }

  // Attempt a second read to get velocities or box coords.
  // In order for this to be a valid read we need at least 36 characters
  // (3 cols * 12 chars, minimum 1 atom case) remaining.
  //mprintf("DEBUG: Characters remaining in buffer= %li\n", tfile.Nremaining());
  //mprintf("DEBUG: [");
  //const char* c = tfile.CurrentLine();
  //for (long int ic = 0; ic < tfile.Nremaining(); ic++)
  //  mprintf("%c", *(c+ic));
  //mprintf("]\n");
  bool hasVel = false;
  boxInfo_.SetNoBox();
  if (tfile.Nremaining() > 35) {
    VEL_.resize( natom3_ );
    eltsRead = tfile.BlockToDoubles( &VEL_[0] );
    if (debug_ > 0)
      mprintf("DEBUG: VEL: %i elts read\n", eltsRead);
    if (eltsRead == natom3_) {
      // Velocity info present
      hasVel = true;
      // Check for box
      tfile.SetupTextBlock(6, 12, 6);
      double xyzabg[6];
      eltsRead = tfile.BlockToDoubles( xyzabg );
      if (eltsRead == 6) {
        boxInfo_.SetupFromXyzAbg(xyzabg);
      } else if (eltsRead != 0) {
        mprinterr("Error: When checking if box info present, expected 6 elements,\n"
                  "Error:  got %i\n", eltsRead);
        return TRAJIN_ERR;
      }
    } else if (eltsRead == 6) {
      // Only box info present
      double xyzabg[6];
      for (unsigned int ii = 0; ii != 6; ii++)
        xyzabg[ii] = VEL_[ii];
      VEL_.clear();
      boxInfo_.SetupFromXyzAbg(xyzabg);
    } else if (eltsRead != 0) {
      mprinterr("Error: When checking if velocity/box info present, expected either\n"
                "Error:  %i or 6 elements, got %i\n", natom3_, eltsRead);
      return TRAJIN_ERR;
    }
  }
  if (boxInfo_.HasBox())
    numBoxCoords_ = 6;

  tfile.CloseFile();
  // Set coordinate info
  SetCoordInfo( CoordinateInfo(boxInfo_, hasVel, hasTemp, hasTime) );
  // Only 1 frame in restart by definition
  return 1;
}

// Traj_AmberRestart::readFrame()
/** Copy buffered coords/velocities/box to input frame. */
int Traj_AmberRestart::readFrame(int set, Frame& frameIn) {
  // Set frame temp
  if (CoordInfo().HasTemp())
    frameIn.SetTemperature( restartTemp_ );
  // Set frame time
  if (CoordInfo().HasTime())
    frameIn.SetTime( restartTime_ );
  // Get coords from buffer
  std::copy(CRD_.begin(), CRD_.end(), frameIn.xAddress());
  // Get velocity from buffer if present
  if (CoordInfo().HasVel()) {
    if (frameIn.HasVelocity()) {
      if (useVelAsCoords_)
        std::copy(VEL_.begin(), VEL_.end(), frameIn.xAddress());
      else
        std::copy(VEL_.begin(), VEL_.end(), frameIn.vAddress());
    }
  }
  // Get box from buffer if present
  if (numBoxCoords_!=0)
    frameIn.SetBox( boxInfo_ ); 
  return 0;
}

// Traj_AmberRestart::readVelocity()
int Traj_AmberRestart::readVelocity(int set, Frame& frameIn) {
  if (CoordInfo().HasVel()) {
    std::copy(VEL_.begin(), VEL_.end(), frameIn.vAddress());
    return 0;
  }
  return 1;
}

// Traj_AmberRestart::writeFrame()
/** Write coords in Frame to file in amber restart format. */
int Traj_AmberRestart::writeFrame(int set, Frame const& frameOut) {
  // If just writing 1 frame dont modify output filename
  if (singleWrite_) {
    if (file_.OpenFile()) return 1;
  } else {
    if (file_.OpenWriteNumbered( set + 1, prependExt_ ) ) return 1;
  }

  // Write out title
  file_.Printf("%-s\n", Title().c_str());
  // Write out atoms
  file_.Printf("%5i", frameOut.Natom());
  // Write out restart time
  if (CoordInfo().HasTime()) {
    if (time0_>=0)
      restartTime_ = (time0_ + (double)set) * dt_;
    else
      restartTime_ = frameOut.Time();
    file_.Printf("%15.7lE",restartTime_);
  }
  // Write out temperature
  if (CoordInfo().HasTemp())
    file_.Printf("%15.7lE",frameOut.Temperature());
  file_.Printf("\n");

  // Write coords to buffer
  file_.BufferBegin();
  file_.DoubleToBuffer(frameOut.xAddress(), natom3_, "%12.7f");
  // Write velocity to buffer. Check V since velocity not known ahead of time
  if (CoordInfo().HasVel() && frameOut.HasVelocity())
    file_.DoubleToBuffer(frameOut.vAddress(), natom3_, "%12.7f");
  // Write box to buffer
  if (numBoxCoords_!=0) {
    if (!frameOut.BoxCrd().Is_X_Aligned())
      mprintf("Warning: Set %i; unit cell is not X-aligned. Box cannot be properly stored as Amber ASCII restart.\n", set+1);
    file_.DoubleToBuffer(frameOut.BoxCrd().XyzPtr(), numBoxCoords_, "%12.7f");
  }

  if (file_.WriteFrame()) return 1;

  file_.CloseFile();

  return 0;
}

// Traj_AmberRestart::Info()
void Traj_AmberRestart::Info() {
  mprintf("is an AMBER restart file");
  if (readAccess_) {
    // If read access we know for sure whether there are velocities.
    if (CoordInfo().HasVel())
      mprintf(" with velocity info");
    else
      mprintf(", no velocities");
    if (useVelAsCoords_) mprintf(" (using velocities as coords)");
  }
}
#ifdef MPI
/// Since files are opened on write this does not need to do anything
int Traj_AmberRestart::parallelOpenTrajout(Parallel::Comm const& commIn) { return 0; }

/** No file access during setupTrajout, so have all threads call it.
  * No need to sync.
  */
int Traj_AmberRestart::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{
  return (setupTrajout(fname, trajParm, cInfoIn, NframesToWrite, append));
}

int Traj_AmberRestart::parallelWriteFrame(int set, Frame const& frameOut) {
  return ( writeFrame(set, frameOut) );
}
#endif
