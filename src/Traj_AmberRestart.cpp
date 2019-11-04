// Traj_AmberRestart
#include <cstdio> // sscanf
#include "Traj_AmberRestart.h"
#include "Topology.h"
#include "ArgList.h"
#include "Frame.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NoTrailingWhitespace

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

// FIXME This needs to be changed if ever used for trajout
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

// Traj_AmberRestart::getBoxAngles()
/** Based on input buffer, determine num box coords and get box angles.
  */
int Traj_AmberRestart::getBoxAngles(std::string const& boxline, Box& trajBox) {
  double box[6];
  if (boxline.empty()) {
    mprinterr("Internal Error: Restart box line is empty.\n");
    return 1;
  }
  numBoxCoords_ = sscanf(boxline.c_str(), "%12lf%12lf%12lf%12lf%12lf%12lf",
                         box, box+1, box+2, box+3, box+4, box+5);
  if (debug_>0) {
    mprintf("DEBUG: Restart BoxLine [%s]\n",boxline.c_str());
    mprintf("       Restart numBoxCoords_=%i\n",numBoxCoords_);
  }
  if (numBoxCoords_==-1) {
    // This can occur if there is an extra newline or whitespace at the end
    // of the restart. Warn the user.
    mprintf("Warning: Restart appears to have an extra newline or whitespace.\n");
    mprintf("         Assuming no box information present.\n");
    trajBox.SetNoBox();
    numBoxCoords_ = 0;
  } else if (numBoxCoords_==6) {
    trajBox.SetBox(box);
  } else {
    mprinterr("Error: Expected 6 box coords in restart box coord line, got %i.\n",
              numBoxCoords_);
    return 1;
  }
  return 0;
}

void Traj_AmberRestart::ReadHelp() {
  mprintf("\tusevelascoords: Use velocities in place of coordinates.\n");
}

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
  BufferedFrame infile;
  if (infile.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (infile.OpenFile()) return TRAJIN_ERR;
  readAccess_ = true;
  // Read in title
  std::string title = infile.GetLine();
  SetTitle( NoTrailingWhitespace(title) );
  // Read in natoms, time, and Replica Temp if present
  std::string nextLine = infile.GetLine();
  if (nextLine.empty()) {
    mprinterr("Error: Could not read restart atoms/time.\n");
    return TRAJIN_ERR;
  }
  int restartAtoms = 0;
  bool hasTemp = false;
  bool hasTime = false;
  int nread = sscanf(nextLine.c_str(),"%i %lE %lE",&restartAtoms,&restartTime_,&restartTemp_);
  if (nread < 1) {
    mprinterr("Error: Unable to read restart atoms/time.\n");
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
              infile.Filename().base(), restartAtoms);
    mprinterr("       match number in associated parmtop (%i)\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  natom3_ = restartAtoms * 3;
  // Calculate the length of coordinate frame in bytes
  infile.SetupFrameBuffer( natom3_, 12, 6 );
  // Read past restart coords
  int bytesRead = infile.AttemptReadFrame();
  //mprintf("DEBUG: %i bytes read, frame size is %zu\n", bytesRead, infile.FrameSize());
  if (bytesRead != (int)infile.FrameSize()) {
    // See if we are just missing EOL
    if (bytesRead + 1 + (int)infile.IsDos() == (int)infile.FrameSize())
      mprintf("Warning: File '%s' missing EOL.\n", infile.Filename().full());
    else {
      mprinterr("Error: Error reading coordinates from Amber restart '%s'.\n",
                 infile.Filename().full());
      return TRAJIN_ERR;
    }
  }
  // Save coordinates
  CRD_.resize( natom3_ );
  infile.BufferBegin();
  infile.BufferToDouble(&CRD_[0], natom3_);
  // Attempt a second read to get velocities or box coords
  bool hasVel = false;
  boxInfo_.SetNoBox();
  nread = infile.AttemptReadFrame();
  if ( nread < 0 ) {
    mprinterr("Error: Error attempting to read box line of Amber restart file.\n");
    return TRAJIN_ERR;
  }
  size_t readSize = (size_t)nread;
  //mprintf("DEBUG: Restart readSize on second read = %i\n",readSize);
  // If 0 no box or velo 
  if (readSize > 0) {
    bool velocitiesRead = (readSize == infile.FrameSize());
    if (readSize + 1 + (unsigned int)infile.IsDos() == infile.FrameSize()) {
      mprintf("Warning: File '%s' missing EOL.\n", infile.Filename().full());
      velocitiesRead = true;
    }
    if (velocitiesRead) {
      // If filled framebuffer again, has velocity info. 
      hasVel = true;
      VEL_.resize( natom3_ );
      infile.BufferBegin();
      infile.BufferToDouble(&VEL_[0], natom3_);
      // If we can read 1 more line after velocity, should be box info.
      nextLine = infile.GetLine();
      if (!nextLine.empty()) {
        if (getBoxAngles(nextLine, boxInfo_)) return TRAJIN_ERR;
      } 
    } else if (readSize<82) {
      // If we read something but didnt fill framebuffer, should have box coords.
      nextLine.assign(infile.Buffer(), readSize);
      if (getBoxAngles(nextLine, boxInfo_)) return TRAJIN_ERR;
    } else {
      // Otherwise, who knows what was read?
      mprinterr("Error: AmberRestart::setupTrajin(): When attempting to read in\n"
                "Error: box coords/velocity info got %lu chars, expected 0, 37,\n"
                "Error: 73, or %lu.\n", readSize, infile.FrameSize());
      mprinterr("Error: This usually indicates a malformed or corrupted restart file.\n");
      return TRAJIN_ERR;
    }
  }
  if (useVelAsCoords_ && !hasVel) {
    mprinterr("Error: 'usevelascoords' specified but no velocities in this restart.\n");
    return TRAJIN_ERR;
  }
  infile.CloseFile();
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
    std::copy(boxInfo_.boxPtr(), boxInfo_.boxPtr()+6, frameIn.bAddress());
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
  if (numBoxCoords_!=0)
    file_.DoubleToBuffer(frameOut.bAddress(), numBoxCoords_, "%12.7f");

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
