// Traj_AmberRestart
#include <cstdio> // sscanf
#include <cstring> // strlen: file detection
#include <cctype> // isdigit, isspace: file detection
#include "Traj_AmberRestart.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename

// CONSTRUCTOR
Traj_AmberRestart::Traj_AmberRestart() :
  debug_(0),
  restartAtoms_(0),
  natom3_(0),
  numBoxCoords_(0),
  restartTime_(0.0),
  restartTemp_(-1.0),
  time0_(1.0),
  dt_(1.0),
  singleWrite_( false),
  HasT_(false),
  HasV_(true)
{}

/** Check for an integer (I5) followed by 0-2 scientific floats (E15.7) */
bool Traj_AmberRestart::ID_TrajFormat(CpptrajFile& fileIn) {
  const size_t BUF_SIZE = 128; // NOTE: Should only need 82 for restart
  char buffer2[BUF_SIZE];
  // Assume file set up for read
  if (fileIn.OpenFile()) return false;
  fileIn.Gets(buffer2, BUF_SIZE); // Title
  fileIn.Gets(buffer2, BUF_SIZE); // natom, [time, temp]
  fileIn.CloseFile();
  if (strlen(buffer2)<=36) {
    int i=0;
    for (; i<5; i++) {
      if (!isspace(buffer2[i]) && !isdigit(buffer2[i])) break;
    }
    if ( i==5 ) {
      if (debug_>0) mprintf("  AMBER RESTART file\n");
      return true;
    }
  }
  return false;
}

// Traj_AmberRestart::closeTraj()
void Traj_AmberRestart::closeTraj() {
  file_.CloseFile();
}

// Traj_AmberRestart::openTraj()
/** Open the restart file. Get title, time, restart atoms, temperature
  */
int Traj_AmberRestart::openTraj() {
  int nread; // Dont declare variables inside a switch block
  std::string nextLine;
  switch (file_.Access()) {
    case CpptrajFile::READ :
      if (file_.OpenFile()) return 1;
      // Read in title
      title_ = file_.GetLine();
      // Read in natoms, time, and Replica Temp if present
      nextLine = file_.GetLine();
      if (nextLine.empty()) {
        mprinterr("Error: AmberRestart::open(): Reading restart atoms/time.\n");
        return 1;
      }
      nread = sscanf(nextLine.c_str(),"%i %lE %lE",&restartAtoms_,&restartTime_,&restartTemp_);
      if (nread<1) {
        mprinterr("Error: AmberRestart::open(): Getting restart atoms/time.\n");
        return 1;
      } else if (nread==1) {
        restartTime_=0.0;
        restartTemp_=-1.0;
        HasT_ = false;
      } else if (nread==2) {
        restartTemp_=-1.0;
        HasT_ = false;
      } else {
        HasT_ = true;
      }
      if (debug_>0) 
        mprintf("  Amber restart: Atoms=%i Time=%lf Temp=%lf\n",restartAtoms_,
                restartTime_, restartTemp_);
      break;
    case CpptrajFile::APPEND :
      mprinterr("Error: Append not supported for Amber Restart files.\n");
      return 1;
      break;
    case CpptrajFile::WRITE : // Nothing to do; open/close handled in writeFrame
      break;
  }

  return 0; 
}

// Traj_AmberRestart::SetNoVelocity()
void Traj_AmberRestart::SetNoVelocity() {
  HasV_ = false;
}

// Traj_AmberRestart::processWriteArgs()
int Traj_AmberRestart::processWriteArgs(ArgList& argIn) {
  // For write, assume we want velocities unless specified
  HasV_ = true;
  if (argIn.hasKey("novelocity")) this->SetNoVelocity();
  time0_ = argIn.getKeyDouble("time0", OUTPUTFRAMESHIFT);
  //if (argIn.hasKey("remdtraj")) this->SetTemperature();
  dt_ = argIn.getKeyDouble("dt",1.0);
  return 0;
}

// Traj_AmberRestart::setupTrajout()
/** Allocate a character buffer based on number of coords and whether 
  * velocities/box info is present.
  */
int Traj_AmberRestart::setupTrajout(std::string const& fname, Topology* trajParm,
                     int NframesToWrite, TrajInfo const& tinfo, bool append)
{
  if (append) {
    mprinterr("Error: Append not supported for Amber Restart.\n");
    return 1;
  }
  if (file_.SetupWrite( fname, debug_ )) return 1;
  // Set trajectory info
  title_ = tinfo.Title;
  HasV_ = tinfo.HasV;
  HasT_ = tinfo.HasT;
  restartAtoms_ = trajParm->Natom();
  natom3_ = restartAtoms_ * 3;
  // Calculate the length of coordinate frame in bytes
  // Dont know ahead of time if velocities will be used, allocate space
  // just in case. Velocity will not be written if V input is NULL.
  file_.SetupFrameBuffer( natom3_ * 2, 12, 6, 0 ); 
  // If box coords are present, allocate extra space for them
  if (tinfo.BoxInfo.Type() != Box::NOBOX) {
    numBoxCoords_ = 6;
    file_.ResizeBuffer( numBoxCoords_ );
  }
  // If number of frames to write == 1 set singleWrite so we dont append
  // frame # to filename.
  if (NframesToWrite==1) singleWrite_=true;
  // Set up title
  if (title_.empty()) {
    title_.assign("Cpptraj Generated Restart");
    title_.resize(80, ' ');
  } else {
    if ( title_.size() > 80) {
      title_.resize(80);
      mprintf("Warning: Amber restart title for %s too long: truncating.\n[%s]\n",
              file_.BaseFileStr(), title_.c_str());
    }
  }

  return 0;
}

// Traj_AmberRestart::getBoxAngles()
/** Based on input buffer, determine num box coords and get box angles.
  */
int Traj_AmberRestart::getBoxAngles(const char *boxline, Box& trajBox) {
  double box[6];
  numBoxCoords_ = sscanf(boxline, "%12lf%12lf%12lf%12lf%12lf%12lf",
                        box,box+1,box+2,box+3,box+4,box+5);
  if (debug_>0) {
    mprintf("DEBUG: Restart BoxLine [%s]\n",boxline);
    mprintf("       Restart numBoxCoords_=%i\n",numBoxCoords_);
  }
  if (numBoxCoords_==-1) {
    // This can occur if there is an extra newline or whitespace at the end
    // of the restart. Warn the user.
    mprintf("Warning: Restart [%s] appears to have an extra newline or whitespace.\n",
            file_.BaseFileStr());
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

// Traj_AmberRestart::setupTrajin()
/** Set up amber restart file for reading. Check that number of atoms matches
  * number of atoms in associated parmtop. Check for box/velocity info.
  */
int Traj_AmberRestart::setupTrajin(std::string const& fname, Topology* trajParm,
                    TrajInfo& tinfo)
{
  if (file_.SetupRead( fname, debug_ )) return -1;
  if (openTraj()) return -1; // Gets title, time, natoms, and temp if present

  // Check that natoms matches parm natoms
  if (restartAtoms_ != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in Amber Restart %s (%i) does not\n",
              file_.BaseFileStr(), restartAtoms_);
    mprinterr("       match number in associated parmtop (%i)\n",trajParm->Natom());
    return -1;
  }
  natom3_ = restartAtoms_ * 3;

  // Calculate the length of coordinate frame in bytes
  file_.SetupFrameBuffer( natom3_, 12, 6, 0 );
  coordSize_ = file_.FrameSize();

  // Read past restart coords 
  if ( file_.ReadFrame() == -1 ) {
    mprinterr("Error: AmberRestart::setupTrajin(): Error reading coordinates.\n");
    return -1; 
  }

  // Attempt a second read to get velocities or box coords
  size_t lineSize = (size_t)file_.ReadFrame();
  //mprintf("DEBUG: Restart lineSize on second read = %i\n",lineSize);
  if (lineSize<=0) {
    // If 0 or -1 no box or velo 
    HasV_ = false;
  } else if (lineSize == file_.FrameSize()) {
    // If filled framebuffer again, has velocity info. 
    HasV_ = true;
    // If we can read 1 more line after velocity, should be box info.
    std::string nextLine = file_.GetLine();
    if (!nextLine.empty()) {
      if (getBoxAngles(nextLine.c_str(), tinfo.BoxInfo)) return -1;
    } 
  } else if (lineSize<82) {
    // If we read something but didnt fill framebuffer, should have box coords.
    HasV_ = false;
    if (getBoxAngles(file_.Buffer(), tinfo.BoxInfo)) return -1;
  } else {
    // Otherwise, who knows what was read?
    mprinterr("Error: AmberRestart::setupTrajin(): When attempting to read in\n");
    mprinterr("Error: box coords/velocity info got %lu chars, expected 0, 37,\n",lineSize);
    mprinterr("Error: 73, or %lu.\n",file_.FrameSize());
    mprinterr("Error: This usually indicates a malformed or corrupted restart file.\n");
    return -1;
  }

  // Recalculate the frame size
  // NOTE: Need to worry about isDos here?
  if (HasV_)
    file_.SetupFrameBuffer( 2*natom3_, 12, 6, 0 );
  if (tinfo.BoxInfo.Type() != Box::NOBOX)
    file_.ResizeBuffer( numBoxCoords_ );
  tinfo.NreplicaDim = 0;
  tinfo.HasV = HasV_;
  tinfo.HasT = HasT_;
  tinfo.IsSeekable = false;
  tinfo.Title = title_;

  /*if (debug_ > 0) {
    mprintf("\tAmber Restart hasBox=%i hasVelocity=%i numBoxCoords=%i\n",
            (int)hasBox_,(int)hasVelocity_,numBoxCoords_);
    mprintf("\tAmber Restart frameSize= %lu\n",frameSize_);
    if (hasBox_) { 
      mprintf("\tAmber Restart box lengths: %lf %lf %lf\n",
              boxLength_[0],boxLength_[1],boxLength_[2]);
      mprintf("\tAmber Restart box angles: %lf %lf %lf\n", 
              boxAngle_[0],boxAngle_[1],boxAngle_[2]);
    }
  }*/
 
  closeTraj();
  // Only 1 frame in restart by definition
  return 1;
}

// Traj_AmberRestart::readFrame()
/** Get the restart file frame. If velocities are present, read those too.
  */
int Traj_AmberRestart::readFrame(int set,double *X,double *V,double *box, double *T) {
  // Read restart coords into frameBuffer_
  if ( file_.ReadFrame()==-1 ) {
    mprinterr("Error: AmberRestart::readFrame(): Error reading coordinates.\n");
    return 1;
  }
  // Set frame temp
  if (HasT_)
    *T = restartTemp_;
  // Get coords from buffer
  file_.BufferBegin();
  file_.BufferToDouble(X, natom3_);
  // Get velocity from buffer if present
  if (HasV_) 
    file_.BufferToDouble(V, natom3_);
  // Get box from buffer if present
  if (numBoxCoords_!=0) 
    file_.BufferToDouble(box, numBoxCoords_);

  return 0;
}

int Traj_AmberRestart::readVelocity(int set, double* V) {
  if (HasV_) {
    if ( file_.ReadFrame()==-1 ) {
      mprinterr("Error: AmberRestart::readVelocity(): Error reading file.\n");
      return 1;
    }
    // Start buffer right after coords.
    file_.BufferBeginAt(coordSize_);
    file_.BufferToDouble(V, natom3_);
    return 0;
  }
  return 1;
}

// Traj_AmberRestart::writeFrame()
/** Write coords in Frame to file in amber restart format. Always calculate the 
  * frame size since coords may have been stripped from Frame.
  */
int Traj_AmberRestart::writeFrame(int set, double *X, double *V, double *box, double T) {
  // If just writing 1 frame dont modify output filename
  if (singleWrite_) {
    if (file_.OpenWriteWithName( file_.FullFileName() )) return 1;
  } else {
    if (file_.OpenWriteWithName( NumberFilename( file_.FullFileName(), set + 1 ) )) return 1;
  }

  // Write out title
  file_.Printf("%-80s\n",title_.c_str());
  // Write out atoms
  file_.Printf("%5i",restartAtoms_);
  // Write out restart time
  if (time0_>=0) {
    restartTime_ = (double) set;
    restartTime_ += time0_;
    restartTime_ *= dt_;
    file_.Printf("%15.7lE",restartTime_);
  }
  // Write out temperature
  if (HasT_)
    file_.Printf("%15.7lE",T);
  file_.Printf("\n");

  // Write coords to buffer
  file_.BufferBegin();
  file_.DoubleToBuffer(X, natom3_, "%12.7lf");
  // Write velocity to buffer. Check V since velocity not known ahead of time
  if (HasV_ && V!=NULL)
    file_.DoubleToBuffer(V, natom3_, "%12.7lf");
  // Write box to buffer
  if (numBoxCoords_!=0)
    file_.DoubleToBuffer(box, numBoxCoords_, "%12.7lf");

  if (file_.WriteFrame()) return 1;

  file_.CloseFile();

  return 0;
}

// Traj_AmberRestart::info()
void Traj_AmberRestart::info() {
  mprintf("is an AMBER restart file");
  // If read access we know for sure whether there are velocities.
  if (file_.Access() != CpptrajFile::WRITE) {
    if (HasV_)
      mprintf(" with velocity info");
    else
      mprintf(", no velocities");
  // If write, not sure yet whether velocities will be written since
  // it also depends on if the frame has velocity info, so only state
  // if novelocity was specified.
  } else {
    if (!HasV_) mprintf(", no velocities");
  }
}

