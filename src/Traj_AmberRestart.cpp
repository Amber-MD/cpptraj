// AmberRestart
#include <cstdio> // sscanf
#include "Traj_AmberRestart.h"
#include "CpptrajStdio.h"
#include "CharBuffer.h"

const size_t AmberRestart::BUF_SIZE = 83;

// CONSTRUCTOR
AmberRestart::AmberRestart() {
  restartAtoms_=0;
  frameSize_=0;
  frameBuffer_=NULL;
  numBoxCoords_=0;
  restartTime_=0;
  restartTemp_=-1.0;
  // Set seekable since only 1 frame read (i.e. we know size) 
  seekable_=true;
  time0_ = OUTPUTFRAMESHIFT;
  dt_ = 1.0;
  singleWrite_ = false;
}

// AmberRestart::closeTraj()
void AmberRestart::closeTraj() {
  CloseFile();
}

// AmberRestart::openTraj()
/** Open the restart file. Get title, time, restart atoms, temperature
  */
int AmberRestart::openTraj() {
  char buffer[BUF_SIZE];
  int nread; // Dont declare variables inside a switch block

  switch (access_) {

    case READ :
      if (OpenFile()) return 1;
      // Read in title
      if (IO->Gets(buffer,BUF_SIZE)) {
         mprinterr("Error: AmberRestart::open(): Reading restart title.\n");
        return 1;
      }
      title_.assign(buffer);
      // Read in natoms, time, and Replica Temp if present
      if (IO->Gets(buffer,BUF_SIZE)) {
        mprinterr("Error: AmberRestart::open(): Reading restart atoms/time.\n");
        return 1;
      }
      nread = sscanf(buffer,"%i %lE %lE",&restartAtoms_,&restartTime_,&restartTemp_);
      if (nread<1) {
        mprinterr("Error: AmberRestart::open(): Getting restart atoms/time.\n");
        return 1;
      } else if (nread==1) {
        restartTime_=0.0;
        restartTemp_=-1.0;
        hasTemperature_=false;
      } else if (nread==2) {
        restartTemp_=-1.0;
        hasTemperature_=false;
      } else {
        hasTemperature_=true;
      }
      if (debug_>0) 
        mprintf("  Amber restart: Atoms=%i Time=%lf Temp=%lf\n",restartAtoms_,
                restartTime_, restartTemp_);
      break;

    case APPEND :
      mprinterr("Error: Append not supported for Amber Restart files.\n");
      return 1;
      break;

    case WRITE :
      // Set up title
      if (title_.empty()) {
        title_.assign("Cpptraj Generated Restart");
        title_.resize(80, ' ');
      } else {
        size_t titleSize = title_.size();
        if (titleSize>80) {
          title_.resize(80);
          mprintf("Warning: Amber restart title for %s too long: truncating.\n[%s]\n",
              BaseName(), title_.c_str());
        }
      }
  }

  return 0; 
}

// AmberRestart::SetNoVelocity()
void AmberRestart::SetNoVelocity() {
  hasVelocity_=false;
}

// AmberRestart::processWriteArgs()
int AmberRestart::processWriteArgs(ArgList *argIn) {
  // For write, assume we want velocities unless specified
  hasVelocity_=true;
  if (argIn->hasKey("novelocity")) this->SetNoVelocity();
  time0_ = argIn->getKeyDouble("time0", OUTPUTFRAMESHIFT);
  if (argIn->hasKey("remdtraj")) this->SetTemperature();
  dt_ = argIn->getKeyDouble("dt",1.0);
  return 0;
}

// AmberRestart::setupTrajout()
/** Allocate a character buffer based on number of coords and whether 
  * velocities/box info is present.
  */
int AmberRestart::setupTrajout(Topology *trajParm) {
  size_t frame_lines;

  restartAtoms_ = trajParm->natom;
  natom3_ = trajParm->natom * 3;
  // Calculate the length of coordinate frame in bytes
  frame_lines = ((size_t)natom3_) / 6;
  if ((natom3_ % 6) > 0)
    ++frame_lines;
  frameSize_ = (((size_t)natom3_ * 12) + frame_lines);

  // Dont know ahead of time if velocities will be used, allocate space
  // just in case. Veolcity will not be written if V input is NULL.
  frameSize_ += frameSize_;

  // If box coords are present, allocate extra space for them
  if (hasBox_) {
    numBoxCoords_ = 6;
    frameSize_ += (((size_t)numBoxCoords_*12)+1);
  }

  // Allocate memory to buffer 1 frame
  // One extra char for NULL
  frameBuffer_ = new char[ frameSize_+1 ];

  // If number of frames to write == 1 set singleWrite so we dont append
  // frame # to filename.
  if (trajParm->parmFrames==1) singleWrite_=true;

  return 0;
}

// AmberRestart::getBoxAngles()
/** Based on input buffer, determine num box coords and get box angles.
  * If successful set hasBox to true.
  */
int AmberRestart::getBoxAngles(char *boxline) {
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
            BaseName());
    mprintf("         Assuming no box information present.\n");
    hasBox_ = false;

  } else if (numBoxCoords_==3) {
    // Lengths read but no angles. Set angles to 0.0, which indicates
    // the prmtop angle should be used in TrajectoryFile
    boxAngle_[0]=0.0;
    boxAngle_[1]=0.0;
    boxAngle_[2]=0.0;
    hasBox_ = true;
  } else if (numBoxCoords_==6) {
    boxAngle_[0]=box[3];
    boxAngle_[1]=box[4];
    boxAngle_[2]=box[5];
    hasBox_ = true;
  } else {
    mprintf("Warning: AmberRestart::getBoxAngles():\n");
    mprintf("         Expect only 3 or 6 box coords in box coord line, got %i.\n",numBoxCoords_);
    hasBox_ = false;
  }
  return 0;
}

// AmberRestart::setupTrajin()
/** Set up amber restart file for reading. Check that number of atoms matches
  * number of atoms in associated parmtop. Check for box/velocity info.
  */
int AmberRestart::setupTrajin(Topology *trajParm) {
  char buffer[BUF_SIZE];
  size_t frame_lines,lineSize;

  if (openTraj()) return -1; // Gets title, time, natoms, and temp if present

  // Check that natoms matches parm natoms
  if (restartAtoms_!=trajParm->natom) {
    mprinterr("Error: Number of atoms in Amber Restart %s (%i) does not\n",
              BaseName(), restartAtoms_);
    mprinterr("       match number in associated parmtop (%i)\n",trajParm->natom);
    return -1;
  }
  natom3_ = restartAtoms_ * 3;

  // Calculate the length of coordinate frame in bytes
  frame_lines = (size_t)natom3_ / 6;
  if ((natom3_ % 6) > 0)
    ++frame_lines;
  // For DOS files CR present before newline
  if (isDos_) frame_lines *= 2;
  frameSize_ = (((size_t)natom3_ * 12) + frame_lines);
  frameBuffer_ = new char[ frameSize_ ];
  //if (debug_>0) mprintf("    Amber Restart frameSize= %i\n",frameSize);

  // Read past restart coords 
  if ( IO->Read(frameBuffer_,frameSize_)==-1 ) {
    mprinterr("Error: AmberRestart::SetupRead(): Error reading coordinates.\n");
    return -1; 
  }

  // Attempt a second read to get velocities or box coords
  lineSize = (size_t)IO->Read(frameBuffer_,frameSize_);
  //mprintf("DEBUG: Restart lineSize on second read = %i\n",lineSize);
  hasVelocity_ = false;
  hasBox_ = false;

  // If 0 or -1 no box or velo 
  if (lineSize<=0) {
    hasVelocity_=false;
    hasBox_ = false; 

  // If filled framebuffer again, has velocity info. 
  } else if (lineSize==frameSize_) {
    hasVelocity_=true;
    // If we can read 1 more line after velocity, should be box info.
    if (IO->Gets(buffer,BUF_SIZE)==0) {
      if (getBoxAngles(buffer)) return -1;
    } 

  // If we read something but didnt fill framebuffer, should have box coords.
  } else if (lineSize<82) {
    hasVelocity_=false;
    // NOTE: OK to NULL terminate frameBuffer here?
    frameBuffer_[lineSize]='\0';
    if (getBoxAngles(frameBuffer_)) return -1;

  // Otherwise, who knows what was read?
  } else {
    mprinterr("Error: AmberRestart::SetupRead(): When attempting to read in\n");
    mprinterr("       box coords/velocity info got %lu chars, expected 0, 37,\n",lineSize);
    mprinterr("       73, or %lu.\n",frameSize_);
    mprinterr("       This usually indicates a malformed or corrupted restart file.\n");
    return -1;
  }

  // Recalculate the frame size
  // NOTE: Need to worry about isDos here?
  if (hasVelocity_) frameSize_ += frameSize_;
  if (hasBox_) frameSize_ += ( ((size_t)numBoxCoords_*12) + 1 );
  delete[] frameBuffer_;
  frameBuffer_ = new char[ frameSize_ ];

  if (debug_ > 0) {
    mprintf("    Amber Restart hasBox=%i hasVelocity=%i numBoxCoords=%i\n",
            (int)hasBox_,(int)hasVelocity_,numBoxCoords_);
    mprintf("    Amber Restart frameSize= %lu\n",frameSize_);
    if (hasBox_) 
      mprintf("    Amber Restart box: %lf %lf %lf\n", boxAngle_[0],boxAngle_[1],boxAngle_[2]);
  }
 
  closeTraj();
  // Only 1 frame in restart by definition
  return 1;
}

// AmberRestart::readFrame()
/** Get the restart file frame. If velocities are present, read those too.
  */
int AmberRestart::readFrame(int set,double *X,double *V,double *box, double *T) {
  // Read restart coords into frameBuffer_
  if ( IO->Read(frameBuffer_,frameSize_)==-1 ) {
    mprinterr("Error: AmberRestart::getFrame(): Error reading coordinates.\n");
    return 1;
  }

  // Set frame temp
  if (hasTemperature_)
    *T = restartTemp_;

  // Get coords from buffer
  BufferBegin();
  BufferToDouble(X, natom3_, 12);
  // Get velocity from buffer if present
  if (hasVelocity_ && V!=NULL) 
    BufferToDouble(V, natom3_, 12);
  // Get box from buffer if present
  if (hasBox_) 
    BufferToDouble(box, numBoxCoords_, 12);

  return 0;
}

// AmberRestart::writeFrame()
/** Write coords in Frame to file in amber restart format. Always calculate the 
  * frame size since coords may have been stripped from Frame.
  */
int AmberRestart::writeFrame(int set, double *X, double *V, double *box, double T) {
  char buffer[1024];

  // If just writing 1 frame dont modify output filename
  if (singleWrite_) {
    if (IO->Open(Name(),"wb")) return 1;
  } else {
    NumberFilename(buffer,(char*)Name(),set + OUTPUTFRAMESHIFT);
    if (IO->Open(buffer,"wb")) return 1;
  }

  // Write out title
  Printf("%-80s\n",title_.c_str());
  // Write out atoms
  Printf("%5i",restartAtoms_);
  // Write out restart time
  if (time0_>=0) {
    restartTime_ = (double) set;
    restartTime_ += time0_;
    restartTime_ *= dt_;
    Printf("%15.7lE",restartTime_);
  }
  // Write out temperature
  if (hasTemperature_)
    Printf("%15.7lE",T);
  Printf("\n");

  // Write coords to buffer
  BufferBegin();
  DoubleToBuffer(X, natom3_, "%12.7lf", 12, 6);
  // Write velocity to buffer. Check V since velocity not known ahead of time
  if (hasVelocity_ && V!=NULL)
    DoubleToBuffer(V, natom3_, "%12.7lf", 12, 6);
  // Write box to buffer
  if (hasBox_)
    BoxToBuffer(box, numBoxCoords_, "%12.7lf",12);

  frameSize_ = (size_t) (bufferPosition_ - frameBuffer_);

  if (IO->Write(frameBuffer_,frameSize_)) return 1;

  // Since when writing amber restarts a number is appended to the filename
  // dont use the CloseFile function in CpptrajFile, just close.
  IO->Close();

  return 0;
}

// AmberRestart::info()
void AmberRestart::info() {
  mprintf("is an AMBER restart file");
  // If read access we know for sure whether there are velocities.
  if (access_!=WRITE) {
    if (hasVelocity_)
      mprintf(" with velocity info");
    else
      mprintf(", no velocities");

  // If write, not sure yet whether velocities will be written since
  // it also depends on if the frame has velocity info, so only state
  // if novelocity was specified.
  } else {
    if (!hasVelocity_) mprintf(", no velocities");
  }
}

