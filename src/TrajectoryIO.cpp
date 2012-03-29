// TrajectoryIO
#include "TrajectoryIO.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
TrajectoryIO::TrajectoryIO() {
  seekable_ = false;
  hasBox_ = false;
  boxAngle_[0] = 0.0;
  boxAngle_[1] = 0.0;
  boxAngle_[2] = 0.0;
  hasTemperature_ = false;
  hasVelocity_ = false;
}

// COPY CONSTRUCTOR
TrajectoryIO::TrajectoryIO(const TrajectoryIO &rhs) :
  CpptrajFile(rhs)
{
  title_ = rhs.title_;
  seekable_ = rhs.seekable_;
  hasBox_ = rhs.hasBox_;
  boxAngle_[0] = rhs.boxAngle_[0];
  boxAngle_[1] = rhs.boxAngle_[1];
  boxAngle_[2] = rhs.boxAngle_[2];
  hasTemperature_ = rhs.hasTemperature_;
  hasVelocity_ = rhs.hasVelocity_;
}

// Assignment
TrajectoryIO &TrajectoryIO::operator=(const TrajectoryIO &rhs) {
  // Self
  if (this == &rhs) return *this;
  // Base class 
  CpptrajFile::operator=(rhs);
  // Deallocate
  // Allocate and copy
  title_ = rhs.title_;
  seekable_ = rhs.seekable_;
  hasBox_ = rhs.hasBox_;
  boxAngle_[0] = rhs.boxAngle_[0];
  boxAngle_[1] = rhs.boxAngle_[1];
  boxAngle_[2] = rhs.boxAngle_[2];
  hasTemperature_ = rhs.hasTemperature_;
  hasVelocity_ = rhs.hasVelocity_;
  return *this;
}

// DESTRUCTOR
TrajectoryIO::~TrajectoryIO() {
}

// TrajectoryIO::SetDebug()
/** Set debug level. */
void TrajectoryIO::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("TrajectoryIO debug level set to %i\n",debug_);
}

// TrajectoryIO::SetTitle()
/** Set title for this trajfile. */
void TrajectoryIO::SetTitle(char *titleIn) {
  //mprintf("DEBUG: Attempting to set title for %s: [%s]\n",trajfilename,titleIn);
  if (titleIn==NULL) return;
  title_.assign(titleIn);
  // Remove terminal newline if one exists
  std::string::iterator lastchar = title_.end();
  --lastchar;
  if (*lastchar=='\n')
    title_.erase( lastchar );
}

// TrajectoryIO::SetTemperature()
/** Indicate temperature information should be read/written. */
void TrajectoryIO::SetTemperature() {
  hasTemperature_ = true;
}

// TrajectoryIO::SetBox()
void TrajectoryIO::SetBox() {
  hasBox_=true;
}

// TrajectoryIO::CheckBoxInfo()
int TrajectoryIO::CheckBoxInfo(AmberParm *trajParm) {
  BoxType boxType;
  if (trajParm==NULL) return 1;
  if (!hasBox_) return 0;
  // If box coords present but no box info in associated parm, print
  // a warning.
  if (trajParm->boxType == NOBOX) {
    mprintf("\tWarning: Box info present in trajectory %s but not in\n",BaseName());
    mprintf("\t         associated parm %s\n",trajParm->parmName);
  }
  boxType = CheckBoxType(boxAngle_, debug_);
  // If box coords present but returned box type is NOBOX then no angles 
  // present in trajectory. Set box angles to parm default. If parm has
  // no box information then assume orthorhombic box.
  // NOTE: Is this only good for amber trajectory?
  if (boxType == NOBOX) {
    if (trajParm->boxType == NOBOX) {
      //mprinterr("Error: No angle information present in trajectory %s\n",trajName);
      //mprinterr("       or parm %s.\n",trajParm->parmName);
      //return 1;
      mprintf("\tWarning: No angle information present in trajectory %s\n",BaseName());
      mprintf("\t         or parm %s - setting angles to 90.0!\n",trajParm->parmName);
      boxAngle_[0] = 90.0;
      boxAngle_[1] = 90.0;
      boxAngle_[2] = 90.0;
    } else {
      if (debug_>0) {
        mprintf("\tWarning: No angle information present in trajectory %s:\n",BaseName());
        mprintf("\t         Using angles from parm %s (beta=%lf).\n",trajParm->parmName,
                trajParm->Box[3]);
      }
      boxAngle_[0] = trajParm->Box[3];
      boxAngle_[1] = trajParm->Box[4];
      boxAngle_[2] = trajParm->Box[5];
    }
    // Set trajectory box type from angles in boxAngle
    boxType = CheckBoxType(boxAngle_, debug_);
  }
  if (debug_>0 || trajParm->boxType == NOBOX) {
    mprintf("\t[%s] Box type is",BaseName());
    if (boxType==NOBOX) mprintf(" None.\n");
    else if (boxType==ORTHO) mprintf(" Orthorhombic.\n");
    else if (boxType==NONORTHO) mprintf(" NonOrthorhombic.\n");
  }
  // If no box info in parm, set it from trajectory
  if (trajParm->boxType == NOBOX) {
    mprintf("\tWarning: Setting parm %s box information from trajectory %s.\n",
            trajParm->parmName,BaseName());
    trajParm->boxType = boxType;
    trajParm->Box[3] = boxAngle_[0];
    trajParm->Box[4] = boxAngle_[1];
    trajParm->Box[5] = boxAngle_[2];
  }
  return 0;
}

bool TrajectoryIO::Seekable() {
  return seekable_;
}

bool TrajectoryIO::HasBox() {
  return hasBox_;
}

bool TrajectoryIO::HasTemperature() {
  return hasTemperature_;
}

bool TrajectoryIO::HasVelocity() {
  return hasVelocity_;
}

