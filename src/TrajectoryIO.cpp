// TrajectoryIO
#include "TrajectoryIO.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

// CONSTRUCTOR
TrajectoryIO::TrajectoryIO() :
  seekable_(false),
  hasBox_(false),
  hasTemperature_(false),
  hasVelocity_(false)
{
  boxAngle_[0] = 0.0;
  boxAngle_[1] = 0.0;
  boxAngle_[2] = 0.0;
  boxLength_[0] = 0.0;
  boxLength_[1] = 0.0;
  boxLength_[2] = 0.0;
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
  boxLength_[0] = rhs.boxAngle_[0];
  boxLength_[1] = rhs.boxAngle_[1];
  boxLength_[2] = rhs.boxAngle_[2];
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
  boxLength_[0] = rhs.boxAngle_[0];
  boxLength_[1] = rhs.boxAngle_[1];
  boxLength_[2] = rhs.boxAngle_[2];
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
int TrajectoryIO::CheckBoxInfo(Topology *trajParm) {
  Box trajbox;
  if (trajParm==NULL) return 1;
  if (!hasBox_) return 0;
  // Check for zero box lengths
  if (boxLength_[0] < SMALL || boxLength_[1] < SMALL || boxLength_[2] < SMALL) {
    mprintf("Warning: Box information present in trajectory but lengths are zero.\n");
    //mprintf("         This will cause imaging to fail unless box information is set.\n");
    mprintf("Warning: DISABLING BOX in parm [%s]!\n",trajParm->c_str());
    trajParm->ParmBox().SetNoBox();
    return 0;
  }
  // If box coords present but no box info in associated parm, print
  // a warning.
  if (trajParm->BoxType() == Box::NOBOX) {
    mprintf("\tWarning: Box info present in trajectory %s but not in\n",BaseName());
    mprintf("\t         associated parm %s\n",trajParm->c_str());
  }
  trajbox.SetAngles( boxAngle_ );
  // If box coords present but returned box type is NOBOX then no angles 
  // present in trajectory. Set box angles to parm default. If parm has
  // no box information then assume orthorhombic box.
  // NOTE: Is this only good for amber trajectory?
  if (trajbox.Type() == Box::NOBOX) {
    if (trajParm->BoxType() == Box::NOBOX) {
      //mprinterr("Error: No angle information present in trajectory %s\n",trajName);
      //mprinterr("       or parm %s.\n",trajParm->parmName);
      //return 1;
      mprintf("\tWarning: No angle information present in trajectory %s\n",BaseName());
      mprintf("\t         or parm %s - setting angles to 90.0!\n",trajParm->c_str());
      boxAngle_[0] = 90.0;
      boxAngle_[1] = 90.0;
      boxAngle_[2] = 90.0;
    } else {
      double temp[6];
      // TODO: Just assign trajbox? Does boxAngle get used anywhere else?
      trajParm->ParmBox().ToDouble( temp );
      if (debug_>0) {
        mprintf("\tWarning: No angle information present in trajectory %s:\n",BaseName());
        mprintf("\t         Using angles from parm %s (beta=%lf).\n",trajParm->c_str(),
                temp[4]);
      }
      boxAngle_[0] = temp[3];
      boxAngle_[1] = temp[4];
      boxAngle_[2] = temp[5];
    }
    trajbox.SetAngles( boxAngle_ );
  }
  if (debug_>0 || trajbox.Type() != Box::NOBOX) {
    mprintf("\t[%s] Box type is %s\n",BaseName(),trajbox.TypeName());
  }
  // If no box info in parm, set it from trajectory
  if (trajParm->BoxType() == Box::NOBOX) {
    mprintf("\tWarning: Setting parm %s box information from trajectory %s.\n",
            trajParm->c_str(),BaseName());
    trajParm->ParmBox().SetAngles( boxAngle_ );
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

