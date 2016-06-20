#include "OutputTrajCommon.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
OutputTrajCommon::OutputTrajCommon() :
  trajParm_(0),
  NframesToWrite_(-1),
  numFramesWritten_(0),
  writeFormat_(TrajectoryFile::UNKNOWN_TRAJ),
  noBox_(false),
  noVel_(false),
  noTemp_(false),
  noTime_(false),
  noFrc_(false),
  noReps_(false),
  append_(false),
  hasRange_(false)
{}

int OutputTrajCommon::CommonTrajoutSetup(FileName const& tnameIn, ArgList& argIn,
                                         TrajectoryFile::TrajFormatType fmtIn)
{
  // Set file name 
  trajName_ = tnameIn;
  // Check for append keyword
  append_ = argIn.hasKey("append");
  // Get specified title if any.
  title_ = argIn.GetStringKey("title");
  // Check for noX arguments - will override info in given CoordinateInfo.
  noBox_ = argIn.hasKey("nobox");
  noVel_ = argIn.hasKey("novelocity");
  noTemp_ = argIn.hasKey("notemperature");
  noTime_ = argIn.hasKey("notime");
  noFrc_ = argIn.hasKey("noforce");
  noReps_ = argIn.hasKey("noreplicadim");
  // If a write format was not specified (UNKNOWN_TRAJ) check the argument
  // list to see if format was specified there.
  writeFormat_ = fmtIn;
  if (writeFormat_==TrajectoryFile::UNKNOWN_TRAJ) {
    writeFormat_ = TrajectoryFile::GetFormatFromArg(argIn);
    // If still UNKNOWN_TRAJ this means no type specified. Check to see if
    // the filename extension is recognized.
    if (writeFormat_ == TrajectoryFile::UNKNOWN_TRAJ) {
      writeFormat_ = TrajectoryFile::GetTypeFromExtension( trajName_.Ext() );
      // Default to Amber trajectory.
      if (writeFormat_ == TrajectoryFile::UNKNOWN_TRAJ) {
        mprintf("Warning: Format not specified and extension '%s' not recognized."
                " Defaulting to Amber Trajectory.\n", trajName_.ext());
        writeFormat_ = TrajectoryFile::AMBERTRAJ;
      }
    }
  }
  // Get a frame range for trajout
  std::string onlyframes = argIn.GetStringKey("onlyframes");
  if (!onlyframes.empty()) {
    if ( FrameRange_.SetRange(onlyframes) )
      mprintf("Warning: trajout %s: onlyframes: %s is not a valid range.\n",
              trajName_.full(), onlyframes.c_str());
    else {
      FrameRange_.PrintRange("\tSaving frames",0);
      mprintf("\n");
    }
    // User frame args start from 1. Start from 0 internally.
    FrameRange_.ShiftBy(-1);
    hasRange_ = true;
  } else {
    if (frameCount_.InitFrameCounter( argIn )) return 1;
    hasRange_ = false;
  }

  return 0;
}

/** If the given file exists and will be appended to, check that the specified
  * format matches the existing format. If not, change to existing format.
  * \return 0 if file exists, 1 if file does not exist.
  */
int OutputTrajCommon::CheckAppendFormat(FileName const& fname,
                                        TrajectoryFile::TrajFormatType& writeFormat)
{
  if (File::Exists(fname)) {
    TrajectoryFile::TrajFormatType appendFormat;
    TrajectoryIO* tio = TrajectoryFile::DetectFormat( fname, appendFormat );
    if (appendFormat ==  TrajectoryFile::UNKNOWN_TRAJ)
      mprintf("Warning: Could not determine file format for 'append'. Using %s\n",
              TrajectoryFile::FormatString( writeFormat ) );
    else {
      if (writeFormat != TrajectoryFile::UNKNOWN_TRAJ && writeFormat != appendFormat)
        mprintf("Warning: Specified format %s for %s but file exists and is %s\n",
                TrajectoryFile::FormatString(writeFormat), fname.full(),
                TrajectoryFile::FormatString(appendFormat));
      writeFormat = appendFormat;
    }
    delete tio;
  } else {
    mprintf("Warning: 'append' specified for non-existent file.\n");
    return 1;
  }
  return 0;
}

int OutputTrajCommon::SetupCoordInfo(Topology* tparmIn, int nFrames, CoordinateInfo const& cInfoIn)
{
  if (tparmIn == 0) return 1;
  trajParm_ = tparmIn;
  // Set coordinate info. Remove box, velocities, etc if specified.
  cInfo_ = cInfoIn;
  if (noBox_ ) cInfo_.SetBox( Box() );
  if (noVel_ ) cInfo_.SetVelocity( false );
  if (noTemp_) cInfo_.SetTemperature( false );
  if (noTime_) cInfo_.SetTime( false );
  if (noFrc_ ) cInfo_.SetForce( false );
  if (noReps_) cInfo_.SetReplicaDims( ReplicaDimArray() );
  // Determine how many frames will be written
  NframesToWrite_ = nFrames;
  if (hasRange_)
    NframesToWrite_ = FrameRange_.Size();
  //trajIsOpen_ = true;
  // If a framerange is defined set it to the beginning of the range
  if (hasRange_)
    rangeframe_ = FrameRange_.begin();
  numFramesWritten_ = 0;
  return 0;
}

void OutputTrajCommon::CommonInfo() const {
  if (trajParm_ != 0) mprintf(", Parm %s", trajParm_->c_str());
  if (noBox_) mprintf(" no box info,");
  if (noVel_) mprintf(" no velocities,");
  if (noTemp_) mprintf(" no temperatures,");
  if (noTime_) mprintf(" no times,");
  if (noFrc_) mprintf(" no forces,");
  if (noReps_) mprintf(" no replica dimensions,");
  if (hasRange_)
    FrameRange_.PrintRange(": Writing frames", 1);
  else if (NframesToWrite_ > 0) {
    mprintf(": Writing %i frames", NframesToWrite_);
    frameCount_.FrameCounterBrief();
  }
  if (append_) mprintf(", appended");
  mprintf("\n");
}
