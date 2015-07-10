#include "OutputTrajCommon.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
OutputTrajCommon::OutputTrajCommon() :
  trajParm_(0),
  numFramesProcessed_(0),
  writeFormat_(TrajectoryFile::UNKNOWN_TRAJ),
  nobox_(false),
  append_(false),
  hasRange_(false)
{}

int OutputTrajCommon::CommonTrajoutSetup(std::string const& tnameIn, ArgList& argIn,
                                         TrajectoryFile::TrajFormatType fmtIn)
{
  // Set file name FIXME: Proper tilde expansion
  trajName_.SetFileName( tnameIn );
  // Check for append keyword
  append_ = argIn.hasKey("append");
  // Get specified title if any.
  title_ = argIn.GetStringKey("title");
  // Check for nobox argument - will override any box info present in parm
  // when trajectory IO is set up.
  nobox_ = argIn.hasKey("nobox");
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
int OutputTrajCommon::CheckAppendFormat(std::string const& fname,
                                        TrajectoryFile::TrajFormatType& writeFormat)
{
  if (fileExists(fname)) {
    TrajectoryFile::TrajFormatType appendFormat;
    TrajectoryIO* tio = TrajectoryFile::DetectFormat( fname, appendFormat );
    if (appendFormat ==  TrajectoryFile::UNKNOWN_TRAJ)
      mprintf("Warning: Could not determine file format for 'append'. Using %s\n",
              TrajectoryFile::FormatString( writeFormat ) );
    else {
      if (writeFormat != TrajectoryFile::UNKNOWN_TRAJ && writeFormat != appendFormat)
        mprintf("Warning: Specified format %s for %s but file exists and is %s\n",
                TrajectoryFile::FormatString(writeFormat), fname.c_str(),
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

