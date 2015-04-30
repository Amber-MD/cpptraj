#include "Trajout.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
Trajout::Trajout() :
  numFramesProcessed_(0),
  trajIsOpen_(false),
  nobox_(false),
  append_(false),
  hasRange_(false),
  rangeframe_(FrameRange_.end())
{}

// Trajout::CommonTrajoutSetup()
/** Write is set up for topology only when first frame written. */
int Trajout::CommonTrajoutSetup(std::string const& tnameIn, ArgList& argIn, Topology* tparmIn, 
                                TrajectoryFile::TrajFormatType& writeFormat)
{
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Mark as not yet open
  trajIsOpen_ = false;
  // Check for append keyword
  append_ = argIn.hasKey("append");
  // Set file name 
  SetTrajFileName( tnameIn, false );
  // If a write format was not specified (UNKNOWN_TRAJ) check the argument
  // list to see if format was specified there.
  if (writeFormat==UNKNOWN_TRAJ) {
    writeFormat = GetFormatFromArg(argIn);
    // If still UNKNOWN_TRAJ this means no type specified. Check to see if
    // the filename extension is recognized.
    if (writeFormat == UNKNOWN_TRAJ) {
      writeFormat = GetTypeFromExtension( TrajFilename().Ext() );
      // Default to Amber trajectory.
      if (writeFormat == UNKNOWN_TRAJ) {
        mprintf("Warning: Format not specified and extension '%s' not recognized."
                " Defaulting to Amber Trajectory.\n", TrajFilename().ext());
        writeFormat = AMBERTRAJ;
      }
    }
  }
  // Get specified title if any.
  title_ = argIn.GetStringKey("title");
  // Get a frame range for trajout
  std::string onlyframes = argIn.GetStringKey("onlyframes");
  if (!onlyframes.empty()) {
    if ( FrameRange_.SetRange(onlyframes) )
      mprintf("Warning: trajout %s: onlyframes: %s is not a valid range.\n",
              TrajFilename().full(), onlyframes.c_str());
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
  // Check for nobox argument - will override any box info present in parm
  // when trajectory IO is set up.
  nobox_ = argIn.hasKey("nobox");

  return 0;
}

int Trajout::FirstFrameSetup(std::string const& trajoutName, TrajectoryIO* trajio, 
                             Topology* tparmIn)
{
  if (debug_>0) rprintf("\tSetting up %s for WRITE, %i atoms, originally %i atoms.\n",
                          TrajFilename().base(),tparmIn->Natom(),TrajParm()->Natom());
  SetTrajParm( tparmIn );
  // Use parm to set up coord info for the traj. If 'nobox' was specified
  // remove any box info.
  CoordinateInfo cInfo = tparmIn->ParmCoordInfo();
  if (nobox_)
    cInfo.SetBox( Box() );
  // Determine how many frames will be written
  int NframesToWrite = TrajParm()->Nframes();
  if (hasRange_)
    NframesToWrite = FrameRange_.Size();
  // Set up write and open for the current parm file 
  if (trajio->setupTrajout(trajoutName, TrajParm(), cInfo, NframesToWrite, append_))
    return 1;
  if (debug_ > 0)
    Frame::PrintCoordInfo(TrajFilename().base(), tparmIn->c_str(), trajio->CoordInfo()); 
  trajIsOpen_ = true;
  // If a framerange is defined set it to the beginning of the range
  if (hasRange_)
    rangeframe_ = FrameRange_.begin();
  numFramesProcessed_ = 0;
  return 0;
}

/** If the given file exists and will be appended to, check that the specified
  * format matches the existing format. If not, change to existing format.
  */
int Trajout::CheckAppendFormat(std::string const& fname, TrajFormatType& writeFormat)
{
  if (fileExists(fname)) {
    TrajectoryFile::TrajFormatType appendFormat;
    TrajectoryIO* tio = DetectFormat( fname, appendFormat );
    if (appendFormat ==  TrajectoryFile::UNKNOWN_TRAJ)
      mprintf("Warning: Could not determine file format for 'append'. Using %s\n",
              FormatString( writeFormat ) );
    else {
      if (writeFormat != TrajectoryFile::UNKNOWN_TRAJ && writeFormat != appendFormat)
        mprintf("Warning: Specified format %s for %s but file exists and is %s\n",
                FormatString(writeFormat), fname.c_str(), FormatString(appendFormat));
      writeFormat = appendFormat;
    }
    delete tio;
  } else {
    mprintf("Warning: 'append' specified for non-existent file.\n");
    append_ = false;
  }
  return 0;
}

Range Trajout::MembersToWrite(std::string const& onlyMembers, int ensembleSize) const {
  Range members;
  int err;
  // Empty String indicates write All members
  if (onlyMembers.empty())
    err = members.SetRange(0, ensembleSize);
  else
    err = members.SetRange(onlyMembers);
  if (err != 0 || members.Empty()) {
    mprinterr("Error: onlymembers: Invalid range (%s)\n", onlyMembers.c_str());
    return Range();
  }
  return members;
}

void Trajout::CommonInfo(TrajectoryIO* trajio) const {
  trajio->Info();
  mprintf(", Parm %s",TrajParm()->c_str());
  if (nobox_) mprintf(" (no box info)");
  if (hasRange_)
    FrameRange_.PrintRange(": Writing frames", 1);
  else {
    mprintf(": Writing %i frames", TrajParm()->Nframes());
    frameCount_.FrameCounterBrief();
  }

  if (append_) mprintf(", appended");
  mprintf("\n");
}
