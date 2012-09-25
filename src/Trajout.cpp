#include "Trajout.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Trajout::Trajout() :
  numFramesProcessed_(0),
  trajio_(0),
  trajIsOpen_(false),
  nobox_(false),
  append_(false),
  hasRange_(false),
  rangeframe_(FrameRange_.end())
{}

// DESTRUCTOR
Trajout::~Trajout() {
  if (trajio_!=0) delete trajio_;
}

// Trajout::SetupTrajWrite()
int Trajout::SetupTrajWrite(std::string const& tnameIn, ArgList *argIn, Topology *tparmIn,
                                   TrajFormatType writeFormatIn)
{
  TrajectoryFile::TrajFormatType writeFormat = writeFormatIn;
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: SetupTrajWrite: No filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Mark as not yet open
  trajIsOpen_ = false;
  // Check for append keyword
  if (argIn!=0) 
    append_ = argIn->hasKey("append");
  // Set up CpptrajFile
  CpptrajFile baseFile;
  if (append_) {
    if (baseFile.SetupAppend( tnameIn, debug_ )) return 1;
  } else {
    if (baseFile.SetupWrite( tnameIn, debug_ )) return 1;
  } 
  // Set file name and base file name
  SetFileNames( tnameIn, baseFile.BaseFileName() );
  // If a write format was not specified (UNKNOWN_TRAJ) check the argument
  // list to see if format was specified there. Defaults to AMBERTRAJ.
  if (writeFormat==UNKNOWN_TRAJ) {
    writeFormat = GetFormatFromArg(*argIn);
    // If still AMBERTRAJ this means no type specified. Check to see if
    // the filename extension is recognized.
    if (writeFormat == AMBERTRAJ) {
      writeFormat = GetTypeFromExtension( baseFile.Extension() );
      if (writeFormat == UNKNOWN_TRAJ) writeFormat = AMBERTRAJ;
    }
  }
  // Set up for the specified format
  trajio_ = SetupTrajectoryIO( writeFormat );
  if (trajio_ == 0) return 1;
  // Set up file information from the baseFile
  trajio_->TrajectoryIO::operator=( baseFile ); // NOTE: Should also set debug 

  // Process additional arguments
  if (argIn != 0) {
    // Get specified title if any - will not set if NULL
    trajio_->SetTitle( argIn->getKeyString("title") );

    // Get a frame range for trajout
    ArgList::ConstArg onlyframes = argIn->getKeyString("onlyframes");
    if (onlyframes!=0) {
      if ( FrameRange_.SetRange(onlyframes) )
        mprintf("Warning: trajout %s: onlyframes: %s is not a valid range.\n",
                FullTrajStr(), onlyframes);
      else
        FrameRange_.PrintRange("\tSaving frames",0);
      // User frame args start from 1. Start from 0 internally.
      FrameRange_.ShiftBy(-1);
      hasRange_ = true;
    }

    // Check for nobox argument - will override any box info present in parm
    // when trajectory IO is set up.
    nobox_ = argIn->hasKey("nobox");

    // Process any write arguments specific to certain formats not related
    // to parm file. Options related to parm file are handled on the first
    // write in WriteFrame.
    if (trajio_->processWriteArgs(argIn)) {
      mprinterr("Error: trajout %s: Could not process arguments.\n",FullTrajStr());
      return 1;
    }
  }

  // No more setup here; Write is set up when first frame written.
  return 0;
}

int Trajout::SetupTrajWriteWithArgs(std::string const& tnameIn, const char *argstring,
                                           Topology *tparmIn, TrajFormatType fmtIn)
{
  ArgList tempArg(argstring, " ");
  return SetupTrajWrite(tnameIn,&tempArg,tparmIn,fmtIn);
}

void Trajout::EndTraj() {
  if (trajIsOpen_) {
    trajio_->closeTraj();
    trajIsOpen_ = false;
  }
}

int Trajout::WriteFrame(int set, Topology *tparmIn, Frame &FrameOut) {
  // Check that input parm matches setup parm - if not, skip
  if (tparmIn->Pindex() != TrajParm()->Pindex()) return 0;

  // First frame setup - set up for the input parm, not necessarily the setup
  // parm; this allows things like atom strippping, etc. A stripped parm will
  // have the same pindex as the original parm.
  if (!trajIsOpen_) {
    if (debug_>0) rprintf("\tSetting up %s for WRITE, %i atoms, originally %i atoms.\n",
                          BaseTrajStr(),tparmIn->Natom(),TrajParm()->Natom());
    SetTrajParm( tparmIn );
    // Use parm to set up box info for the traj unless nobox was specified.
    // If box angles are present in traj they will be used instead.
    if (!nobox_) {
      if (TrajParm()->BoxType()!=Box::NOBOX)
        trajio_->SetBox();
    }
    // Determine how many frames will be written
    int NframesToWrite = TrajParm()->Nframes();
    if (hasRange_)
      NframesToWrite = FrameRange_.Size();
    // Set up write for the current parm file 
    if (trajio_->setupTrajout(TrajParm(), NframesToWrite)) return 1;
    // Open output traj and mark as set up.
    if (trajio_->openTraj()) return 1;
    trajIsOpen_ = true;
    // If a framerange is defined set it to the begining of the range
    if (hasRange_)
      rangeframe_ = FrameRange_.begin();
    numFramesProcessed_ = 0;
  }
  // If there is a framerange defined, check if this frame matches. If so,
  // write this frame and increment to the next frame in the range.
  if (hasRange_) {
    // If no more frames in the framerange, skip.
    if (rangeframe_ == FrameRange_.end()) return 0;
    // If this frame is not the next in the range, skip.
    if ( *rangeframe_ != set ) return 0;
    // This frame is next in the range. Advance FrameRange iterator.
    ++rangeframe_;
  }

  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (trajio_->writeFrame(set, FrameOut.xAddress(), FrameOut.vAddress(),
                               FrameOut.bAddress(), FrameOut.Temperature())) return 1;
  ++numFramesProcessed_;

  return 0;
}

void Trajout::PrintInfo(int showExtended) {
  mprintf("  [%s] ",BaseTrajStr());
  trajio_->info();
  mprintf(", Parm %s",TrajParm()->c_str());
  if (trajio_->HasBox()) mprintf(" (with box info)");
  if (hasRange_)
    FrameRange_.PrintRange(": Writing frames",OUTPUTFRAMESHIFT);
  else
    mprintf(": Writing %i frames", TrajParm()->Nframes());
  if (append_) mprintf(", appended");
}

