#include "Trajin.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

// CONSTRUCTOR
Trajin::Trajin() :
  start_(0),
  stop_(-1),
  offset_(1),
  total_frames_(0),
  total_read_frames_(-1),
  currentFrame_(0),
  targetSet_(0),
  frameskip_(1),
  numFramesProcessed_(0),
  useProgress_(true)
{}

// Trajin::SetupTrajIO()
int Trajin::SetupTrajIO( std::string const& fname, TrajectoryIO& trajio, ArgList* argIn ) {
  // -1 indicates an error.
  // -2 indicates the number of frames could not be determined, read to EOF.
  total_frames_ = trajio.setupTrajin(fname, TrajParm());
  if (total_frames_ == TrajectoryIO::TRAJIN_ERR) {
    mprinterr("Error: Could not set up %s for reading.\n", fname.c_str());
    return 1;
  }
  if (total_frames_>-1)
    mprintf("\t[%s] contains %i frames.\n", TrajFilename().base(), total_frames_);
  else
    mprintf("\t[%s] contains an unknown number of frames.\n",TrajFilename().base());
  // Set stop based on calcd number of frames.
  if (total_frames_==0) {
    mprinterr("Error: trajectory %s contains no frames.\n",TrajFilename().base());
    return 1;
  }
  if (total_frames_>0)
    stop_ = total_frames_;
  else
    stop_ = -1;
  // Set the start, stop, and offset args based on user input. Do some bounds
  // checking. 
  if (argIn != 0) {
    // Frame args start at 1. Internal frame #s start at 0. 
    // So for a traj with 10 frames:
    // - Internal #: 0 1 2 3 4 5 6 7 8 9
    // - Frame Arg#: 1 2 3 4 5 6 7 8 9 10
    // - Defaults: startArg=1, stopArg=-1, offsetArg=1
    int startArg, stopArg, offsetArg;
    // lastframe is a special case where only the last frame will be selected
    if (argIn->hasKey("lastframe")) {
      if (total_frames_>0) {
        startArg = total_frames_;
        stopArg = total_frames_;
        offsetArg = 1;
      } else {
        mprinterr("Error: [%s] lastframe specified but # frames could not be determined.\n",
                  TrajFilename().base());
        return 1;
      }
    } else {
      startArg = argIn->getNextInteger(1);
      if (argIn->hasKey("last"))
        stopArg = -1;
      else
        stopArg = argIn->getNextInteger(-1);
      offsetArg = argIn->getNextInteger(1);
    }
    if (startArg!=1) {
      if (startArg<1) {
        mprintf("    Warning: %s start argument %i < 1, setting to 1.\n",TrajFilename().base(),startArg);
        start_ = 0; // cpptraj = ptraj - 1
      } else if (total_frames_>=0 && startArg>total_frames_) {
        // If startArg==stopArg and is greater than # frames, assume we want
        // the last frame (useful when reading for reference structure).
        if (startArg==stopArg) {
          mprintf("    Warning: %s start %i > #Frames (%i), setting to last frame.\n",
                  TrajFilename().base(),startArg,total_frames_);
          start_ = total_frames_ - 1;
        } else {
          mprinterr("Error: [%s] start %i > #Frames (%i), no frames will be processed.\n",
                  TrajFilename().base(),startArg,total_frames_);
          //start=startArg - 1;
          return 1;
        }
      } else
        start_ = startArg - 1;
    }
    if (stopArg!=-1) {
      if ((stopArg - 1)<start_) { // cpptraj = ptraj - 1
        mprinterr("Error: [%s] stop %i < start, no frames will be processed.\n",
                TrajFilename().base(),stopArg);
        //stop = start;
        return 1;
      } else if (total_frames_>=0 && stopArg>total_frames_) {
        mprintf("    Warning: %s stop %i >= #Frames (%i), setting to max.\n",
                TrajFilename().base(),stopArg,total_frames_);
        stop_ = total_frames_;
      } else
        stop_ = stopArg;
    }
    if (offsetArg!=1) {
      if (offsetArg<1) {
        mprintf("    Warning: %s offset %i < 1, setting to 1.\n",
                TrajFilename().base(),offsetArg);
        offset_ = 1;
      } else if (stop_!=-1 && offsetArg > stop_ - start_) {
        mprintf("    Warning: %s offset %i is so large that only 1 set will be processed.\n",
                TrajFilename().base(),offsetArg);
        offset_ = offsetArg;
      } else
        offset_ = offsetArg;
    }
    if (debug_>0)
      mprintf("DEBUG [%s] SetArgs: Start %i Stop %i  Offset %i\n",
              TrajFilename().base(),start_,stop_,offset_);
  }
  return 0;
}

int Trajin::CheckBoxInfo(const char* parmName, Box& parmBox, Box const& trajBox) {
  if (!trajBox.HasBox()) {
    if ( parmBox.HasBox()) {
      // No box in traj but box in parm - disable parm box.
      mprintf("Warning: Box information present in parm but not in trajectory.\n");
      mprintf("Warning: DISABLING BOX in parm [%s]!\n", parmName);
      parmBox.SetNoBox();
      return 0;
    } else
      // No box in traj or parm - exit.
      return 0;
  }
  // Check for zero box lengths
  if ( trajBox.BoxX() < SMALL || trajBox.BoxY() < SMALL || trajBox.BoxZ() < SMALL ) {
    mprintf("Warning: Box information present in trajectory but lengths are zero.\n");
    mprintf("Warning: DISABLING BOX in parm [%s]!\n", parmName);
    parmBox.SetNoBox();
    return 0;
  }
  // If box coords present but no box info in associated parm, print
  // a warning. Set parm box from trajectory
  if (!parmBox.HasBox()) {
    mprintf("Warning: Box info present in trajectory %s but not in\n", TrajFilename().base());
    mprintf("Warning: associated parm %s\n", parmName);
    mprintf("Warning: Setting parm box information from trajectory.\n");
    parmBox = trajBox;
  } else if (parmBox.Type() != trajBox.Type()) {
    mprintf("Warning: Trajectory box is %s but parm box is %s\n", 
            trajBox.TypeName(), parmBox.TypeName());
    mprintf("Warning: Setting parm box information from trajectory.\n");
    parmBox = trajBox;
  }
  return 0;
}

// Trajin::setupFrameInfo()
/** Calculate number of frames that will be read based on start, stop, and
  * offset (total_read_frames). 
  * \return the total number of frames that will be read for this traj.
  * \return -1 if the number of frames could not be determined.
  */
int Trajin::setupFrameInfo() {
  int Nframes;
  int ptraj_start_frame, ptraj_end_frame;
  int traj_start_frame, traj_end_frame;
  // DEBUG - No mpi for now
  int worldrank = 0;
  int worldsize = 1;

  //mprintf("DEBUG: Calling setupFrameInfo for %s with %i %i %i\n",trajName,
  //        start,stop,offset);
  if (stop_==-1) return -1;

  // Calc total frames that will be read
  Nframes = stop_ - start_;
  total_read_frames_ = Nframes / offset_;
  // Round up
  if ( (Nframes % offset_) > 0 )
    ++total_read_frames_;
  //divresult = div( (stop - start), offset);
  //total_read_frames = divresult.quot;
  //if (divresult.rem!=0) total_read_frames++;

  // Calc min num frames read by any given thread
  // last thread gets leftovers
  // In case of 0, last thread gets the frame
  Nframes = total_read_frames_ / worldsize;
  int leftover_frames = total_read_frames_ % worldsize;
  //divresult = div(total_read_frames,worldsize);
  //Nframes=divresult.quot;

  // Ptraj (local) start and end frame
  ptraj_start_frame = (worldrank*Nframes);
  ptraj_end_frame = ptraj_start_frame + Nframes;
  // Last thread gets the leftovers
  if (worldrank == worldsize-1)
    ptraj_end_frame += leftover_frames;

  // Actual Traj start and end frame (for seeking)
  traj_start_frame = (ptraj_start_frame * offset_) + start_;
  traj_end_frame = ((ptraj_end_frame-1) * offset_) + start_;

  start_ = traj_start_frame;
  stop_ = traj_end_frame;

  if ( total_read_frames_ == 0) {
    mprinterr("Error: No frames will be read from %s based on start, stop,\n",TrajFilename().base());
    mprinterr("         and offset values (%i, %i, %i)\n",start_+1,stop_+1,offset_);
  }

  return total_read_frames_;
}

// Trajin::SingleFrame()
/** Tell the trajectory to set up stop and offset so that only start frame
  * will be processed.
  */
void Trajin::SingleFrame() {
  stop_ = start_ + 1;
  offset_ = 1;
  // Call setupFrameInfo to recalc total_read_frames. Since setupFrameInfo 
  // should have already been called in SetupRead (and thus any errors 
  // handled there) dont check for an error here. It should return 1.
  if ( setupFrameInfo() != 1 ) {
    mprintf("  Warning: Single frame requested for %s but not calcd!\n",TrajFilename().base());
    mprintf("           start/stop/offset (%i, %i, %i)\n",start_+1,stop_+1,offset_);
  }
}

// Trajin::PrepareForRead()
void Trajin::PrepareForRead(bool useIn, bool seekable) {
  numFramesProcessed_ = 0;
  // Setup progress bar
  useProgress_ = useIn;
  if (useProgress_) progress_.SetupProgress( total_read_frames_ );
  // Determine what frames will be read
  targetSet_ = start_;
  if (seekable) {
    frameskip_ = offset_;
    currentFrame_ = start_;
  } else {
    frameskip_ = 1;
    currentFrame_ = 0;
  }
}

// Trajin::PrintInfoLine()
void Trajin::PrintInfoLine() {
  if (stop_ != -1)
    rprintf( "----- [%s] (%i-%i, %i) -----\n",TrajFilename().base(),start_+1,stop_+1,offset_);
  else
    rprintf( "----- [%s] (%i-EOF, %i) -----\n",TrajFilename().base(),start_+1,offset_);
}

// Trajin::PrintFrameInfo()
void Trajin::PrintFrameInfo() {
  if (stop_!=-1 && total_frames_>0)
    //mprintf(": %i-%i, %i (reading %i of %i)",start,stop,offset,total_read_frames,total_frames);
    mprintf(" (reading %i of %i)",total_read_frames_,total_frames_);
  else if (stop_!=-1 && total_frames_ < 0)
    mprintf(" (reading %i)",total_read_frames_);
  else
    mprintf(", unknown #frames, start=%i offset=%i",start_,offset_);
}
