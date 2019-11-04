#include "TrajFrameCounter.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

TrajFrameCounter::TrajFrameCounter() :
  start_(0),
  stop_(-1),
  offset_(1),
  total_frames_(0),
  total_read_frames_(-1),
  current_(0),
  numFramesProcessed_(0)
{}

/** Parse argument list for trajectory-related frame args. Frame args start at
  * 1, internal frame #s start at 0. So for a traj with 10 frames:
  * - Internal #: 0 1 2 3 4 5 6 7 8 9
  * - Frame Arg#: 1 2 3 4 5 6 7 8 9 10
  * - Defaults: start_=1, stop_=-1, offset_=1
  */
int TrajFrameCounter::CheckFrameArgs(int nframes, ArgList &argIn) {
  total_frames_ = nframes;
  if (total_frames_==0) {
    mprinterr("Error: trajectory contains no frames.\n");
    return 1;
  }
  //if (total_frames_>0)
  //  stop_ = total_frames_;
  //else
  //  stop_ = -1;

  if (argIn.hasKey("lastframe")) {
    // lastframe is a special case where only the last frame will be selected
    if (total_frames_ > 0) {
      start_ = total_frames_;
      stop_ = total_frames_;
      offset_ = 1;
    } else {
      mprinterr("Error: lastframe specified but # frames could not be determined.\n");
      return 1;
    }
  } else {
    start_ = argIn.getNextInteger(1);
    // Last explicitly selects final frame as stop arg.
    if (argIn.hasKey("last"))
      stop_ = -1;
    else
      stop_ = argIn.getNextInteger(-1);
    offset_ = argIn.getNextInteger(1);
  }
  // Check that start argument is valid.
  if (start_ != 1) {
    if (start_ == 0) {
      mprintf("Warning: start argument is 0, setting to 1.\n", start_);
      start_ = 1; //start_ = 0;
    } else if (start_ < 0) {
      // Negative start means we want that many frames before stop
      if (stop_ == -1) {
        if (total_frames_ >=0)
          stop_ = total_frames_;
        else {
          mprinterr("Error: For start < 0, stop argument must be specified when # frames unknown.\n");
          return 1;
        }
      }
      mprintf("\tStarting %i frames before frame %i\n", -start_, stop_);
      start_ = stop_ + start_;
      if (start_ < 1) {
        mprintf("Warning: would start before frame 1, setting start to 1.\n");
        start_ = 1;
      }
    } else if (total_frames_ >= 0 && start_ > total_frames_) {
      // start_==stop_ and greater than # frames, archaic 'lastframe'.
      if (start_ == stop_) {
        mprintf("Warning: start %i > #Frames (%i), setting to last frame.\n",
                start_, total_frames_);
        start_ = total_frames_; //start_ = total_frames_ - 1;
      } else {
        mprinterr("Error: start %i > #Frames (%i), no frames will be processed.\n",
                  start_, total_frames_);
        //start=start_ - 1;
        return 1;
      }
    }
  }
  start_--; // Internal frame nums start from 0.
  // Check that stop argument is valid
  if (stop_ != -1) {
    if ( (stop_ - 1) < start_) { // Internal frame nums start from 0.
      mprinterr("Error: stop %i < start, no frames will be processed.\n", stop_);
      //stop = start;
      return 1;
    } else if (total_frames_ >= 0 && stop_ > total_frames_) {
      mprintf("Warning: stop %i > #Frames (%i), setting to max.\n", stop_, total_frames_);
      stop_ = total_frames_;
    } 
  } else if (total_frames_ >= 0) // -1 means use last frame
    stop_ = total_frames_;
  // Check that offset argument is valid.
  if (offset_ != 1) {
    if (offset_ < 1) {
      mprintf("Warning: offset %i < 1, setting to 1.\n", offset_);
      offset_ = 1;
    } else if (stop_ != -1 && offset_ >= (stop_ - start_)) {
      mprintf("Warning: offset %i is so large that only 1 set will be processed.\n",
              offset_);
    }
  }
  //mprintf("DEBUG SetArgs: Start %i Stop %i  Offset %i\n", start_, stop_, offset_);
  // Calculate actual number of frames that will be read based on start,
  // stop, and offset.
  total_read_frames_ = -1;
  if (stop_ != -1) {
    int Nframes = stop_ - start_;
    total_read_frames_ = Nframes / offset_;
    // Round up
    if ( (Nframes % offset_) > 0 )
      ++total_read_frames_;
    if (total_read_frames_ == 0) {
      mprinterr("Error: No frames will be read based on start, stop, "
                "and offset values (%i, %i, %i)\n", start_+1, stop_, offset_);
      return 1;
    }
  }
  return 0;
}

void TrajFrameCounter::PrintFrameInfo() const {
  if (stop_!=-1 && total_frames_>0)
    //mprintf(": %i-%i, %i (reading %i of %i)",start,stop,offset,total_read_frames,total_frames);
    mprintf(" (reading %i of %i)",total_read_frames_,total_frames_);
  else if (stop_!=-1 && total_frames_ < 0)
    mprintf(" (reading %i)",total_read_frames_);
  else
    mprintf(", unknown #frames, start=%i offset=%i",start_,offset_);
}

void TrajFrameCounter::PrintInfoLine(const char* fname) const {
  if (stop_ != -1)
    mprintf( "----- %s (%i-%i, %i) -----\n", fname, start_+1, stop_, offset_);
  else
    mprintf( "----- %s (%i-EOF, %i) -----\n", fname,start_+1,offset_);
}
