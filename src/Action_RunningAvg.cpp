// Action_RunningAvg
#include "Action_RunningAvg.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_RunningAvg::Action_RunningAvg() : 
  Nwindow_(0),
  d_Nwindow_(0),
  frameThreshold_(0),
  currentWindow_(0),
  windowNatom_(0)
{} 

// Action_RunningAvg::init()
/// Expected call: runningaverage [window <value>] 
int Action_RunningAvg::init( ) {
  // Get Keywords
  Nwindow_ = actionArgs.getKeyInt("window",5);
  if (Nwindow_ < 1 ) {
    mprinterr("Error: RunningAvg: window must be >= 1.\n");
    return 1;
  }

  // Reserve space for Nwindow frames
  Window_.resize( Nwindow_ );
  // Frame above which averaging will start
  frameThreshold_ = Nwindow_ - 1;
  currentWindow_ = 0;
  windowNatom_ = 0;
  // For division of frames, cast Nwindow to double
  d_Nwindow_ = (double)Nwindow_;
  // Get Masks
  // Dataset: Set up by adding to the main data set list.

  mprintf("    RUNNINGAVG: Running average of size %i will be performed over input coords.\n",
          Nwindow_);

  return 0;
}

// Action_RunningAvg::setup()
int Action_RunningAvg::setup() {
  // If windowNatom is 0, this is the first setup.
  // If windowNatom is not 0, setup has been called for another parm.
  // Check if the number of atoms has changed. If so the running average
  // will break.
  if ( currentParm->Natom() != windowNatom_ ) {
    if (windowNatom_!=0) {
      mprintf("Warning: # atoms in parm %s different than previous parm.\n",
              currentParm->c_str());
      mprintf("         Running average will NOT be carried over between parms!\n");
      return 1;
    }
    windowNatom_ = currentParm->Natom();
    // Set up a frame for each window, no masses
    for (int i = 0; i < Nwindow_; i++)
      Window_[i].SetupFrame( windowNatom_ );
    // Setup frame to hold average coords
    avgFrame_.SetupFrame( windowNatom_ );
    // Zero avgFrame
    avgFrame_.ZeroCoords();
    // Set up result
    resultFrame_.SetupFrame( windowNatom_ );
  } 

  // Print info for this parm
  mprintf("\tRunning average set up for %i atoms.\n",windowNatom_);
  return 0;  
}

// Action_RunningAvg::action()
int Action_RunningAvg::action() {
  // If frameNum is >= Nwindow, subtract from avgFrame. currentWindow is at
  // the frame that should be subtracted.
  if (frameNum > frameThreshold_) { 
    //mprintf("DBG:\tSubtracting Window[%i] from avgFrame.\n",currentWindow_);
    avgFrame_ -= Window_[currentWindow_];
  }

  // Add current coordinates to avgFrame
  //mprintf("DBG:\tAdding frame %i to avgFrame.\n",frameNum);
  avgFrame_ += *currentFrame;

  // Store current coordinates in Window
  //mprintf("DBG:\tAssigning frame %i to window %i (%i = %i)\n",frameNum,currentWindow_,
  //        Window_[currentWindow_].natom, currentFrame->natom);
  Window_[currentWindow_] = *currentFrame;
  ++currentWindow_;
  // If currentWindow is out of range, reset
  if (currentWindow_==Nwindow_) currentWindow_=0;

  // If not enough frames to average yet return 3 to indicate further
  // processing should be suppressed.
  if (frameNum < frameThreshold_)
    return (int)ACTION_SUPPRESSCOORDOUTPUT;
  // Otherwise there are enough frames to start processing the running average
  else {
    //mprintf("DBG:\tCalculating average for frame %i\n",frameNum); 
    resultFrame_.Divide( avgFrame_, d_Nwindow_ );
    // Set frame
    currentFrame = &resultFrame_;
  }

  return (int)ACTION_OK;
} 

