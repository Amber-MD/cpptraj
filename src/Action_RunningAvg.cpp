// RunningAvg
#include "Action_RunningAvg.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
RunningAvg::RunningAvg() {
  Nwindow = 0;
  d_Nwindow = 0;
  frameThreshold = 0;
  currentWindow = 0;
  windowNatom = 0;
} 

// RunningAvg::init()
/// Expected call: runningaverage [window <value>] 
int RunningAvg::init( ) {
  // Get Keywords
  Nwindow = actionArgs.getKeyInt("window",5);
  if (Nwindow < 1 ) {
    mprinterr("Error: RunningAvg: window must be >= 1.\n");
    return 1;
  }

  // Reserve space for Nwindow frames
  Window.resize( Nwindow );
  // Frame above which averaging will start
  frameThreshold = Nwindow - 1;
  currentWindow = 0;
  windowNatom = 0;
  // For division of frames, cast Nwindow to double
  d_Nwindow = (double)Nwindow;
  // Get Masks
  // Dataset: Set up by adding to the main data set list.

  mprintf("    RUNNINGAVG: Running average of size %i will be performed over input coords.\n",
          Nwindow);

  return 0;
}

// RunningAvg::SeparateInit()
void RunningAvg::SeparateInit(int windowIn, int debugIn) {
  isSeparate = true;
  debug = debugIn;
  Nwindow = windowIn;
  // Reserve space for Nwindow frames
  Window.resize( Nwindow );
  // Frame above which averaging will start
  frameThreshold = Nwindow - 1;
  currentWindow = 0;
  windowNatom = 0;
  // For division of frames, cast Nwindow to double
  d_Nwindow = (double)Nwindow;
}

// RunningAvg::setup()
int RunningAvg::setup() {
  // If windowNatom is 0, this is the first setup.
  // If windowNatom is not 0, setup has been called for another parm.
  // Check if the number of atoms has changed. If so the running average
  // will break.
  if ( currentParm->natom != windowNatom ) {
    if (windowNatom!=0) {
      mprintf("Warning: # atoms in parm %s different than previous parm.\n",
              currentParm->parmName);
      mprintf("         Running average will NOT be carried over between parms!\n");
      return 1;
    }
    windowNatom = currentParm->natom;
    // Set up a frame for each window, no masses
    for (int i = 0; i < Nwindow; i++)
      Window[i].SetupFrame( windowNatom, NULL );
    // Setup frame to hold average coords
    avgFrame.SetupFrame( windowNatom, NULL );
    // Zero avgFrame
    avgFrame.ZeroCoords();
    // Set up result
    resultFrame.SetupFrame( windowNatom, NULL );
  } 

  // Print info for this parm
  if (!isSeparate)
    mprintf("\tRunning average set up for %i atoms.\n",windowNatom);
  return 0;  
}

// RunningAvg::action()
int RunningAvg::action() {
  // If frameNum is >= Nwindow, subtract from avgFrame. currentWindow is at
  // the frame that should be subtracted.
  if (frameNum > frameThreshold) { 
    //mprintf("DBG:\tSubtracting Window[%i] from avgFrame.\n",currentWindow);
    avgFrame -= Window[currentWindow];
  }

  // Add current coordinates to avgFrame
  //mprintf("DBG:\tAdding frame %i to avgFrame.\n",frameNum);
  avgFrame += *currentFrame;

  // Store current coordinates in Window
  //mprintf("DBG:\tAssigning frame %i to window %i (%i = %i)\n",frameNum,currentWindow,
  //        Window[currentWindow].natom, currentFrame->natom);
  Window[currentWindow] = *currentFrame;
  ++currentWindow;
  // If currentWindow is out of range, reset
  if (currentWindow==Nwindow) currentWindow=0;

  // If not enough frames to average yet return 3 to indicate further
  // processing should be suppressed.
  if (frameNum < frameThreshold)
    return (int)ACTION_SUPPRESSCOORDOUTPUT;
  // Otherwise there are enough frames to start processing the running average
  else {
    //mprintf("DBG:\tCalculating average for frame %i\n",frameNum); 
    resultFrame.Divide( avgFrame, d_Nwindow );
    // Set frame
    currentFrame = &resultFrame;
  }

  return (int)ACTION_OK;
} 

