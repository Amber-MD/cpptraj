// RmsAvgCorr
#include "Action_RmsAvgCorr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
RmsAvgCorr::RmsAvgCorr() {
  parmNatom = 0;
} 

// DESTRUCTOR
RmsAvgCorr::~RmsAvgCorr() {
}

// RmsAvgCorr::init()
/** Expected call: rmsavgcorr [<mask>] [<refmask>]
  *                [ref <refname> | refindex <ref#> | reference | first ]
  */
int RmsAvgCorr::init( ) {
  std::string refcmd="first";
  std::string mask1="";
  std::string mask2="";
  std::string rmsdargs;
  // Get Keywords
  char *outfilename = actionArgs.getKeyString("out",NULL);
  // Get keywords that should be passed to RMSD action
  if (actionArgs.hasKey("reference")) 
    refcmd="reference";
  else if (actionArgs.hasKey("first"))
    refcmd="first";
  else {
    char *referenceName=actionArgs.getKeyString("ref",NULL);
    char *refindex=actionArgs.getKeyString("refindex",NULL);
    if (referenceName!=NULL) {
      refcmd.assign(referenceName);
      refcmd = "ref " + refcmd;
    } else if (refindex!=NULL) {
      refcmd.assign(refindex);
      refcmd = "refindex " + refcmd;
    }
  }
    
  // Get Masks - will be passed to RMSD action
  char *m1 = actionArgs.getNextMask();
  if (m1!=NULL) mask1.assign(m1);
  char *m2 = actionArgs.getNextMask();
  if (m2!=NULL) mask2.assign(m2);

  // Set up arg list for RMSD action
  rmsdargs = mask1 + " " + mask2 + " " + refcmd;
  rmsdArglist.SetList((char*)rmsdargs.c_str()," ");
  rmsdaction.SetArg(rmsdArglist);

  // Set up dataset to hold correlation 
  Ct = DSL->Add(DOUBLE, actionArgs.getNextString(),"RACorr");
  if (Ct==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(outfilename,Ct);

  mprintf("    RMSAVGCORR:");
  if (outfilename!=NULL) mprintf(" Output to %s",outfilename);
  mprintf("\n\tRmsd Args: [%s]\n",rmsdArglist.ArgLine());
  // Initialize RMSD action, first value for Ct is just the avg RMSD 
  rmsdaction.Init(&localdata, FL, DFL, PFL, debug);

  return 0;
}

// RmsAvgCorr::setup()
int RmsAvgCorr::setup() {
  // Check that # atoms never changes. Currently this action will not
  // work with multiple parms.
  if (parmNatom==0)
    parmNatom = currentParm->natom;
  else {
    if (parmNatom != currentParm->natom) {
      mprinterr("Error: RmsAvgCorr: Currently cannot be used with multiple topologies.\n");
      return 1;
    }
  }
 
  // Setup RMSD action 
  return (rmsdaction.Setup(&currentParm));
}

// RmsAvgCorr::action()
/** Store input frames for later calculation of RMSDs using running 
  * averages with different window sizes. Also calculate the initial
  * RMSD since this is correlation with running average window size
  * == 1 (i.e. no averaging).
  */
int RmsAvgCorr::action() {
  // Store frame
  Frame *fCopy = currentFrame->FrameCopy();
  ReferenceFrames.AddFrame(fCopy,currentParm);
  // Calc initial RMSD
  rmsdaction.DoAction(&currentFrame, frameNum);
  return 0;
} 

// RmsAvgCorr::print()
/** Calculate the RMSD using running averages of frames in ReferenceFrames
  * with different window sizes. The average RMSD for each window size
  * will be the "correlation" value.
  */ 
void RmsAvgCorr::print() {
  double ZeroData = 0; // Used to blank frames in RMSD data set < window size
  double avg;
  RunningAvg runavgaction;

  // First data point (window==1) is just the normal rmsd, 
  // already calcd in action.
  avg = rmsdaction.rmsd->Avg(NULL);
  Ct->Add(0, &avg);
  // setup enforces one parm, so this is safe. 
  currentParm = ReferenceFrames.GetFrameParm(0);
  // LOOP OVER DIFFERENT RUNNING AVG WINDOW SIZES 
  for (int window = 2; window < ReferenceFrames.NumFrames(); window++ ) {
    // Set up arg list for running average
    runavgaction.ArgumentList("window %i",window); 
    // Init running average
    runavgaction.Init(&localdata, FL, DFL, PFL, debug);
    // Setup running average 
    runavgaction.Setup(&currentParm);
    // Re-initialize Rmsd action
    rmsdaction.SetArg( rmsdArglist );
    rmsdaction.Init(&localdata, FL, DFL, PFL, debug);
    // Setup Rmsd action
    rmsdaction.Setup(&currentParm);
    // LOOP OVER FRAMES
    for (int frame = 0; frame < ReferenceFrames.NumFrames(); frame++) {
      currentFrame = ReferenceFrames.GetFrame( frame );
      // If running average indicates its good to go, perform RMSD.
      // Otherwise set corresponding point in RMSD dataset to zero.
      if (runavgaction.DoAction(&currentFrame, frame) == ACTION_OK) {
        rmsdaction.DoAction(&currentFrame, frame);
      } else {
        rmsdaction.rmsd->Add(frame, &ZeroData);
      }
    } // END LOOP OVER FRAMES
    // Take average rmsd for this window size
    avg = rmsdaction.rmsd->Avg(NULL);
    Ct->Add(window-1, &avg);
  } // END LOOP OVER WINDOWS
}

