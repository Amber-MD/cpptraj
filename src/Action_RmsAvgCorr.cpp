// RmsAvgCorr
#include "Action_RmsAvgCorr.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
// DEBUG
//#include <cstdio> //sprintf

// CONSTRUCTOR
RmsAvgCorr::RmsAvgCorr() {
  parmNatom = 0;
  rmsmask = NULL;
} 

// DESTRUCTOR
RmsAvgCorr::~RmsAvgCorr() {
}

// RmsAvgCorr::init()
/** Expected call: rmsavgcorr [<mask>] 
  */
int RmsAvgCorr::init( ) {
  // Get Keywords
  char *outfilename = actionArgs.getKeyString("out",NULL);
  // Get Masks
  rmsmask = actionArgs.getNextMask();

  // Set up dataset to hold correlation 
  Ct = DSL->Add(DOUBLE, actionArgs.getNextString(),"RACorr");
  if (Ct==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(outfilename,Ct);

  mprintf("    RMSAVGCORR:");
  if (rmsmask!=NULL)
    mprintf(" Mask [%s]",rmsmask);
  else
    mprintf(" All atoms");
  if (outfilename!=NULL) mprintf(", Output to %s",outfilename);

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
  return 0; 
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
  int window, frame;
  Frame *referenceFrame;
  AmberParm *referenceParm;
  Rmsd rmsdaction;
  RunningAvg runavgaction;
  // DEBUG
/*
  DataFile *debug_data;
  char debugname[32];
*/
  // setup enforces one parm, so this is safe. 
  referenceParm = ReferenceFrames.GetFrameParm(0);

  mprintf("    RMSAVGCORR: Performing RMSD calcs over running avg of coords with window\n");
  mprintf("                sizes ranging from 1 to %i\n",ReferenceFrames.NumFrames()-1);
# ifdef _OPENMP
  // Currently DataSet is not thread-safe. Use temp storage.
  double *Ct_openmp = new double[ ReferenceFrames.NumFrames() - 1 ];
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    mprintf("                Parallelizing calculation with %i threads.\n",omp_get_num_threads());
}
#endif

  // Initialize RMSD action, first value for Ct (window==1) is 
  // just the avg RMSD with no running averaging.
  rmsdaction.SeparateInit(rmsmask, debug);
   // Setup RMSD action 
  rmsdaction.Setup(&referenceParm);
  // Calc initial RMSD
  Frame *tempFrame = new Frame(); 
  for (frame = 0; frame < ReferenceFrames.NumFrames(); frame++) {
    referenceFrame = ReferenceFrames.GetFrame( frame );
    // Protect referenceFrame from being rotated
    *tempFrame = *referenceFrame;
    rmsdaction.DoAction(&tempFrame, frame);
  }
  delete tempFrame;
  // DEBUG
/*
  sprintf(debugname,"dbgrmsd.dat");
  debug_data = new DataFile(debugname);
  debug_data->AddSet( rmsdaction.rmsd );
  debug_data->Write();
  delete debug_data;
*/
  avg = rmsdaction.rmsd->Avg(NULL);
  Ct->Add(0, &avg);

  // LOOP OVER DIFFERENT RUNNING AVG WINDOW SIZES 
# ifdef _OPENMP
#pragma omp parallel private(window, frame, referenceFrame,avg,runavgaction,rmsdaction) firstprivate(referenceParm)
{
  //mythread = omp_get_thread_num();
#pragma omp for
#endif
  for (window = 2; window < ReferenceFrames.NumFrames(); window++ ) {
    // Initialize running average and rmsd for this window size
    runavgaction.SeparateInit(window,debug);
    rmsdaction.SeparateInit(rmsmask, debug);
    // Setup running average and rmsd for this window size
    runavgaction.Setup(&referenceParm);
    rmsdaction.Setup(&referenceParm);
    // LOOP OVER FRAMES
    for (frame = 0; frame < ReferenceFrames.NumFrames(); frame++) {
      referenceFrame = ReferenceFrames.GetFrame( frame );
      // If running average indicates its good to go, perform RMSD.
      // Otherwise set corresponding point in RMSD dataset to zero.
      if (runavgaction.DoAction(&referenceFrame, frame) == ACTION_OK) 
        rmsdaction.DoAction(&referenceFrame, frame);
      else 
        rmsdaction.rmsd->Add(frame, &ZeroData);
    } // END LOOP OVER FRAMES
    // Take average rmsd for this window size
    avg = rmsdaction.rmsd->Avg(NULL);
#   ifdef _OPENMP
    Ct_openmp[window-1] = avg;
#   else 
    Ct->Add(window-1, &avg);
#   endif
  } // END LOOP OVER WINDOWS
#ifdef _OPENMP
} // END pragma omp parallel
  // Put Ct_openmp into Ct dataset
  for (window = 1; window < ReferenceFrames.NumFrames()-1; window++)
    Ct->Add(window, Ct_openmp+window);
  delete[] Ct_openmp;
#endif
}

