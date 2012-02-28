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
  useMass = false;
  maxframes = -1;
  separateName = NULL;
} 

// RmsAvgCorr::init()
/** Expected call: rmsavgcorr [<mask>] 
  */
int RmsAvgCorr::init( ) {
  // Get Keywords
  char *outfilename = actionArgs.getKeyString("out",NULL);
# ifdef _OPENMP
  if (actionArgs.hasKey("output")) {
    mprinterr("Error: 'output' keyword not supported in OpenMP version of rmsavgcorr.\n");
    return 1;
  }
# else
  separateName = actionArgs.getKeyString("output",NULL);
# endif
  useMass = actionArgs.hasKey("mass");
  maxframes = actionArgs.getKeyInt("stop",-1);
  // Get Masks
  rmsmask = actionArgs.getNextMask();

  // Set up dataset to hold correlation 
  Ct = DSL->Add(DOUBLE, actionArgs.getNextString(),"RACorr");
  if (Ct==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(outfilename,Ct);

  mprintf("    RMSAVGCORR:");
  if (rmsmask!=NULL)
    mprintf(" RMSD Mask [%s]",rmsmask);
  else
    mprintf(" All atoms");
  if (useMass) mprintf(" (mass-weighted)");
  if (outfilename!=NULL) mprintf(", Output to %s",outfilename);
  if (maxframes!=-1) mprintf(", max frames %i",maxframes);
  mprintf(".\n");
  if (separateName != NULL)
    mprintf("\tSeparate datafile will be written to %s\n",separateName);

  return 0;
}

// RmsAvgCorr::setup()
/** Check that # atoms never changes. Currently this action will not
  * work with multiple parms.
  */
int RmsAvgCorr::setup() {
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
  * averages with different window sizes. 
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
  int window, frame, FrameMax;
  Frame *referenceFrame;
  AmberParm *referenceParm;
  Rmsd rmsdaction;
  RunningAvg runavgaction;
  CpptrajFile separateDatafile;

  // If 'output' specified open up separate datafile that will be written
  // to as correlation is calculated; useful for very long runs.
  int error = 0;
  if (separateName!=NULL) {
    error += separateDatafile.SetupFile(separateName,WRITE,debug);
    error += separateDatafile.OpenFile();
    if (error>0) {
      mprinterr("Error: Could not set up separate data file %s\n",separateName);
      return;
    }
  }
  // DEBUG
/*
  DataFile *debug_data;
  char debugname[32];
*/
  // setup enforces one parm, so this is safe. 
  referenceParm = ReferenceFrames.GetFrameParm(0);

  if (maxframes==-1)
    FrameMax = ReferenceFrames.NumFrames();
  else
    FrameMax = maxframes+1;

  mprintf("    RMSAVGCORR: Performing RMSD calcs over running avg of coords with window\n");
  mprintf("                sizes ranging from 1 to %i",FrameMax-1);
  if (useMass) mprintf(", mass-weighted");
  mprintf(".\n");
# ifdef _OPENMP
  // Currently DataSet is not thread-safe. Use temp storage.
  double *Ct_openmp = new double[ FrameMax - 1 ];
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    mprintf("                Parallelizing calculation with %i threads.\n",omp_get_num_threads());
}
#endif

  // Initialize RMSD action, first value for Ct (window==1) is 
  // just the avg RMSD with no running averaging.
  rmsdaction.SeparateInit(rmsmask, useMass, debug);
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
  if (separateName!=NULL)
    separateDatafile.IO->Printf("%8i %lf\n",1,avg);

  // LOOP OVER DIFFERENT RUNNING AVG WINDOW SIZES 
# ifdef _OPENMP
#pragma omp parallel private(window, frame, referenceFrame,avg,runavgaction,rmsdaction) firstprivate(referenceParm)
{
  //mythread = omp_get_thread_num();
#pragma omp for
#endif
  for (window = 2; window < FrameMax; window++ ) {
    // Initialize running average and rmsd for this window size
    runavgaction.SeparateInit(window,debug);
    rmsdaction.SeparateInit(rmsmask, useMass, debug);
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
    if (separateName!=NULL)
      separateDatafile.IO->Printf("%8i %lf\n",window, avg);
#   endif
  } // END LOOP OVER WINDOWS
#ifdef _OPENMP
} // END pragma omp parallel
  // Put Ct_openmp into Ct dataset
  for (window = 1; window < FrameMax-1; window++)
    Ct->Add(window, Ct_openmp+window);
  delete[] Ct_openmp;
#endif
  if (separateName!=NULL)
    separateDatafile.CloseFile();
}

