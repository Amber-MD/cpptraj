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
  useMass = false;
  maxwindow = -1;
  separateName = NULL;
  ReferenceParm = NULL;
} 

// RmsAvgCorr::init()
/** Expected call: rmsavgcorr [<mask>] [out <filename>] [output <separatename>] 
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
  maxwindow = actionArgs.getKeyInt("stop",-1);
  // Get Masks
  char *rmsmask = actionArgs.getNextMask();
  Mask0.SetMaskString( rmsmask );

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
  if (maxwindow!=-1) mprintf(", max window %i",maxwindow);
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
  if (parmNatom==0) {
    parmNatom = currentParm->natom;
    ReferenceParm = currentParm;
  } else {
    if (parmNatom != currentParm->natom) {
      mprinterr("Error: RmsAvgCorr: Currently cannot be used with multiple topologies.\n");
      return 1;
    }
  }

  // Set up mask
  if ( currentParm->SetupIntegerMask(Mask0, activeReference) ) {
    mprinterr("Error: RmsAvgCorr::setup: Could not set up mask [%s] for parm %s\n",
              Mask0.MaskString(), currentParm->parmName);
    return 1;
  }
  return 0; 
}

// RmsAvgCorr::action()
/** Store input frames for later calculation of RMSDs using running 
  * averages with different window sizes. 
  */
int RmsAvgCorr::action() {
  // Store coordinates to be used in RMS fit. 
  if (ReferenceCoords.AddCoordsByMask(currentFrame->X, &Mask0)) return 1; 
  return 0;
} 

// RmsAvgCorr::print()
/** Calculate the RMSD using running averages of coordinates in 
  * ReferenceCoords with different window sizes. The average RMSD for each 
  * window size will be the "correlation" value.
  */ 
void RmsAvgCorr::print() {
  double avg;
  int window, frame, WindowMax;
  AmberParm *strippedParm;
  float *coord;
  int natom_coord;
  CpptrajFile separateDatafile;

  double U[9], Trans[6];
  bool first;
  Frame refFrame;
  Frame tgtFrame;

  Frame sumFrame;
  int frameThreshold, subtractWindow;
  double d_Nwindow;

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
  // Strip the reference parm to match everything in ReferenceCoords.
  // NOTE: setup() sets ReferenceParm to the first parm read, also ensures 
  //       that this routine is always called with the same # atoms. This 
  //       could cause trouble if very different parms are called that happen 
  //       to have the same # of atoms.
  strippedParm = ReferenceParm->modifyStateByMask(Mask0.Selected,NULL);
  if (strippedParm==NULL) {
    mprinterr("Error: RmsAvgCorr: Could not create topology to match RMS mask.\n");
    return;
  }
  // Set up frames for reference and target atoms for RMSD calc that match the 
  // stripped parm; this also sets masses in case useMass specified.
  refFrame.SetupFrame(strippedParm->natom, strippedParm->mass);
  tgtFrame.SetupFrame(strippedParm->natom, strippedParm->mass);
  // Set up frame for holding sum of coordindates over window frames. 
  // No need for mass. 
  sumFrame.SetupFrame(strippedParm->natom, NULL);

  // Determine max window size to average over
  if (maxwindow==-1)
    WindowMax = ReferenceCoords.Ncoords();
  else {
    WindowMax = maxwindow+1;
    if (WindowMax > ReferenceCoords.Ncoords()) {
      WindowMax = ReferenceCoords.Ncoords();
      mprintf("Warning: RmsAvgCorr: stop (%i) > max # frames (%s), using max.\n",
              maxwindow,WindowMax);
    }
  }

  // Print calc summary
  mprintf("    RMSAVGCORR: Performing RMSD calcs over running avg of coords with window\n");
  mprintf("                sizes ranging from 1 to %i",WindowMax-1);
  if (useMass) mprintf(", mass-weighted");
  mprintf(".\n");
# ifdef _OPENMP
  // Currently DataSet is not thread-safe. Use temp storage.
  double *Ct_openmp = new double[ WindowMax - 1 ];
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    mprintf("                Parallelizing calculation with %i threads.\n",omp_get_num_threads());
}
#endif

  // First value for Ct (window==1) is just the avg RMSD with no running 
  // averaging. Since all non-rms atoms have been stripped, all atoms in
  // ReferenceCoords will be used.
  // Get coords of first frame for use as reference.
  coord = ReferenceCoords.Coord(0,&natom_coord);
  // Sanity check - #atoms must match parm!
  if (natom_coord != strippedParm->natom) {
    mprinterr("Internal Error: # atoms in coord (%i) does not match # atoms in parm (%i).\n",
              natom_coord, strippedParm->natom);
    delete strippedParm;
    return;
  }
  // Set up first frame as reference
  refFrame.SetupFrameFromCoords( coord );
  // Pre-center reference
  refFrame.CenterReference(Trans+3, useMass);
  // Calc initial RMSD
  avg = 0;
  for (frame = 0; frame < ReferenceCoords.Ncoords(); frame++) {
    coord = ReferenceCoords[frame];
    tgtFrame.SetupFrameFromCoords( coord );
    avg += tgtFrame.RMSD_CenteredRef(refFrame, U, Trans, useMass);
  }
  // DEBUG
/*
  sprintf(debugname,"dbgrmsd.dat");
  debug_data = new DataFile(debugname);
  debug_data->AddSet( rmsdaction.rmsd );
  debug_data->Write();
  delete debug_data;
*/
  avg /= ReferenceCoords.Ncoords(); 
  Ct->Add(0, &avg);
  if (separateName!=NULL)
    separateDatafile.IO->Printf("%8i %lf\n",1,avg);

  // LOOP OVER DIFFERENT RUNNING AVG WINDOW SIZES 
# ifdef _OPENMP
#pragma omp parallel private(window, frame, avg, frameThreshold, subtractWindow, d_Nwindow, coord, first, U, Trans) firstprivate(strippedParm,refFrame,tgtFrame,sumFrame)
{
  //mythread = omp_get_thread_num();
#pragma omp for schedule(dynamic)
#endif
  for (window = 2; window < WindowMax; window++ ) {
    // Initialize and set up running average for this window size
    frameThreshold = window - 2;
    subtractWindow = 0;
    d_Nwindow = (double) window;
    sumFrame.ZeroCoords();
    // LOOP OVER FRAMES
    avg = 0;
    first = true;
    for (frame = 0; frame < ReferenceCoords.Ncoords(); frame++) {
      coord = ReferenceCoords[frame];
      tgtFrame.SetupFrameFromCoords( coord );
      // Add current coordinates to sumFrame
      sumFrame += tgtFrame;
      // Do we have enough frames to start calculating a running avg?
      // If so, start calculating RMSDs
      if ( frame > frameThreshold ) {
        tgtFrame.Divide( sumFrame, d_Nwindow );
        // If first, this is the first running-avgd frame, use as reference
        // for RMSD calc for this window size.
        if (first) {
          // Set corods only for speed (everything else is same anyway)
          refFrame.SetFrameCoords( tgtFrame.X );
          refFrame.CenterReference(Trans+3, useMass);
          first = false;
        }
        avg += tgtFrame.RMSD_CenteredRef(refFrame, U, Trans, useMass);
        // Subtract frame at subtractWindow from sumFrame 
        coord = ReferenceCoords[subtractWindow];
        tgtFrame.SetupFrameFromCoords( coord );
        sumFrame -= tgtFrame;
        ++subtractWindow;
      }
    } // END LOOP OVER FRAMES
    // Take average rmsd for this window size. The number of frames for which
    // RMSD was calculated is (Total#frames) - (window size) + 1
    d_Nwindow = (double)ReferenceCoords.Ncoords() - (double)window + 1;
    avg /= d_Nwindow;
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
  for (window = 1; window < WindowMax-1; window++)
    Ct->Add(window, Ct_openmp+window);
  delete[] Ct_openmp;
#endif
  if (separateName!=NULL)
    separateDatafile.CloseFile();
  delete strippedParm;
}

