// Action_RmsAvgCorr
#include "Action_RmsAvgCorr.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_RmsAvgCorr::Action_RmsAvgCorr() :
  ReferenceParm_(NULL),
  Ct_(NULL),
  parmNatom_(0),
  maxwindow_(-1),
  useMass_(false)
{ } 

void Action_RmsAvgCorr::Help() {
  mprintf("rmsavgcorr [<mask>] [out <filename>] [output <separatename>]\n");
}

// Action_RmsAvgCorr::init()
Action::RetType Action_RmsAvgCorr::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  ArgList::ConstArg outfilename = actionArgs.getKeyString("out");
# ifdef _OPENMP
  if (actionArgs.hasKey("output")) {
    mprinterr("Error: 'output' keyword not supported in OpenMP version of rmsavgcorr.\n");
    return 1;
  }
# else
  separateName_ = actionArgs.GetStringKey("output");
# endif
  useMass_ = actionArgs.hasKey("mass");
  maxwindow_ = actionArgs.getKeyInt("stop",-1);
  // Get Masks
  ArgList::ConstArg rmsmask = actionArgs.getNextMask();
  Mask0_.SetMaskString( rmsmask );

  // Set up dataset to hold correlation 
  Ct_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"RACorr");
  if (Ct_==NULL) return Action::ERR;
  // Add dataset to data file list
  DFL->Add(outfilename,Ct_);

  mprintf("    RMSAVGCORR:");
  if (rmsmask!=NULL)
    mprintf(" RMSD Mask [%s]",rmsmask);
  else
    mprintf(" All atoms");
  if (useMass_) mprintf(" (mass-weighted)");
  if (outfilename!=NULL) mprintf(", Output to %s",outfilename);
  if (maxwindow_!=-1) mprintf(", max window %i",maxwindow_);
  mprintf(".\n");
  if (!separateName_.empty())
    mprintf("\tSeparate datafile will be written to %s\n",separateName_.c_str());

  return Action::OK;
}

// Action_RmsAvgCorr::setup()
/** Check that # atoms never changes. Currently this action will not
  * work with multiple parms.
  */
Action::RetType Action_RmsAvgCorr::Setup(Topology* currentParm, Topology** parmAddress) {
  if (parmNatom_==0) {
    parmNatom_ = currentParm->Natom();
    ReferenceParm_ = currentParm;
  } else {
    if (parmNatom_ != currentParm->Natom()) {
      mprinterr("Error: RmsAvgCorr: Currently cannot be used with multiple topologies.\n");
      return Action::ERR;
    }
  }

  // Set up mask
  if ( currentParm->SetupIntegerMask(Mask0_) ) {
    mprinterr("Error: RmsAvgCorr::setup: Could not set up mask [%s] for parm %s\n",
              Mask0_.MaskString(), currentParm->c_str());
    return Action::ERR;
  }
  return Action::OK; 
}

// Action_RmsAvgCorr::action()
/** Store input frames for later calculation of RMSDs using running 
  * averages with different window sizes. 
  */
Action::RetType Action_RmsAvgCorr::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Store coordinates to be used in RMS fit. 
  if (ReferenceCoords_.AddFrameByMask(*currentFrame, Mask0_))
    return Action::ERR;
  return Action::OK;
} 

// Action_RmsAvgCorr::print()
/** Calculate the RMSD using running averages of coordinates in 
  * ReferenceCoords with different window sizes. The average RMSD for each 
  * window size will be the "correlation" value.
  */ 
void Action_RmsAvgCorr::Print() {
  double avg;
  int window, frame, WindowMax;
  Topology *strippedParm;
  CpptrajFile separateDatafile;

  double U[9], Trans[6];
  bool first;

  int frameThreshold, subtractWindow;
  double d_Nwindow;

  // If 'output' specified open up separate datafile that will be written
  // to as correlation is calculated; useful for very long runs.
  if (!separateName_.empty()) {
    if (separateDatafile.OpenWrite(separateName_)) {
      mprinterr("Error: Could not set up separate data file %s\n",separateName_.c_str());
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
  strippedParm = ReferenceParm_->modifyStateByMask(Mask0_);
  if (strippedParm==NULL) {
    mprinterr("Error: RmsAvgCorr: Could not create topology to match RMS mask.\n");
    return;
  }
  // Ensure that the max # atoms in ReferenceCoords does not exceed # atoms
  // in stripped topology. Sanity check only; setup() ensures this should
  // not happen.
  if (ReferenceCoords_.MaxNatom() > strippedParm->Natom()) {
    mprinterr("Internal Error: Max # atoms in ref coords (%i) != stripped # atoms (%i).\n",
              ReferenceCoords_.MaxNatom(), strippedParm->Natom());
    delete strippedParm;
    return;
  }

  // Set up frames for reference and target atoms for RMSD calc that match the 
  // stripped parm; this also sets masses in case useMass specified.
  Frame refFrame(strippedParm->Atoms());
  Frame tgtFrame(strippedParm->Atoms());
  // Set up frame for holding sum of coordindates over window frames. 
  // No need for mass. 
  Frame sumFrame(strippedParm->Natom());

  // Determine max window size to average over
  if (maxwindow_==-1)
    WindowMax = ReferenceCoords_.Ncoords();
  else {
    WindowMax = maxwindow_+1;
    if (WindowMax > ReferenceCoords_.Ncoords()) {
      WindowMax = ReferenceCoords_.Ncoords();
      mprintf("Warning: RmsAvgCorr: stop (%i) > max # frames (%s), using max.\n",
              maxwindow_,WindowMax);
    }
  }

  // Print calc summary
  mprintf("    RMSAVGCORR: Performing RMSD calcs over running avg of coords with window\n");
  mprintf("                sizes ranging from 1 to %i",WindowMax-1);
  if (useMass_) mprintf(", mass-weighted");
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
  refFrame = ReferenceCoords_[0];
  // Pre-center reference
  refFrame.CenterReference(Trans+3, useMass_);
  // Calc initial RMSD
  avg = 0;
  for (frame = 0; frame < ReferenceCoords_.Ncoords(); frame++) {
    tgtFrame = ReferenceCoords_[frame];
    avg += tgtFrame.RMSD_CenteredRef(refFrame, U, Trans, useMass_);
  }
  // DEBUG
/*
  sprintf(debugname,"dbgrmsd.dat");
  debug_data = new DataFile(debugname);
  debug_data->AddSet( rmsdaction.rmsd );
  debug_data->Write();
  delete debug_data;
*/
  avg /= ReferenceCoords_.Ncoords(); 
  Ct_->Add(0, &avg);
  if (!separateName_.empty())
    separateDatafile.Printf("%8i %lf\n",1,avg);

  // LOOP OVER DIFFERENT RUNNING AVG WINDOW SIZES 
# ifdef _OPENMP
#pragma omp parallel private(window, frame, avg, frameThreshold, subtractWindow, d_Nwindow, first, U, Trans) firstprivate(strippedParm,refFrame,tgtFrame,sumFrame)
{
  //mythread = omp_get_thread_num();
#pragma omp for schedule(dynamic)
#endif
  for (window = 2; window < WindowMax; window++ ) {
    // Initialize and set up running average for this window size
    frameThreshold = window - 2;
    // TODO: Make subtractWindow a const iterator to CoordList
    subtractWindow = 0;
    d_Nwindow = (double) window;
    sumFrame.ZeroCoords();
    // LOOP OVER FRAMES
    avg = 0;
    first = true;
    for (frame = 0; frame < ReferenceCoords_.Ncoords(); frame++) {
      tgtFrame = ReferenceCoords_[frame];
      // Add current coordinates to sumFrame
      sumFrame += tgtFrame;
      // Do we have enough frames to start calculating a running avg?
      // If so, start calculating RMSDs
      if ( frame > frameThreshold ) {
        tgtFrame.Divide( sumFrame, d_Nwindow );
        // If first, this is the first running-avgd frame, use as reference
        // for RMSD calc for this window size.
        if (first) {
          // Set coords only for speed (everything else is same anyway)
          refFrame.SetCoordinates( tgtFrame );
          // Pre-center reference
          refFrame.CenterReference(Trans+3, useMass_);
          first = false;
        }
        avg += tgtFrame.RMSD_CenteredRef(refFrame, U, Trans, useMass_);
        // Subtract frame at subtractWindow from sumFrame 
        tgtFrame = ReferenceCoords_[subtractWindow];
        sumFrame -= tgtFrame;
        ++subtractWindow;
      }
    } // END LOOP OVER FRAMES
    // Take average rmsd for this window size. The number of frames for which
    // RMSD was calculated is (Total#frames) - (window size) + 1
    d_Nwindow = (double)ReferenceCoords_.Ncoords() - (double)window + 1;
    avg /= d_Nwindow;
#   ifdef _OPENMP
    Ct_openmp[window-1] = avg;
#   else 
    Ct_->Add(window-1, &avg);
    if (!separateName_.empty())
      separateDatafile.Printf("%8i %lf\n",window, avg);
#   endif
  } // END LOOP OVER WINDOWS
#ifdef _OPENMP
} // END pragma omp parallel
  // Put Ct_openmp into Ct dataset
  for (window = 1; window < WindowMax-1; window++)
    Ct_->Add(window, Ct_openmp+window);
  delete[] Ct_openmp;
#endif
  if (!separateName_.empty())
    separateDatafile.CloseFile();
  delete strippedParm;
}

