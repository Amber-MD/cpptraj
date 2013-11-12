#include <cmath> // exp
#include "Analysis_Rms2d.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Analysis_Rms2d::Analysis_Rms2d() :
  mode_(NORMAL),
  coords_(0),
  nofit_(false),
  useMass_(false),
  RefParm_(0),
  rmsdataset_(0),
  Ct_(0)
{ } 

void Analysis_Rms2d::Help() {
  mprintf("\t[crdset <crd set>] [<name>] [<mask>] [out <filename>]\n");
  mprintf("\t[dme] [mass] [nofit]\n");
  mprintf("\t[reftraj <traj> [parm <parmname> | parmindex <#>] [<refmask>]]\n");
  mprintf("\t[corr <corrfilename>]\n");
  mprintf("\tCalculate RMSD between all frames in <crd set>, or between frames in\n");
  mprintf("\t<crd set> and frames in <traj>.\n");
  mprintf("\n\t<crd set> can be created with the 'createcrd' command.\n");
}

Analysis::RetType Analysis_Rms2d::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)datasetlist->FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: rms2d: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Get keywords
  nofit_ = analyzeArgs.hasKey("nofit");
  useMass_ = analyzeArgs.hasKey("mass");
  std::string outfilename = analyzeArgs.GetStringKey("out");
  if (outfilename.empty()) outfilename = analyzeArgs.GetStringKey("rmsout");
  DataFile* rmsdFile = DFLin->AddDataFile(outfilename, analyzeArgs);
  std::string reftrajname = analyzeArgs.GetStringKey("reftraj");
  if (!reftrajname.empty()) {
    RefParm_ = PFLin->GetParm(analyzeArgs);
    if (RefParm_==0) {
      mprinterr("Error: Rms2d: Could not get parm for reftraj %s.\n",reftrajname.c_str());
      mprinterr("Error:        Ensure parm has been previously loaded.\n");
      return Analysis::ERR;
    }
    mode_ = REFTRAJ; 
  }
  if (analyzeArgs.hasKey("dme")) mode_ = DME;
  // Check for correlation; if so, reftraj not supported
  DataFile* corrfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("corr"), analyzeArgs);
  if (corrfile != 0 && mode_ == REFTRAJ) {
    mprinterr("Error: Rms2d: Keyword 'corr' not supported with 'reftraj'\n");
    return Analysis::ERR;
  }
  // Require an output filename if corr not specified
  if (rmsdFile == 0 && corrfile == 0) {
    mprinterr("Error: Rms2d: No output filename specified.\n");
    Help();
    return Analysis::ERR;
  }
  // Get mask
  TgtMask_.SetMaskString( analyzeArgs.GetMaskNext() );
  // If reference is trajectory, get mask and open traj
  if (mode_ == REFTRAJ) {
    // Get RMS mask string for reference trajectory
    std::string refmaskexpr = analyzeArgs.GetMaskNext();
    if (refmaskexpr.empty())
      refmaskexpr = TgtMask_.MaskExpression();
    RefMask_.SetMaskString( refmaskexpr );
    // Attempt to set up reference trajectory
    if (RefTraj_.SetupTrajRead(reftrajname, analyzeArgs, RefParm_)) {
      mprinterr("Error: Rms2d: Could not set up reftraj %s.\n", reftrajname.c_str());
      return Analysis::ERR;
    }
  } else {
    // Set up DataSet and DataFile for correlate if specified
    if (corrfile != 0 && rmsdataset_ != 0) {
      Ct_ = datasetlist->AddSetAspect( DataSet::DOUBLE, rmsdataset_->Name(), "Corr" );
      if (Ct_ == 0) return Analysis::ERR;
      corrfile->AddSet( Ct_ );
    }
  }
  // Set up output DataSet
  rmsdataset_ = (DataSet_MatrixFlt*)datasetlist->AddSet( DataSet::MATRIX_FLT, 
                                                         analyzeArgs.GetStringNext(), "Rms2d" );
  if (rmsdataset_ == 0) {
    mprinterr("Error: Could not set up DataSet for calculating 2DRMS.\n");
    return Analysis::ERR;
  }
  // Format DataSet 
  rmsdataset_->SetPrecision(8,3);
  // Add to output file
  if (rmsdFile != 0) {
    rmsdFile->AddSet( rmsdataset_ );
    rmsdFile->ProcessArgs("square2d");
  }

  mprintf("    RMS2D: COORDS set [%s], mask [%s]", coords_->Legend().c_str(),
          TgtMask_.MaskString());
  if (mode_ == REFTRAJ) {
      mprintf(", ref traj %s (mask [%s]) %i frames", RefTraj_.TrajFilename().base(),
              RefMask_.MaskString(), RefTraj_.TotalReadFrames());
  }
  if (mode_ == DME) 
    mprintf(", using DME");
  else {
    if (nofit_  ) mprintf(" (no fitting)");
    if (useMass_) mprintf(" (mass-weighted)");
  }
  if (rmsdFile != 0) 
    mprintf(", output to %s",rmsdFile->DataFilename().base());
  mprintf("\n");
  if (corrfile != 0)
    mprintf("           RMSD auto-correlation will be calculated and output to %s\n",
            corrfile->DataFilename().base());

  return Analysis::OK;
}

/** Calculate the RMSD of each frame in ReferenceCoords to each other frame.
  * Since this results in a symmetric matrix use TriangleMatrix to store
  * results.
  */
int Analysis_Rms2d::Calc2drms() 
{
  float R;
  int nref, nframe; 
  int totalref = coords_->Size();
  rmsdataset_->AllocateTriangle( coords_->Size() );
  mprintf("  RMS2D: Calculating RMSDs between each frame (%lu total).\n  ", rmsdataset_->Size());
  // Set up target and reference frames basd on mask
  Frame RefFrame;
  RefFrame.SetupFrameFromMask( TgtMask_, coords_->Top().Atoms() );
  Frame TgtFrame = RefFrame;
  int endref = totalref - 1;
  ParallelProgress progress( endref );
  // LOOP OVER REFERENCE FRAMES
# ifdef _OPENMP
# pragma omp parallel private(nref, nframe, R) firstprivate(TgtFrame, RefFrame, progress)
  {
    progress.SetThread(omp_get_thread_num());
#   pragma omp for schedule(dynamic)
# endif
    for (nref=0; nref < endref; nref++) {
      progress.Update(nref);
      // Get the current reference frame - no box crd
      // TODO: Use coords_->GetFrame instead?
      RefFrame.SetFromCRD( (*coords_)[nref], TgtMask_);
      // Select and pre-center reference atoms (if fitting)
      if (!nofit_)
        RefFrame.CenterOnOrigin(useMass_);
      // LOOP OVER TARGET FRAMES
      for (nframe = nref + 1; nframe < totalref; nframe++) {
        // Get the current target frame
        TgtFrame.SetFromCRD( (*coords_)[nframe], TgtMask_);
        if (nofit_) // Perform no fit RMS calculation
          R = (float)TgtFrame.RMSD_NoFit(RefFrame, useMass_);
        else         // Perform fit RMS calculation
          R = (float)TgtFrame.RMSD_CenteredRef(RefFrame, useMass_);
#       ifdef _OPENMP
        rmsdataset_->SetElement(nframe, nref, R);
#       else
        rmsdataset_->AddElement( R );
#       endif
        // DEBUG
        //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      } // END loop over target frames
    } // END loop over reference frames
# ifdef _OPENMP
  }
# endif
  progress.Finish();
  return 0;
}

// Analysis_Rms2d::CalcDME()
int Analysis_Rms2d::CalcDME() 
{
  int nref, nframe;
  int totalref = coords_->Size();
  rmsdataset_->AllocateTriangle( coords_->Size() );
  mprintf("  RMS2D: Calculating DMEs between each frame (%lu total).\n  ", rmsdataset_->Size());
  // Set up target and reference frames basd on mask
  Frame RefFrame;
  RefFrame.SetupFrameFromMask( TgtMask_, coords_->Top().Atoms() );
  Frame TgtFrame = RefFrame;
  int endref = totalref - 1;
  ParallelProgress progress(endref);
  // LOOP OVER REFERENCE FRAMES
# ifdef _OPENMP
# pragma omp parallel private(nref, nframe) firstprivate(TgtFrame, RefFrame, progress)
  {
    progress.SetThread(omp_get_thread_num());
#   pragma omp for schedule(dynamic)
#   endif
    for (nref=0; nref < endref; nref++) {
      progress.Update(nref);
      // Get the current reference frame - no box crd
      RefFrame.SetFromCRD( (*coords_)[nref], TgtMask_);
      // LOOP OVER TARGET FRAMES
      for (nframe=nref+1; nframe < totalref; nframe++) {
        // Get the current target frame
        TgtFrame.SetFromCRD( (*coords_)[nframe], TgtMask_);
        // Perform DME calc
#       ifdef _OPENMP
        rmsdataset_->SetElement( nframe, nref, TgtFrame.DISTRMSD(RefFrame) );
#       else
        rmsdataset_->AddElement( TgtFrame.DISTRMSD(RefFrame) );
#       endif
        // DEBUG
        //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      } // END loop over target frames
    } // END loop over reference frames
# ifdef _OPENMP
  }
# endif
  progress.Finish();
  return 0;
}

/** Calculate the autocorrelation of the RMSDs. For proper weighting
  * exp[ -RMSD(framei, framei+lag) ] is used. This takes advantage of
  * the fact that 0.0 RMSD essentially means perfect correlation (1.0).
  */
void Analysis_Rms2d::CalcAutoCorr() {
  int N = (int)rmsdataset_->Nrows();
  int lagmax = N;
  // By definition for lag == 0 RMS is 0 for all frames,
  // translates to correlation of 1.
  double ct = 1;
  Ct_->Add(0, &ct); 
  for (int i = 1; i < lagmax; i++) {
    ct = 0;
    int jmax = N - i;
    for (int j = 0; j < jmax; j++) {
      // When i == j GetElement returns 0.0
      double rmsval = rmsdataset_->GetElement(j, j+i);
      ct += exp( -rmsval );
    }
    ct /= jmax;
    Ct_->Add(i, &ct);
  }
}

/** Calc RMSD of every frame in reference traj to every frame in 
  * ReferenceCoords.
  */
int Analysis_Rms2d::CalcRmsToTraj() {
  double R;

  // Set up reference mask for reference parm
  if (RefParm_->SetupIntegerMask(RefMask_)) {
    mprinterr("Error: Could not set up reference mask [%s] for parm %s\n",
              RefMask_.MaskString(), RefParm_->c_str());
    return 1;
  }
  RefMask_.MaskInfo();
  // Ensure # ref atoms == # tgt atoms
  if (RefMask_.Nselected() != TgtMask_.Nselected()) {
    mprinterr("Error: # Selected atoms in ref traj not equal to selected # atoms in\n");
    mprinterr("Error: %s (%i)\n", coords_->Legend().c_str(), TgtMask_.Nselected());
    return 1;
  }
  // Setup reference frame for selected reference atoms
  Frame RefFrame( RefParm_->Atoms() );
  Frame SelectedRef( RefFrame, RefMask_ );
  int totalref = RefTraj_.TotalReadFrames();
  // Setup target from from Coords
  Frame SelectedTgt;
  SelectedTgt.SetupFrameFromMask( TgtMask_, coords_->Top().Atoms() );
  int totaltgt = coords_->Size();
  // Set up output matrix
  int max = totalref * totaltgt;
  mprintf("  RMS2D: Calculating RMSDs between each input frame and each reference\n"); 
  mprintf("         trajectory %s frame (%i total).\n  ",
          RefTraj_.TrajFilename().base(), max);
  rmsdataset_->Allocate2D( totalref, totaltgt );
  if (RefTraj_.BeginTraj(true)) {
    mprinterr("Error: Rms2d: Could not open reference trajectory.\n");
    return 1;
  }

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref; nref++) {
    // Get the current reference frame from trajectory
    RefTraj_.GetNextFrame(RefFrame);
    // Set reference atoms and pre-center if fitting
    SelectedRef.SetCoordinates(RefFrame, RefMask_);
    if (!nofit_)
      SelectedRef.CenterOnOrigin(useMass_);
    // LOOP OVER TARGET FRAMES
    for (int nframe=0; nframe < totaltgt; nframe++) {
      // Get selected atoms of the current target frame
      SelectedTgt.SetFromCRD( (*coords_)[nframe], TgtMask_);
      if (nofit_) {
        // Perform no fit RMS calculation
        R = SelectedTgt.RMSD_NoFit(SelectedRef, useMass_);
      } else {
        // Perform fit RMS calculation
        R = SelectedTgt.RMSD_CenteredRef(SelectedRef, useMass_);
      }
      rmsdataset_->SetElement( nref, nframe, R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
    } // END loop over target frames
  } // END loop over reference frames
  RefTraj_.EndTraj();
  return 0;
}

/** Perform the rms calculation of each frame to each other frame. */
Analysis::RetType Analysis_Rms2d::Analyze() {
  int err = 0;
  // Set up target mask
  if (coords_->Top().SetupIntegerMask( TgtMask_ )) return Analysis::ERR;
  TgtMask_.MaskInfo();
  if (TgtMask_.None()) { 
    mprinterr("Error: No atoms selected for [%s]\n", TgtMask_.MaskString());
    return Analysis::ERR;
  }
  switch (mode_) {
    case NORMAL:
      err = Calc2drms( );
      if (Ct_ != 0) CalcAutoCorr( );
      break;
    case REFTRAJ: 
      err = CalcRmsToTraj( ); 
      break;
    case DME: 
      err = CalcDME( ); 
      if (Ct_ != 0) CalcAutoCorr( );
      break;
  }
  if (err != 0) return Analysis::ERR;
  return Analysis::OK;
}
