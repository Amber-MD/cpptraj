#include <cmath> // exp
#include "Analysis_Rms2d.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "Matrix_2D.h"

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
  mprintf("rms2d <crd set> [<mask>] rmsout <filename> [nofit] [mass] [dme]\n");
  mprintf("      [reftraj <traj> [parm <parmname> | parmindex <#>] [<refmask>]]\n");
  mprintf("      [corr <corrfilename>]\n");
  mprintf("\tCalculate RMSD between all frames in <crd set>, or between frames in\n");
  mprintf("<crd set> and frames in <traj>.\n");
  mprintf("<crd set> can be created with the 'createcrd' command.\n");
}

Analysis::RetType Analysis_Rms2d::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                                             TopologyList* PFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: rms2d: Specify crd set name.\n");
    Help();
    return Analysis::ERR;
  }
  coords_ = (DataSet_Coords*)datasetlist->FindSetOfType( setname, DataSet::COORDS );
  if (coords_ == 0) {
    mprinterr("Error: rms2d: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Get keywords
  nofit_ = analyzeArgs.hasKey("nofit");
  useMass_ = analyzeArgs.hasKey("mass"); 
  rmsdFile_ = analyzeArgs.GetStringKey("rmsout");
  std::string reftrajname = analyzeArgs.GetStringKey("reftraj");
  if (!reftrajname.empty()) {
    RefParm_ = PFLin->GetParm(analyzeArgs);
    if (RefParm_==0) {
      mprinterr("Error: Rms2d: Could not get parm for reftraj %s.\n",reftrajname.c_str());
      return Analysis::ERR;
    }
    mode_ = REFTRAJ; 
  }
  if (analyzeArgs.hasKey("dme")) mode_ = DME;
  // Check for correlation; if so, reftraj not supported
  corrfilename_ = analyzeArgs.GetStringKey("corr");
  if (!corrfilename_.empty() && mode_ == REFTRAJ) {
    mprinterr("Error: Rms2d: Keyword 'corr' not supported with 'reftraj'\n");
    return Analysis::ERR;
  }
  // Require an output filename if corr not specified
  if (rmsdFile_.empty() && corrfilename_.empty()) {
    mprinterr("Error: Rms2d: No output filename specified.\n");
    Help();
    return Analysis::ERR;
  }
  // Get mask
  maskexpr_ =  analyzeArgs.GetMaskNext();
  // If reference is trajectory, get mask and open traj
  if (mode_ == REFTRAJ) {
    // Get RMS mask string for reference trajectory
    std::string refmaskexpr = analyzeArgs.GetMaskNext();
    if (refmaskexpr.empty()) refmaskexpr = maskexpr_;
    RefMask_.SetMaskString( refmaskexpr );
    // Attempt to set up reference trajectory
    if (RefTraj_.SetupTrajRead(reftrajname, 0, RefParm_)) {
      mprinterr("Error: Rms2d: Could not set up reftraj %s.\n", reftrajname.c_str());
      return Analysis::ERR;
    }
    // Set up DataSet for RefTraj, full 2d matrix
    rmsdataset_ = datasetlist->AddSet( DataSet::MATRIX2D, analyzeArgs.GetStringNext(), "Rms2d" );
  } else {
    // Set up DataSet for normal rms2d, half matrix
    rmsdataset_ = datasetlist->AddSet( DataSet::TRIMATRIX, analyzeArgs.GetStringNext(), "Rms2d" );
    // Set up DataSet and DataFile for correlate
    if (!corrfilename_.empty() && rmsdataset_ != 0) {
      Ct_ = datasetlist->AddSetAspect( DataSet::DOUBLE, rmsdataset_->Name(), "Corr" );
      if (Ct_ == 0) return Analysis::ERR;
    }
  }
  if (rmsdataset_ == 0) {
    mprinterr("Error: Could not set up DataSet for calculating 2DRMS.\n");
    return Analysis::ERR;
  }
  // Format DataSet 
  rmsdataset_->SetPrecision(8,3);

  mprintf("    RMS2D: COORDS set [%s]", coords_->Legend().c_str());
  if (!maskexpr_.empty())
    mprintf(", mask [%s]", maskexpr_.c_str());
  else
    mprintf(", all atoms");
  switch (mode_) {
    case REFTRAJ:
      mprintf(", ref traj %s (mask [%s]) %i frames", RefTraj_.FullTrajStr(),
              RefMask_.MaskString(), RefTraj_.TotalReadFrames());
      break;
    case DME: mprintf(", using DME"); break;
    case NORMAL: // RMSD
      if (nofit_)
        mprintf(" (no fitting)");
      if (useMass_)
        mprintf(" (mass-weighted)");
  }
  if (!rmsdFile_.empty()) 
    mprintf(" output to %s",rmsdFile_.c_str());
  mprintf("\n");
  if (!corrfilename_.empty())
    mprintf("           RMSD auto-correlation will be calculated and output to %s\n",
            corrfilename_.c_str());

  return Analysis::OK;
}

/** Calculate the RMSD of each frame in ReferenceCoords to each other frame.
  * Since this results in a symmetric matrix use TriangleMatrix to store
  * results.
  */
int Analysis_Rms2d::Calc2drms(DataSet_Coords& coordsIn, TriangleMatrix& Distances,
                              bool nofitIn, bool useMassIn, std::string const& maskexpr) 
{
  float R;
 
  int totalref = coordsIn.Size();
  Distances.Setup( totalref );
  int max = Distances.Nelements();
  mprintf("  RMS2D: Calculating RMSDs between each frame (%i total).\n  ",max);

  // Set up progress Bar
  ProgressBar progress(totalref - 1);
  // Set up mask
  AtomMask tgtmask( maskexpr );
  if (coordsIn.Top().SetupIntegerMask( tgtmask )) return 1;
  tgtmask.MaskInfo();
  if (tgtmask.None()) return 1;
  // Set up target and reference frames basd on mask
  Frame RefFrame;
  RefFrame.SetupFrameFromMask( tgtmask, coordsIn.Top().Atoms() );
  Frame TgtFrame = RefFrame;

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref - 1; nref++) {
    progress.Update(nref);
    // Get the current reference frame - no box crd
    RefFrame.SetFromCRD(coordsIn[nref], 0, tgtmask);
    // Select and pre-center reference atoms (if fitting)
    if (!nofitIn)
      RefFrame.CenterReference(useMassIn);
  
    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      // Get the current target frame
      TgtFrame.SetFromCRD(coordsIn[nframe], 0, tgtmask);
      if (nofitIn) {
        // Perform no fit RMS calculation
        R = (float) TgtFrame.RMSD_NoFit(RefFrame, useMassIn);
      } else {
        // Perform fit RMS calculation
        R = (float) TgtFrame.RMSD_CenteredRef(RefFrame, useMassIn);
      }
      Distances.AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
    } // END loop over target frames
  } // END loop over reference frames

  return 0;
}

// Analysis_Rms2d::CalcDME()
int Analysis_Rms2d::CalcDME(DataSet_Coords& coordsIn, TriangleMatrix& Distances,
                            std::string const& maskexpr) 
{
  int totalref = coordsIn.Size();
  Distances.Setup( totalref );
  int max = Distances.Nelements();
  mprintf("  RMS2D: Calculating DMEs between each frame (%i total).\n  ",max);

  // Set up progress Bar
  ProgressBar progress(totalref - 1);
  // Set up mask
  AtomMask tgtmask( maskexpr );
  if (coordsIn.Top().SetupIntegerMask( tgtmask )) return 1;
  tgtmask.MaskInfo();
  if (tgtmask.None()) return 1;
  // Set up target and reference frames basd on mask
  Frame RefFrame;
  RefFrame.SetupFrameFromMask( tgtmask, coordsIn.Top().Atoms() );
  Frame TgtFrame = RefFrame;

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref - 1; nref++) {
    progress.Update(nref);
    // Get the current reference frame - no box crd
    RefFrame.SetFromCRD(coordsIn[nref], 0, tgtmask);
    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      // Get the current target frame
      TgtFrame.SetFromCRD(coordsIn[nframe], 0, tgtmask);
      // Perform DME calc
      Distances.AddElement( (float)TgtFrame.DISTRMSD(RefFrame) );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
    } // END loop over target frames
  } // END loop over reference frames

  return 0;
}

/** Calculate the autocorrelation of the RMSDs. For proper weighting
  * exp[ -RMSD(framei, framei+lag) ] is used. This takes advantage of
  * the fact that 0.0 RMSD essentially means perfect correlation (1.0).
  */
void Analysis_Rms2d::CalcAutoCorr(TriangleMatrix& Distances) {
  int N = Distances.Nrows();
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
      double rmsval = Distances.GetElement(j, j+i);
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
  float R;

  Matrix_2D* rmsdata = (Matrix_2D*)rmsdataset_;
  // Set up mask
  AtomMask TgtMask( maskexpr_ );
  if (coords_->Top().SetupIntegerMask( TgtMask )) return 1;
  TgtMask.MaskInfo();
  // Set up reference mask for reference parm
  if (RefParm_->SetupIntegerMask(RefMask_)) {
    mprinterr("Error: Could not set up reference mask [%s] for parm %s\n",
              RefMask_.MaskString(), RefParm_->c_str());
    return 1;
  }
  RefMask_.MaskInfo();
  // Ensure # ref atoms == # tgt atoms
  if (RefMask_.Nselected() != TgtMask.Nselected()) {
    mprinterr("Error: # Selected atoms in ref traj not equal to selected # atoms in\n");
    mprinterr("Error: %s (%i)\n", coords_->Legend().c_str(), TgtMask.Nselected());
    return 1;
  }
  // Setup reference frame for selected reference atoms
  Frame RefFrame( RefParm_->Atoms() );
  Frame SelectedRef( RefFrame, RefMask_ );
  int totalref = RefTraj_.TotalReadFrames();
  // Setup target from from Coords
  Frame SelectedTgt;
  SelectedTgt.SetupFrameFromMask( TgtMask, coords_->Top().Atoms() );
  int totaltgt = coords_->Size();
  // Set up output matrix
  int max = totalref * totaltgt;
  mprintf("  RMS2D: Calculating RMSDs between each input frame and each reference\n"); 
  mprintf("         trajectory %s frame (%i total).\n  ",
          RefTraj_.BaseTrajStr(), max);
  rmsdata->Setup( totalref, totaltgt );
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
      SelectedRef.CenterReference(useMass_);
    // LOOP OVER TARGET FRAMES
    for (int nframe=0; nframe < totaltgt; nframe++) {
      // Get selected atoms of the current target frame
      SelectedTgt.SetFromCRD( (*coords_)[nframe], 0, TgtMask);
      if (nofit_) {
        // Perform no fit RMS calculation
        R = (float) SelectedTgt.RMSD_NoFit(SelectedRef, useMass_);
      } else {
        // Perform fit RMS calculation
        R = (float) SelectedTgt.RMSD_CenteredRef(SelectedRef, useMass_);
      }
      rmsdata->AddElement( R );
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
  switch (mode_) {
    case NORMAL:
      err = Calc2drms( *coords_, *((TriangleMatrix*)rmsdataset_), nofit_, useMass_, maskexpr_ );
      if (Ct_ != 0) CalcAutoCorr( *((TriangleMatrix*)rmsdataset_) );
      break;
    case REFTRAJ: 
      err = CalcRmsToTraj(); 
      break;
    case DME: 
      err = CalcDME( *coords_, *((TriangleMatrix*)rmsdataset_), maskexpr_ ); 
      if (Ct_ != 0) CalcAutoCorr( *((TriangleMatrix*)rmsdataset_) );
      break;
  }
  if (err != 0) return Analysis::ERR;
  return Analysis::OK;
}

void Analysis_Rms2d::Print(DataFileList* DFL) {
  if (!corrfilename_.empty())
    DFL->AddSetToFile(corrfilename_, Ct_);
  if (!rmsdFile_.empty()) {
    DataFile* rmsdatafile = DFL->AddSetToFile( rmsdFile_, rmsdataset_ );
    rmsdatafile->ProcessArgs("square2d");
  }
}
