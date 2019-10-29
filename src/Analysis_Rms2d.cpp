#include <cmath> // exp
#include "Analysis_Rms2d.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "DataSet_Coords_TRJ.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

// CONSTRUCTOR
Analysis_Rms2d::Analysis_Rms2d() :
  mode_(RMS_FIT),
  TgtTraj_(0),
  useReferenceTraj_(false),
  useMass_(false),
  RefTraj_(0),
  RefParm_(0),
  rmsdataset_(0),
  Ct_(0)
{ } 

void Analysis_Rms2d::Help() const {
  mprintf("\t[crdset <crd set>] [<name>] [<mask>] [out <filename>]\n"
          "\t[dme | nofit | srmsd] [mass]\n"
          "\t[reftraj <traj> [parm <parmname> | parmindex <#>] [<refmask>]]\n"
          "\t[corr <corrfilename>]\n"
          "  Calculate RMSD between all frames in <crd set>, or between frames in\n"
          "  <crd set> and frames in <traj>.\n"
          "  <crd set> can be created with the 'createcrd' command.\n");
}

const char* Analysis_Rms2d::ModeStrings_[] = {
  "RMSD (fit)", "RMSD (no fitting)", "DME", "symmetry-corrected RMSD"
};

// Analysis_Rms2d::Setup()
Analysis::RetType Analysis_Rms2d::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (TgtTraj_ == 0) {
    mprinterr("Error: rms2d: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Get keywords
  if (analyzeArgs.hasKey("nofit"))
    mode_ = RMS_NOFIT;
  else if (analyzeArgs.hasKey("dme"))
    mode_ = DME;
  else if (analyzeArgs.hasKey("srmsd"))
    mode_ = SRMSD;
  else
    mode_ = RMS_FIT;
  useMass_ = analyzeArgs.hasKey("mass");
  std::string outfilename = analyzeArgs.GetStringKey("out");
  if (outfilename.empty()) outfilename = analyzeArgs.GetStringKey("rmsout"); // DEPRECATED
  DataFile* rmsdFile = setup.DFL().AddDataFile(outfilename, "square2d", analyzeArgs);
  std::string reftrajname = analyzeArgs.GetStringKey("reftraj");
  if (!reftrajname.empty()) {
    RefParm_ = setup.DSL().GetTopology(analyzeArgs); // TODO Use coords set
    useReferenceTraj_ = true;
  } else
    useReferenceTraj_ = false;
  // Check for correlation. 
  DataFile* corrfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("corr"), analyzeArgs);
  // Get target mask.
  if (TgtMask_.SetMaskString( analyzeArgs.GetMaskNext() )) return Analysis::ERR;
  // Get reference mask. 
  std::string refmaskexpr = analyzeArgs.GetMaskNext();
  if (refmaskexpr.empty())
    refmaskexpr = TgtMask_.MaskExpression();
  if (RefMask_.SetMaskString( refmaskexpr )) return Analysis::ERR;
  // When 'corr', reftraj or coords with different mask not supported
  if (corrfile != 0 ) {
    if (useReferenceTraj_) {
      mprinterr("Error: Keyword 'corr' not supported with 'reftraj'\n");
      return Analysis::ERR;
    } else {
      if (TgtMask_.MaskExpression() != RefMask_.MaskExpression()) {
        mprinterr("Error: Keyword 'corr' not supported when masks differ if not using 'reftraj'\n");
        return Analysis::ERR;
      }
    }
  }
  // If reference is trajectory, open traj
  if (useReferenceTraj_) {
    // Find out if reference traj is already present
    RefTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( reftrajname );
    if (RefTraj_ == 0) {
      // Reference traj not yet present. Load it; requires parm.
      if (RefParm_==0) {
        mprinterr("Error: Could not get parm for reftraj %s.\n"
                  "Error:   Ensure parm has been previously loaded.\n",reftrajname.c_str());
        return Analysis::ERR;
      }
      DataSet_Coords_TRJ* DCT = (DataSet_Coords_TRJ*)
                                setup.DSL().AddSet(DataSet::TRAJ, reftrajname, "RmsRefTraj");
      if (DCT == 0) return Analysis::ERR;
      if (DCT->AddSingleTrajin( reftrajname, analyzeArgs, RefParm_ )) return Analysis::ERR;
      RefTraj_ = (DataSet_Coords*)DCT;
    } else
      RefParm_ = RefTraj_->TopPtr();
  }
  // Set up output DataSet
  rmsdataset_ = (DataSet_MatrixFlt*)setup.DSL().AddSet( DataSet::MATRIX_FLT, 
                                                         analyzeArgs.GetStringNext(), "Rms2d" );
  if (rmsdataset_ == 0) {
    mprinterr("Error: Could not set up DataSet for calculating 2DRMS.\n");
    return Analysis::ERR;
  }
  // Format DataSet 
  rmsdataset_->SetupFormat().SetFormatWidthPrecision(8,3);
  // Add to output file
  if (rmsdFile != 0)
    rmsdFile->AddDataSet( rmsdataset_ );
  // Set up DataSet for corr if specified
  if (corrfile != 0) {
    Ct_ = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(rmsdataset_->Meta().Name(), "Corr") );
    if (Ct_ == 0) return Analysis::ERR;
    corrfile->AddDataSet( Ct_ );
  }

  mprintf("    RMS2D: COORDS set [%s], mask [%s]", TgtTraj_->legend(),
          TgtMask_.MaskString());
  if ( TgtMask_.MaskExpression() != RefMask_.MaskExpression() )
    mprintf(" ref mask [%s]", RefMask_.MaskString());
  mprintf(", using %s", ModeStrings_[mode_]);
  if (useMass_) mprintf(", mass-weighted");
  mprintf("\n");
  if (useReferenceTraj_)
    mprintf("\tReference trajectory '%s', %zu frames\n",
            RefTraj_->legend(), RefTraj_->Size());
  if (rmsdFile != 0) 
    mprintf("\tOutput to '%s'\n",rmsdFile->DataFilename().full());
  if (corrfile != 0)
    mprintf("\tRMSD auto-correlation will be calculated and output to '%s'\n",
            corrfile->DataFilename().full());

  return Analysis::OK;
}

// Analysis_Rms2d::Analyze()
Analysis::RetType Analysis_Rms2d::Analyze() {
  int err = 0;
  // Set up target mask
  if (TgtTraj_->Top().SetupIntegerMask( TgtMask_ )) return Analysis::ERR;
  TgtMask_.MaskInfo();
  if (TgtMask_.None()) { 
    mprinterr("Error: No atoms selected for [%s]\n", TgtMask_.MaskString());
    return Analysis::ERR;
  }
  // Set up reference mask. If no reference parm use target parm.
  if (RefParm_ == 0)
    err = TgtTraj_->Top().SetupIntegerMask( RefMask_ );
  else
    err = RefParm_->SetupIntegerMask(RefMask_);
  if (err != 0) {
    mprinterr("Error: Could not set up reference mask [%s] for parm %s\n",
              RefMask_.MaskString(), RefParm_->c_str());
    return Analysis::ERR;
  }
  if (TgtMask_.MaskExpression() != RefMask_.MaskExpression())
    RefMask_.MaskInfo();
  // Ensure # ref atoms == # tgt atoms
  if (RefMask_.Nselected() != TgtMask_.Nselected()) {
    mprinterr("Error: # Selected atoms in '%s' not equal to selected # atoms in\n"
              "Error:   '%s' (%i)\n", TgtMask_.MaskString(), RefMask_.MaskString(),
              RefMask_.Nselected());
    return Analysis::ERR;
  }
  // Set up symmetry-corrected RMSD calc if necessary - always fit!
  if (mode_ == SRMSD) {
    SRMSD_.InitSymmRMSD(true, useMass_, 0);
    if (SRMSD_.SetupSymmRMSD( TgtTraj_->Top(), TgtMask_, false)) // No remap
      return Analysis::ERR;
  }
  // Actually do the 2D rms calculation
  err = Calculate_2D();
  if (err != 0) return Analysis::ERR;
  return Analysis::OK;
}

// Analysis_Rms2d::Calculate_2D()
/** Perform specified calculation between each frame in Coords to each other
  * frame. If reference and target masks are equal this results in a symmetric
  * matrix, otherwise a full matrix is needed.
  */
int Analysis_Rms2d::Calculate_2D() {
  float R = 0.0;
  int nref, ntgt, tgtstart;
  if (!useReferenceTraj_)
    RefTraj_ = TgtTraj_;
  int totalref = RefTraj_->Size();
  int totaltgt = TgtTraj_->Size();
  // Determine if target mask and reference mask are equial.
  bool masksAreEqual = (TgtMask_ == RefMask_);
  // If tgt and ref masks are the same and ref traj is tgt traj, only need an
  // upper-triangle matrix, otherwise need a full one.
  bool calculateFullMatrix = (!masksAreEqual || useReferenceTraj_);
  if (calculateFullMatrix)
    rmsdataset_->Allocate2D( totalref, totaltgt );
  else
    rmsdataset_->AllocateTriangle( totaltgt );
  // Print Info
  mprintf("\tCalculating %s", ModeStrings_[mode_]);
  if (masksAreEqual) {
    mprintf(" using mask [%s]", TgtMask_.MaskString());
    if (useReferenceTraj_)
      mprintf(" between frames in '%s' and frames in '%s'", TgtTraj_->legend(), RefTraj_->legend());
    else
      mprintf(" between each frame in '%s'", TgtTraj_->legend());
  } else {
    if (useReferenceTraj_)
      mprintf(" between frames in '%s' [%s] and frames in '%s' [%s]",
              TgtTraj_->legend(), TgtMask_.MaskString(),
              RefTraj_->legend(), RefMask_.MaskString());
    else
      mprintf(" between frames in '%s' [%s] to [%s]",
              TgtTraj_->legend(), TgtMask_.MaskString(), RefMask_.MaskString());
  }
  mprintf(" (%zu total).\n", rmsdataset_->Size());
  // Set up target and reference frames based on mask. Both have same topology.
  Frame SelectedRef, SelectedTgt;
  SelectedRef.SetupFrameFromMask( RefMask_, RefTraj_->Top().Atoms() );
  SelectedTgt.SetupFrameFromMask( TgtMask_, TgtTraj_->Top().Atoms() );
# ifdef _OPENMP
  Frame RefFrame;
  bool isTRAJ = false;
  if (RefTraj_->Type() == DataSet::TRAJ) {
    RefFrame = RefTraj_->AllocateFrame();
    isTRAJ = true;
  }
# endif
  ParallelProgress progress( totalref );
  // LOOP OVER REFERENCE FRAMES
# ifdef _OPENMP
  SymmetricRmsdCalc SRMSD_OMP = SRMSD_;
# define SRMSD_ SRMSD_OMP
# pragma omp parallel private(nref, ntgt, tgtstart) firstprivate(SRMSD_OMP, R, SelectedTgt, SelectedRef, RefFrame, progress)
  {
    if (omp_get_thread_num()==0)
      mprintf("\tParallelizing calculation with %i OpenMP threads.\n", omp_get_num_threads());
    progress.SetThread(omp_get_thread_num());
#   pragma omp for schedule(dynamic)
# endif
    for (nref=0; nref < totalref; nref++) {
      progress.Update(nref);
      // Get the current reference frame
#     ifdef _OPENMP
      if (isTRAJ) {
        RefTraj_->GetFrame( nref, RefFrame );
        SelectedRef.SetCoordinates( RefFrame, RefMask_ );
      } else
#     endif
        RefTraj_->GetFrame( nref, SelectedRef, RefMask_ );
      // Select and pre-center reference atoms (if fitting)
      if (mode_ == RMS_FIT || mode_ == SRMSD)
        SelectedRef.CenterOnOrigin(useMass_);
      // LOOP OVER TARGET FRAMES
      if (calculateFullMatrix)
        tgtstart = 0;
      else
        tgtstart = nref + 1;
      for (ntgt = tgtstart; ntgt < totaltgt; ntgt++) {
        // Get the current target frame
        TgtTraj_->GetFrame( ntgt, SelectedTgt, TgtMask_ );
        switch (mode_) {
          case RMS_FIT:   R = (float)SelectedTgt.RMSD_CenteredRef(SelectedRef, useMass_); break;
          case RMS_NOFIT: R = (float)SelectedTgt.RMSD_NoFit(SelectedRef, useMass_); break;
          case DME:       R = (float)SelectedTgt.DISTRMSD(SelectedRef); break;
          case SRMSD:     R = (float)SRMSD_.SymmRMSD_CenteredRef(SelectedTgt, SelectedRef); break;
        }
        rmsdataset_->SetElement(nref, ntgt, R);
        // DEBUG
        //mprinterr("%12i %12i %12.4lf\n",nref,ntgt,R);
      } // END loop over target frames
    } // END loop over reference frames
# ifdef _OPENMP
# undef SRMSD_
  }
# endif
  progress.Finish();
  if (Ct_ != 0) CalcAutoCorr( );
  return 0;
}

/** Calculate the pseudo-autocorrelation of the RMSDs. For proper weighting
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
