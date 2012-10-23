#include <cmath> // exp
#include "Action_Rms2d.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "Matrix_2D.h"

// CONSTRUCTOR
Action_Rms2d::Action_Rms2d() :
  nofit_(false),
  useMass_(false),
  RefTraj_(NULL),
  RefParm_(NULL),
  mass_ptr_(NULL),
  rmsdataset_(0),
  Ct_(0)
{ } 

void Action_Rms2d::Help() {
  mprintf("rms2d <mask> <refmask> rmsout <filename> [nofit]\n");
  mprintf("      [reftraj <traj> [parm <parmname> | parmindex <#>]]\n");
  mprintf("      [corr <corrfilename>]\n");
}

// DESTRUCTOR
Action_Rms2d::~Action_Rms2d() {
  if (RefTraj_!=NULL) delete RefTraj_; 
}

// Action_Rms2d::init()
Action::RetType Action_Rms2d::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  nofit_ = actionArgs.hasKey("nofit");
  useMass_ = actionArgs.hasKey("mass"); 
  std::string rmsdFile = actionArgs.GetStringKey("rmsout");
  std::string reftraj = actionArgs.GetStringKey("reftraj");
  if (!reftraj.empty()) {
    RefParm_ = PFL->GetParm(actionArgs);
    if (RefParm_==NULL) {
      mprinterr("Error: Rms2d: Could not get parm for reftraj %s.\n",reftraj.c_str());
      return Action::ERR;
    }
  }
  // Check for correlation; if so, reftraj not supported
  std::string corrfilename = actionArgs.GetStringKey("corr");
  if (!corrfilename.empty() && !reftraj.empty()) {
    mprinterr("Error: Rms2d: Keyword 'corr' not supported with 'reftraj'\n");
    return Action::ERR;
  }
  // Require an output filename if corr not specified
  if (rmsdFile.empty() && corrfilename.empty()) {
    mprinterr("Error: Rms2d: No output filename specified; use 'rmsout' keyword.\n");
    return Action::ERR;
  }

  // Get the RMS mask string for frames
  ArgList::ConstArg mask0 = actionArgs.getNextMask();
  FrameMask_.SetMaskString(mask0);

  // Check if reference will be a series of frames from a trajectory
  if (!reftraj.empty()) {
    // Get RMS mask string for reference trajectory
    ArgList::ConstArg maskRef = actionArgs.getNextMask();
    // If no reference mask specified, make same as RMS mask
    if (maskRef==NULL) maskRef=mask0;
    RefMask_.SetMaskString(maskRef);
    // Attempt to set up reference trajectory
    RefTraj_ = new TrajectoryFile();
    if (RefTraj_->SetupTrajRead(reftraj, NULL, RefParm_)) {
      mprinterr("Error: Rms2d: Could not set up reftraj %s.\n",reftraj.c_str());
      return Action::ERR;
    }
    // Set up DataSet for RefTraj, full 2d matrix
    rmsdataset_ = DSL->AddSet( DataSet::MATRIX2D, rmsdFile, "Rms2d" );
  } else {
    // Set up DataSet for normal rms2d, half matrix
    rmsdataset_ = DSL->AddSet( DataSet::TRIMATRIX, rmsdFile, "Rms2d" );
    // Set up DataSet and DataFile for correlate
    if (!corrfilename.empty()) {
      Ct_ = DSL->AddSetAspect( DataSet::DOUBLE, rmsdFile, "Corr" );
      if (Ct_ == NULL) return Action::ERR;
      DFL->AddSetToFile(corrfilename, Ct_);
    }
  }
  if (rmsdataset_ == NULL) {
    mprinterr("Error: Could not set up DataSet for calculating 2DRMS.\n");
    return Action::ERR;
  }
  // Format DataSet and set up output file
  rmsdataset_->SetPrecision(8,3);
  DataFile* rmsdatafile = DFL->AddSetToFile( rmsdFile, rmsdataset_ );
  rmsdatafile->ProcessArgs("square2d");

  mprintf("    RMS2D: Mask [%s]",FrameMask_.MaskString());
  if (!reftraj.empty()) {
    // Set up reference trajectory and open
    mprintf(", ref traj %s (mask [%s]) %i frames", RefTraj_->FullTrajStr(),
            RefMask_.MaskString(), RefTraj_->Total_Read_Frames());
  }
  if (nofit_)
    mprintf(" (no fitting)");
  if (useMass_)
    mprintf(" (mass-weighted)");
  if (!rmsdFile.empty()) 
    mprintf(" output to %s",rmsdFile.c_str());
  mprintf("\n");
  if (!corrfilename.empty())
    mprintf("           RMSD auto-correlation will be calculated and output to %s\n",
            corrfilename.c_str());

  return Action::OK;
}

// Action_Rms2d::setup()
/** Set up frame mask so that only selected atoms in frames will be stored.
  */
Action::RetType Action_Rms2d::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask(FrameMask_) ) {
    mprinterr("Error: Rms2d::setup: Could not set up mask [%s] for parm %s\n",
              FrameMask_.MaskString(), currentParm->c_str());
    return Action::ERR;
  }
  if (FrameMask_.None()) {
    mprinterr("Error: Rms2d::setup: No atoms selected for mask [%s], parm %s\n",
              FrameMask_.MaskString(), currentParm->c_str());
    return Action::ERR;
  }
  // If useMass, store mass information here. If mass information changes
  // (i.e. a different parm than the first is loaded) this makes calculating
  // the 2drms much more complicated since each frame could have a different
  // set of masses. To simplify things, only allow mass-weighting when reading
  // one parm.
  if (useMass_) {
    if (mass_ptr_ == NULL) {
      mass_ptr_ = currentParm;
    } else {
      mprintf("Warning: Rms2d::Setup: 'mass' only allowed with one parm. Disabling 'mass'.\n");
      mass_ptr_=NULL;
    }
  }
  return Action::OK;  
}

// Action_Rms2d::action()
/** Store current frame coords according to mask.
  */
Action::RetType Action_Rms2d::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  if (ReferenceCoords_.AddFrameByMask(*currentFrame, FrameMask_)) 
    return Action::ERR;

  return Action::OK;
} 

// Action_Rms2d::Calc2drms()
/** Calculate the RMSD of each frame in ReferenceCoords to each other frame.
  * Since this results in a symmetric matrix use TriangleMatrix to store
  * results.
  */
int Action_Rms2d::Calc2drms() {
  Frame RefFrame;
  Frame TgtFrame;
  double U[9], Trans[6];
  float R;
 
  TriangleMatrix* Distances = (TriangleMatrix*) rmsdataset_;
  if ( Distances == NULL ) { // Sanity check
    mprinterr("Error: Could not set up DataSet for calculating 2DRMS.\n");
    return 1;
  }

  int totalref = ReferenceCoords_.Ncoords();
  Distances->Setup( totalref );
  int max = Distances->Nelements();
  mprintf("  RMS2D: Calculating RMSDs between each frame (%i total).\n  ",max);

  // Set up progress Bar
  ProgressBar progress(totalref - 1);

  if (mass_ptr_!=NULL && useMass_) {
    // Set up mass info. If mass info present this means only 1 parm used,
    // so tgt and ref will always have same # of atoms. Use first frame
    // in ref coords to set up target and reference frames.
    RefFrame.SetupFrameFromMask(FrameMask_, mass_ptr_->Atoms());
    TgtFrame.SetupFrameFromMask(FrameMask_, mass_ptr_->Atoms());
    // If no mass info ensure that RefFrame and TgtFrame are large enough to
    // hold the largest set of coords in ReferenceCoords. No mass.
  } else {
    int maxrefnatom = ReferenceCoords_.MaxNatom();
    RefFrame.SetupFrame( maxrefnatom );
    TgtFrame.SetupFrame( maxrefnatom );
  }

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref - 1; nref++) {
    progress.Update(nref);
    // Get the current reference frame
    RefFrame = ReferenceCoords_[nref];
    // Select and pre-center reference atoms (if fitting)
    if (!nofit_)
      RefFrame.CenterReference(Trans+3, useMass_);
  
    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      // Get the current target frame
      TgtFrame = ReferenceCoords_[nframe];
      // Ensure # ref atoms == # tgt atoms
      if (RefFrame.Natom() != TgtFrame.Natom()) {
        mprintf("\tWarning: Rms2d: # atoms in ref %i (%i) != # atoms in tgt %i (%i)\n",
                nref+1,RefFrame.Natom(),nframe+1,TgtFrame.Natom());
        mprintf("\t         Assigning this pair RMSD of -1.0\n");
        R = -1.0;
      } else if (nofit_) {
        // Perform no fit RMS calculation
        R = (float) TgtFrame.RMSD(RefFrame, useMass_);
      } else {
        // Perform fit RMS calculation
        R = (float) TgtFrame.RMSD_CenteredRef(RefFrame, U, Trans, useMass_);
      }
      Distances->AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
    } // END loop over target frames
  } // END loop over reference frames
  // Calculate correlation if specified
  if (Ct_ != 0) AutoCorrelate( *Distances );

  return 0;
}

// Action_Rms2d::CalcRmsToTraj()
/** Calc RMSD of every frame in reference traj to every frame in 
  * ReferenceCoords.
  */
int Action_Rms2d::CalcRmsToTraj() {
  Frame RefFrame;
  Frame SelectedTgt;
  Frame SelectedRef;
  double U[9], Trans[6];
  float R;

  Matrix_2D* rmsdata = (Matrix_2D*)rmsdataset_;
  if ( rmsdata == NULL ) { // Sanity check
    mprinterr("Error: Could not set up DataSet for calculating 2DRMS to traj.\n");
    return 1;
  }

  // Set up reference mask for reference parm
  if (RefParm_->SetupIntegerMask(RefMask_)) {
    mprinterr("Error: Could not set up reference mask [%s] for parm %s\n",
              RefMask_.MaskString(), RefParm_->c_str());
    return 1;
  }
  // Setup frame for selected reference atoms
  SelectedRef.SetupFrameFromMask(RefMask_, RefParm_->Atoms()); 
  RefFrame.SetupFrameM(RefParm_->Atoms());
  int totalref = RefTraj_->Total_Read_Frames();

  int totaltgt = ReferenceCoords_.Ncoords();
  int max = totalref * totaltgt;
  mprintf("  RMS2D: Calculating RMSDs between each input frame and each reference\n"); 
  mprintf("         trajectory %s frame (%i total).\n  ",
          RefTraj_->BaseTrajStr(), max);
  rmsdata->Setup( totalref, totaltgt );
  if (RefTraj_->BeginTraj(false)) {
    mprinterr("Error: Rms2d: Could not open reference trajectory.\n");
    return 1;
  }
  // Set up progress Bar
  ProgressBar progress(totalref);

  if (mass_ptr_!=NULL && useMass_) {
    // Set up selected target mass info
    SelectedTgt.SetupFrameFromMask(FrameMask_, mass_ptr_->Atoms());
  } else {
    // If no mass, ensure SelectedTgt can hold max #atoms in ReferenceCoords
    int maxtgtnatom = ReferenceCoords_.MaxNatom();
    SelectedTgt.SetupFrame( maxtgtnatom );
  }
  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref; nref++) {
    progress.Update(nref);
    // Get the current reference frame from trajectory
    RefTraj_->GetNextFrame(RefFrame);
  
    // Set up dataset for this reference frame
    /*std::string setname = "Frame_" + integerToString( nref+1 );
    DataSet* rmsdata = DSL->Add( DataSet::FLOAT, setname.c_str(), "Rms2d" );
    if (rmsdata==NULL) {
      mprinterr("Error: rms2d: Could not add dataset for ref frame %i\n", nref+1);
      return;
    }
    DFL->Add(rmsdFile_, rmsdata);*/
    // Set reference atoms and pre-center if fitting
    SelectedRef.SetCoordinates(RefFrame, RefMask_);
    if (!nofit_)
      SelectedRef.CenterReference(Trans+3, useMass_);

    // LOOP OVER TARGET FRAMES
    for (int nframe=0; nframe < totaltgt; nframe++) {
      // Get selected atoms of the current target frame
      SelectedTgt = ReferenceCoords_[nframe];
      // Ensure # ref atoms == # tgt atoms
      if (SelectedRef.Natom() != SelectedTgt.Natom()) {
        mprintf("\tWarning: Rms2d: Selected # atoms in ref %i (%i) != selected # atoms\n",
                nref+1, SelectedRef.Natom());
        mprintf("\t         in tgt %i (%i).Assigning this pair RMSD of -1.0\n",
                nframe+1, SelectedTgt.Natom());
        R = -1.0;
      } else if (nofit_) {
        // Perform no fit RMS calculation
        R = (float) SelectedTgt.RMSD(SelectedRef, useMass_);
      } else {
        // Perform fit RMS calculation
        R = (float) SelectedTgt.RMSD_CenteredRef(SelectedRef, U, Trans, useMass_);
      }
      //rmsdata->Add(nframe, &R);
      rmsdata->AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
    } // END loop over target frames
  } // END loop over reference frames
  RefTraj_->EndTraj();
  return 0;
}

// Action_Rms2d::AutoCorrelate()
/** Calculate the autocorrelation of the RMSDs. For proper weighting
  * exp[ -RMSD(framei, framei+lag) ] is used. This takes advantage of
  * the fact that 0.0 RMSD essentially means perfect correlation (1.0).
  */
int Action_Rms2d::AutoCorrelate(TriangleMatrix& Distances) {
  double ct;
  int lagmax = ReferenceCoords_.Ncoords();
  int N = ReferenceCoords_.Ncoords();

  // By definition for lag == 0 RMS is 0 for all frames,
  // translates to correlation of 1.
  ct = 1;
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

  return 0;
}

// Action_Rms2d::print()
/** Perform the rms calculation of each frame to each other frame.
  */
void Action_Rms2d::Print() {
  if (RefTraj_==NULL) 
    Calc2drms();
  else
    CalcRmsToTraj();
}

