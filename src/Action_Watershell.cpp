#include "Action_Watershell.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

// CONSTRUCTOR
Action_Watershell::Action_Watershell() :
  lowerCutoff_(0),
  upperCutoff_(0),
  CurrentParm_(0),
  lower_(0),
  upper_(0)
{}

// Action_Watershell::Help()
void Action_Watershell::Help() const {
  mprintf("\t<solutemask> [out <filename>] [lower <lower cut>] [upper <upper cut>]\n"
          "\t[noimage] [<solventmask>] [<set name>]\n"
          "  Calculate # of waters in 1st and 2nd solvation shells (defined by\n"
          "  <lower cut> (default 3.4 Ang.) and <upper cut> (default 5.0 Ang.)\n"
          "  distance cut-offs respectively.\n");
}

// Action_Watershell::Init()
Action::RetType Action_Watershell::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  image_.InitImaging( !actionArgs.hasKey("noimage") );
  // Get keywords
  std::string filename = actionArgs.GetStringKey("out");
  lowerCutoff_ = actionArgs.getKeyDouble("lower", 3.4);
  upperCutoff_ = actionArgs.getKeyDouble("upper", 5.0);
  // Get solute mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: Solute mask must be specified.\n");
    return Action::ERR;
  }
  soluteMask_.SetMaskString( maskexpr );
  // Check for solvent mask
  std::string solventmaskexpr = actionArgs.GetMaskNext();
  if (!solventmaskexpr.empty())
    solventMask_.SetMaskString( solventmaskexpr );
  // For backwards compat., if no 'out' assume next string is file name 
  if (filename.empty() && actionArgs.Nargs() > 2 && !actionArgs.Marked(2))
    filename = actionArgs.GetStringNext();
  DataFile* outfile = init.DFL().AddDataFile( filename, actionArgs );

  // Set up datasets
  std::string dsname = actionArgs.GetStringNext();
  if (dsname.empty())
    dsname = init.DSL().GenerateDefaultName("WS");
  lower_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, "lower"));
  upper_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, "upper"));
  if (lower_ == 0 || upper_ == 0) return Action::ERR;
  if (outfile != 0) {
    outfile->AddDataSet(lower_);
    outfile->AddDataSet(upper_);
  }
# ifdef _OPENMP
  // Determine number of parallel threads
  int numthreads = 0;
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    numthreads = omp_get_num_threads();
}
  shellStatus_thread_.resize( numthreads );
# endif
  mprintf("    WATERSHELL:");
  if (outfile != 0) mprintf(" Output to %s", outfile->DataFilename().full());
  mprintf("\n");
  if (!image_.UseImage())
    mprintf("\tImaging is disabled.\n");
  mprintf("\tThe first shell will contain solvent < %.3f angstroms from\n",
          lowerCutoff_);
  mprintf("\t  the solute; the second shell < %.3f angstroms...\n",
          upperCutoff_);
  mprintf("\tSolute atoms will be specified by [%s]\n",soluteMask_.MaskString());
  if (solventMask_.MaskStringSet())
    mprintf("\tSolvent atoms will be specified by [%s]\n", solventMask_.MaskString());
# ifdef _OPENMP
  if (shellStatus_thread_.size() > 1)
    mprintf("\tParallelizing calculation with %zu threads.\n", shellStatus_thread_.size());
# endif
  mprintf("\t# solvent molecules in 'lower' shell stored in set '%s'\n", lower_->legend());
  mprintf("\t# solvent molecules in 'upper' shell stored in set '%s'\n", upper_->legend());

  // Pre-square upper and lower cutoffs
  lowerCutoff_ *= lowerCutoff_;
  upperCutoff_ *= upperCutoff_;

  return Action::OK;
}

// Action_Watershell::Setup()
/** Set up solute and solvent masks. If no solvent mask was specified use 
  * solvent information in the current topology.
  */
Action::RetType Action_Watershell::Setup(ActionSetup& setup) {
  // Set up solute mask
  if (setup.Top().SetupIntegerMask( soluteMask_ )) return Action::ERR;
  if ( soluteMask_.None() ) {
    mprintf("Warning: No atoms in solute mask [%s].\n",soluteMask_.MaskString());
    return Action::SKIP;
  }
  if (solventMask_.MaskStringSet()) {
    // Set up solvent mask
    if (setup.Top().SetupIntegerMask( solventMask_ )) return Action::ERR;
  } else {
    // Use all solvent atoms.
    solventMask_.ResetMask();
    for (Topology::mol_iterator mol = setup.Top().MolStart();
                                mol != setup.Top().MolEnd(); ++mol)
      if ( mol->IsSolvent() )
        solventMask_.AddAtomRange( mol->BeginAtom(), mol->EndAtom() );
  }
  if ( solventMask_.None() ) {
    if ( solventMask_.MaskStringSet() )
      mprintf("Warning: No solvent atoms selected by mask [%s]\n", solventMask_.MaskString());
    else
      mprintf("Warning: No solvent atoms in topology %s\n", setup.Top().c_str());
    return Action::SKIP;
  }
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf("\tImaging is on.\n");
  else
    mprintf("\tImaging is off.\n");
  // Allocate space to record status of each solvent molecule.
  // NOTE: Doing this by residue instead of by molecule does waste some memory,
  //       but it means watershell can be used even if no molecule info present. 
# ifdef _OPENMP
  // Each thread needs space to record residue status to avoid clashes
  for (std::vector<Iarray>::iterator it = shellStatus_thread_.begin();
                                     it != shellStatus_thread_.end(); ++it)
    it->assign( setup.Top().Nres(), 0 );
# else
  shellStatus_.assign( setup.Top().Nres(), 0 );
# endif
  // Store current topology
  CurrentParm_ = setup.TopAddress();
  return Action::OK;    
}

// Action_Watershell::action()
Action::RetType Action_Watershell::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 ucell, recip;
  int nlower = 0;
  int nupper = 0;

  if (image_.ImageType()==NONORTHO) frm.Frm().BoxCrd().ToRecip(ucell,recip);

  int Vidx, Uidx, Vat, currentRes;
  int NU = soluteMask_.Nselected();
  int NV = solventMask_.Nselected();
  double dist;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(dist,Vidx,Vat,currentRes,Uidx,mythread)
  {
  mythread = omp_get_thread_num();
# pragma omp for
# endif
  // Assume solvent mask is the larger one.
  // Loop over solvent atoms
  for (Vidx = 0; Vidx < NV; Vidx++) {
    // Figure out which solvent residue this is
    Vat = solventMask_[Vidx];
    currentRes = (*CurrentParm_)[ Vat ].ResNum();
    // Loop over solute atoms
    for (Uidx = 0; Uidx < NU; Uidx++) {
      // If residue is not yet marked as 1st shell, calc distance
#     ifdef _OPENMP
      if ( shellStatus_thread_[mythread][currentRes] < 2 )
#     else
      if ( shellStatus_[currentRes] < 2 )
#     endif
      {
        dist = DIST2(frm.Frm().XYZ(soluteMask_[Uidx]), frm.Frm().XYZ(Vat),
                     image_.ImageType(), frm.Frm().BoxCrd(), ucell, recip );
        // Less than upper, 2nd shell
        if (dist < upperCutoff_) 
        {
#         ifdef _OPENMP
          shellStatus_thread_[mythread][currentRes] = 1;
#         else
          shellStatus_[currentRes] = 1;
#         endif
          // Less than lower, 1st shell
          if (dist < lowerCutoff_)
#           ifdef _OPENMP
            shellStatus_thread_[mythread][currentRes] = 2;
#           else
            shellStatus_[currentRes] = 2;
#           endif
        }
      }
    } // End loop over solute atoms
  } // End loop over solvent atoms
# ifdef _OPENMP
  } // END parallel
  // Combine results from each thread.
  for (unsigned int res = 0; res < shellStatus_thread_[0].size(); res++) {
    int shell = 0;
    for (unsigned int thread = 0; thread < shellStatus_thread_.size(); thread++) {
      shell = std::max( shellStatus_thread_[thread][res], shell );
      // Dont break here so we can reset counts. Could also do with a fill 
      shellStatus_thread_[thread][res] = 0;
    }
    if (shell > 0) {
      ++nupper;
      if (shell > 1) ++nlower;
    }
  }
# else
  // Now each residue is marked 0 (no shell), 1 (second shell), 2 (first shell)
  for (Iarray::iterator shell = shellStatus_.begin(); shell != shellStatus_.end(); ++shell)
  {
    if ( *shell > 0 ) {
      ++nupper;
      if ( *shell > 1 ) ++nlower;
    }
    // Reset for next pass
    *shell = 0;
  }
# endif
  lower_->Add(frameNum, &nlower);
  upper_->Add(frameNum, &nupper);

  return Action::OK;
}
