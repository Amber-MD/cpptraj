#include <cmath> // floor
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "Action_Watershell.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"

#ifdef CUDA
// CUDA Kernel wrappers
extern void Action_Closest_NoCenter(const double*,double*,const double*,double,int,int,int,ImagingType,const double*,const double*,const double*);
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
# ifndef CUDA
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
#ifdef CUDA
  mprintf("\tDistance calculations will be GPU-accelerated with CUDA.\n");
#else
# ifdef _OPENMP
  if (shellStatus_thread_.size() > 1)
    mprintf("\tParallelizing calculation with %zu threads.\n", shellStatus_thread_.size());
# endif
#endif
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
    // Set number of atoms; needed for CharMask conversion.
    solventMask_.SetNatoms( setup.Top().Natom() );
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
#ifdef CUDA
  // Since we are using the 'closest' kernels under the hood, all solvent mols
  // must have the same size.
  int first_solvent_mol = setup.Top()[ solventMask_[0] ].MolNum();
  NAtoms_ = setup.Top().Mol( first_solvent_mol ).NumAtoms();
  for (AtomMask::const_iterator atm = solventMask_.begin(); atm != solventMask_.end(); ++atm) {
    int mol = setup.Top()[*atm].MolNum();
    if (NAtoms_ != setup.Top().Mol( mol ).NumAtoms()) {
      mprinterr("Error: CUDA version of 'watershell' requires all solvent mols be same size.\n");
      return Action::ERR;
    }
  }
  // Determine how many solvent molecules are selected
  NsolventMolecules_ = 0;
  CharMask cMask( solventMask_.ConvertToCharMask(), solventMask_.Nselected() );
  for (Topology::mol_iterator mol = setup.Top().MolStart(); mol != setup.Top().MolEnd(); ++mol)
    if ( cMask.AtomsInCharMask( mol->BeginAtom(), mol->EndAtom() ) )
      NsolventMolecules_++;
  // Sanity check
  if ( (NsolventMolecules_ * NAtoms_) != solventMask_.Nselected() ) {
    mprinterr("Error: CUDA version of 'watershell' requires all atoms in solvent mols be selected.\n");
    return Action::ERR;
  }
  // Allocate space for selected solvent atom coords and distances
  V_atom_coords_.resize( NsolventMolecules_ * NAtoms_ * 3, 0.0 );
  V_distances_.resize( NsolventMolecules_ );
#else /* CUDA */
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
#endif /* CUDA */
  // Set up imaging
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf("\tImaging is on.\n");
  else
    mprintf("\tImaging is off.\n");
  // Allocate temp space for selected solute atom coords.
  soluteCoords_.resize( soluteMask_.Nselected() * 3 );
  // Store current topology
  CurrentParm_ = setup.TopAddress();
  return Action::OK;    
}

// Action_Watershell::DoAction()
Action::RetType Action_Watershell::DoAction(int frameNum, ActionFrame& frm) {
  int nlower = 0;
  int nupper = 0;
# ifdef CUDA
  // Copy solvent atom coords to array
  unsigned int idx = 0; // Index into V_atom_coords_
  for (AtomMask::const_iterator atm = solventMask_.begin();
                                atm != solventMask_.end(); ++atm, idx += 3)
  {
    const double* xyz = frm.Frm().XYZ( *atm );
    V_atom_coords_[idx  ] = xyz[0];
    V_atom_coords_[idx+1] = xyz[1];
    V_atom_coords_[idx+2] = xyz[2];
  }
  Matrix_3x3 ucell, recip;
  if (image_.ImageType() == NONORTHO)
    frm.Frm().BoxCrd().ToRecip(ucell, recip);
  // Copy solute atom coords to array
  idx = 0;
  for (AtomMask::const_iterator atm = soluteMask_.begin();
                                atm != soluteMask_.end(); ++atm, idx += 3)
  {
    const double* XYZ = frm.Frm().XYZ( *atm );
    soluteCoords_[idx  ] = XYZ[0];
    soluteCoords_[idx+1] = XYZ[1];
    soluteCoords_[idx+2] = XYZ[2];
  }

  Action_Closest_NoCenter( &V_atom_coords_[0], &V_distances_[0], &soluteCoords_[0],
                           9999999999999.0, NsolventMolecules_, NAtoms_, soluteMask_.Nselected(),
                           image_.ImageType(),
                           frm.Frm().BoxCrd().boxPtr(), ucell.Dptr(), recip.Dptr() );

  // V_distances_ now has the closest distance of each solvent molecule to
  // solute. Determine shell status of each.
  for (Darray::const_iterator dist2 = V_distances_.begin(); dist2 != V_distances_.end(); ++dist2)
    if (*dist2 < upperCutoff_)
    {
      ++nupper;
      // Less than lower, 1st shell
      if (*dist2 < lowerCutoff_)
        ++nlower;
    }

# else
  // ---------------------------------------------------------------------------
  int NV = solventMask_.Nselected();
  int Vidx;
# ifdef _OPENMP
  int mythread;
# endif
  int* status = 0;

  if (image_.ImageType() == NONORTHO) {
    // ----- NON-ORTHORHOMBIC IMAGING ------------
    Matrix_3x3 ucell, recip;
    frm.Frm().BoxCrd().ToRecip(ucell, recip);
    // Wrap all solute atoms back into primary cell, save coords.
    Image::WrapToCell0( soluteCoords_, frm.Frm(), soluteMask_, ucell, recip );
    // Calculate every imaged distance of all solvent atoms to solute
#   ifdef _OPENMP
#   pragma omp parallel private(Vidx, mythread, status)
    {
    mythread = omp_get_thread_num();
    status = &(shellStatus_thread_[mythread][0]);
#   pragma omp for
#   else
    status = &shellStatus_[0];
#   endif
    for (Vidx = 0; Vidx < NV; Vidx++)
    {
      int Vat = solventMask_[Vidx];
      int currentRes = (*CurrentParm_)[ Vat ].ResNum();
      // Convert to frac coords
      Vec3 vFrac = recip * Vec3( frm.Frm().XYZ( Vat ) );
      // Wrap to primary unit cell
      vFrac[0] = vFrac[0] - floor(vFrac[0]);
      vFrac[1] = vFrac[1] - floor(vFrac[1]);
      vFrac[2] = vFrac[2] - floor(vFrac[2]);
      // Loop over all images of this solvent atom
      for (int ix = -1; ix != 2; ix++)
        for (int iy = -1; iy != 2; iy++)
          for (int iz = -1; iz != 2; iz++)
          {
            // Convert image back to Cartesian
            Vec3 vCart = ucell.TransposeMult( vFrac + Vec3(ix, iy, iz) );
            // Loop over all solute atoms
            for (unsigned int idx = 0; idx < soluteCoords_.size(); idx += 3)
            {
              if (status[currentRes] < 2)
              {
                // Residue is not yet marked as 1st shell, calc distance
                double x = vCart[0] - soluteCoords_[idx  ];
                double y = vCart[1] - soluteCoords_[idx+1];
                double z = vCart[2] - soluteCoords_[idx+2];
                double dist2 = x*x + y*y + z*z;
                // Less than upper, 2nd shell
                if (dist2 < upperCutoff_) 
                {
                  status[currentRes] = 1;
                  // Less than lower, 1st shell
                  if (dist2 < lowerCutoff_)
                    status[currentRes] = 2;
                }
              } else {
                // Residue already in first shell. No need for more distance calcs
                ix = iy = iz = 1;
                break;
              }
            } // END loop over solute atoms 
          } // END loop over images (Z)
    } // END loop over solvent atoms
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  } else {
    // ----- ORTHORHOMBIC/NO IMAGING -------------
    // Store selected solute coordinates.
    int idx = 0;
    for (AtomMask::const_iterator atm = soluteMask_.begin();
                                  atm != soluteMask_.end(); ++atm, idx += 3)
    {
      const double* XYZ = frm.Frm().XYZ( *atm );
      soluteCoords_[idx  ] = XYZ[0];
      soluteCoords_[idx+1] = XYZ[1];
      soluteCoords_[idx+2] = XYZ[2];
    }
    // Calculate distance of all solvent atoms to solute
#   ifdef _OPENMP
#   pragma omp parallel private(Vidx, mythread, status)
    {
    mythread = omp_get_thread_num();
    status = &(shellStatus_thread_[mythread][0]);
#   pragma omp for
#   else
    status = &shellStatus_[0];
#   endif
    for (Vidx = 0; Vidx < NV; Vidx++)
    {
      int Vat = solventMask_[Vidx];
      int currentRes = (*CurrentParm_)[ Vat ].ResNum();
      Vec3 Vcoord( frm.Frm().XYZ( Vat ) );
      // Loop over all solute atoms
      for (unsigned int idx = 0; idx < soluteCoords_.size(); idx += 3)
      {
        if (status[currentRes] < 2)
        {
          // Residue is not yet marked as 1st shell, calc distance
          Vec3 Ucoord( soluteCoords_[idx], soluteCoords_[idx+1], soluteCoords_[idx+2] );
          double dist2;
          if (image_.ImageType() == ORTHO)
            dist2 = DIST2_ImageOrtho( Vcoord, Ucoord, frm.Frm().BoxCrd() );
          else
            dist2 = DIST2_NoImage( Vcoord, Ucoord );
          // Less than upper, 2nd shell
          if (dist2 < upperCutoff_)
          {
            status[currentRes] = 1;
            // Less than lower, 1st shell
            if (dist2 < lowerCutoff_)
              status[currentRes] = 2;
          }
        } else {
          // Residue already in first shell. No need for more distance calcs.
          break;
        }
      } // END loop over solute atoms
    } // END loop over solvent atoms
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  }
# ifdef _OPENMP
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
# endif /* CUDA */
  lower_->Add(frameNum, &nlower);
  upper_->Add(frameNum, &nupper);

  return Action::OK;
}
