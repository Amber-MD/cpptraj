#include <cmath>
#include <algorithm> // sort
#include <cfloat> // DBL_MAX
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "Action_Closest.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"
#ifdef CUDA
#include "cuda_kernels/kernel_wrappers.cuh"
#else
#include "DistRoutines.h"
#endif

// CONSTRUCTOR
Action_Closest::Action_Closest() :
# ifdef CUDA
  GPU_MEM_(0),
  V_atom_coords_(0),
  U_atom_coords_(0),
  V_distances_(0),
# endif
  outFile_(0),
  framedata_(0),
  moldata_(0),
  distdata_(0),
  atomdata_(0),
  Nclosest_(0),
  closestWaters_(0),
  targetNclosest_(0),
  firstAtom_(false),
  useMaskCenter_(false),
  newParm_(0),
  NsolventMolecules_(0)
{} 

void Action_Closest::Help() const {
  mprintf("\t<# to keep> <mask> [solventmask <solvent mask>] [noimage]\n"
          "\t[first | oxygen] [center] [closestout <filename> [name <setname>]]\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Keep only the closest <# to keep> solvent molecules to atoms in <mask>.\n"
          "  Molecules can be marked as solvent with the 'solvent' command.\n"
          "  If 'center' specified use geometric center of atoms in <mask>.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// DESTRUCTOR
Action_Closest::~Action_Closest() {
  if (newParm_!=0) delete newParm_;
# ifdef CUDA
  if (GPU_MEM_ != 0) delete[] GPU_MEM_;
# endif
}

// Action_Closest::Init()
Action::RetType Action_Closest::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  // Get Keywords
  closestWaters_ = actionArgs.getNextInteger(-1);
  if (closestWaters_ < 0) {
    mprinterr("Error: Invalid # solvent molecules to keep (%i).\n",
              closestWaters_);
    return Action::ERR;
  }
  // Save target # closest in case it is changed in Setup().
  targetNclosest_ = closestWaters_;
  if ( actionArgs.hasKey("oxygen") || actionArgs.hasKey("first") )
    firstAtom_=true;
  useMaskCenter_ = actionArgs.hasKey("center");
  imageOpt_.InitImaging( !(actionArgs.hasKey("noimage")) );
  topWriter_.InitTopWriter(actionArgs, "closest", debugIn);
  // Setup output file and sets if requested.
  // Will keep track of Frame, Mol#, Distance, and first solvent atom
  std::string filename = actionArgs.GetStringKey("closestout");
  if (!filename.empty()) {
    std::string dsetName = actionArgs.GetStringKey("name");
    if (dsetName.empty())
      dsetName = init.DSL().GenerateDefaultName("CLOSEST");
    // Set up datasets
    framedata_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsetName, "Frame"));
    moldata_   = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsetName, "Mol"));
    distdata_  = init.DSL().AddSet(DataSet::DOUBLE,  MetaData(dsetName, "Dist"));
    atomdata_  = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsetName, "FirstAtm"));
    if (framedata_==0 || moldata_==0 || distdata_==0 || atomdata_==0) {
      mprinterr("Error: Could not setup data sets for output file %s\n",
                filename.c_str());
      return Action::ERR;
    }
    // Add sets to datafile in list.
    outFile_ = init.DFL().AddDataFile( filename, "noxcol", actionArgs );
    if (outFile_ == 0) {
      mprinterr("Error: Could not set up output file %s\n", filename.c_str());
      return Action::ERR;
    }
    outFile_->AddDataSet(framedata_);
    outFile_->AddDataSet(moldata_);
    outFile_->AddDataSet(distdata_);
    outFile_->AddDataSet(atomdata_);
  }

  // Get Masks
  std::string mask1 = actionArgs.GetStringKey("solventmask");
  if (!mask1.empty()) {
    if (solventMask_.SetMaskString( mask1 )) {
      mprinterr("Error: Could not set solvent mask string.\n");
      return Action::ERR;
    }
  }
  mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: No mask specified.\n");
    return Action::ERR;
  }
  if (distanceMask_.SetMaskString(mask1)) return Action::ERR;

  mprintf("    CLOSEST: Finding closest %i solvent molecules to atoms in mask %s\n",
          closestWaters_, distanceMask_.MaskString());
  if (useMaskCenter_)
    mprintf("\tGeometric center of atoms in mask will be used.\n");
  if (!imageOpt_.UseImage()) 
    mprintf("\tImaging will be turned off.\n");
  if (solventMask_.MaskStringSet())
    mprintf("\tSolvent will be selected by mask '%s'\n", solventMask_.MaskString());
  if (firstAtom_)
    mprintf("\tOnly first atom of solvent molecule used for distance calc.\n");
  if (outFile_!=0)
    mprintf("\tClosest molecules will be saved to %s\n",outFile_->DataFilename().base());
  topWriter_.PrintOptions();
# ifdef CUDA
  mprintf("\tDistance calculations will be GPU-accelerated with CUDA.\n");
# endif
  return Action::OK;
}

// Action_Closest::Setup()
/** Like the strip action, closest will modify the current parm keeping info
  * for atoms in mask plus the closestWaters solvent molecules. Set up the
  * vector of MolDist objects, one for every solvent molecule in the original
  * parm file. Atom masks for each solvent molecule will be set up.
  */
Action::RetType Action_Closest::Setup(ActionSetup& setup) {
  if (setup.Top().Nmol() < 1) {
    mprintf("Warning: 'closest' requires molecule information. Skipping.\n");
    return Action::SKIP;
  }
  closestWaters_ = targetNclosest_;
  // Determine solvent molecules
  std::vector<bool> IsSolventMol;
  IsSolventMol.reserve( setup.Top().Nmol() );
  int nSolvent = 0;
  if (solventMask_.MaskStringSet()) {
    // Use only atoms selected by solventMask
    if (setup.Top().SetupCharMask( solventMask_ )) return Action::ERR;
    solventMask_.MaskInfo();
    if (solventMask_.None()) {
      mprintf("Warning: No solvent selected by mask '%s'. Skipping.\n", solventMask_.MaskString());
      return Action::SKIP;
    }
    for (Topology::mol_iterator Mol = setup.Top().MolStart();
                                Mol != setup.Top().MolEnd(); ++Mol)
    {
      if ( solventMask_.AtomsInCharMask(Mol->MolUnit()) ) {
        IsSolventMol.push_back( true );
        nSolvent++;
      } else
        IsSolventMol.push_back( false );
    }
  } else {
    // Select everything; all solvent atoms are fair game.
    solventMask_.ResetMask();
    solventMask_.InitCharMask(setup.Top().Natom(), true);
    for (Topology::mol_iterator Mol = setup.Top().MolStart();
                                Mol != setup.Top().MolEnd(); ++Mol)
      IsSolventMol.push_back( Mol->IsSolvent() );
    nSolvent = setup.Top().Nsolvent();
  }
  mprintf("\t%i molecules out of %i selected as solvent.\n", nSolvent, setup.Top().Nmol());
  // If there are no solvent molecules this action is not valid.
  if (nSolvent == 0) {
    mprintf("Warning: Topology %s does not contain solvent.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  if (closestWaters_ == nSolvent) {
    mprintf("Warning: # solvent to keep (%i) == # solvent molecules in '%s' (%i)\n",
            closestWaters_, setup.Top().c_str(), nSolvent);
  } else if (closestWaters_ > nSolvent) {
    mprintf("Warning: # solvent to keep (%i) > # solvent molecules in '%s' (%i)\n",
              closestWaters_, setup.Top().c_str(), nSolvent);
    closestWaters_ = nSolvent;
    mprintf("Warning:  Keeping %i solvent molecules.\n", closestWaters_);
  }
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );
  if (imageOpt_.ImagingEnabled())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tImaging off.\n"); 
  // LOOP OVER MOLECULES
  // 1: Check that all solvent molecules contain same # atoms. Solvent 
  //    molecules must be identical for the command to work properly; 
  //    the prmtop strip occurs only once so the solvent params become fixed.
  // 2: Set up a mask for all solvent molecules.
  SolventMols_.clear();
  // NOTE: May not be necessary to init 'solvent'
  MolDist solvent;
  solvent.D = 0.0;
  solvent.mol = 0;
  SolventMols_.resize(nSolvent, solvent);
  std::vector<MolDist>::iterator mdist = SolventMols_.begin();
  std::vector<bool>::iterator isSolvent = IsSolventMol.begin();
  // 3: Set up the soluteMask for all non-solvent molecules.
  stripMask_.ResetMask();
  int molnum = 1;
  int newnatom = 0;
  int nclosest = 0;
  unsigned int NsolventAtoms = 0;
  keptWaterAtomNum_.resize(closestWaters_);
  for (Topology::mol_iterator Mol = setup.Top().MolStart();
                              Mol != setup.Top().MolEnd(); ++Mol, ++isSolvent)
  {
    if ( !(*isSolvent) ) {
      // Not solvent, add to solute mask.
      stripMask_.AddUnit( Mol->MolUnit() );
      newnatom += Mol->NumAtoms();
    } else {
      // Solvent, check for same # of atoms.
      if (NsolventAtoms == 0)
        NsolventAtoms = Mol->NumAtoms();
      else if ( NsolventAtoms != Mol->NumAtoms() ) {
        mprinterr("Error: Solvent molecules in '%s' are not of uniform size.\n"
                  "Error:   First solvent mol = %u atoms, solvent mol %i = %u atoms.\n",
                  setup.Top().c_str(), NsolventAtoms, molnum, Mol->NumAtoms());
        return Action::ERR;
      }
      // mol here is the output molecule number which is why it starts from 1.
      mdist->mol = molnum;
      // Entire solvent molecule mask
      mdist->mask.AddUnit( Mol->MolUnit() );
      // Atoms in the solvent molecule to actually calculate distances to.
      if (firstAtom_) {
        mdist->solventAtoms.assign(1, Mol->MolUnit().Front() );
      } else {
        mdist->solventAtoms.clear();
        mdist->solventAtoms.reserve( Mol->NumAtoms() );
        for (Unit::const_iterator seg = Mol->MolUnit().segBegin();
                                  seg != Mol->MolUnit().segEnd(); ++seg)
          for (int svatom = seg->Begin(); svatom < seg->End(); svatom++)
            if (solventMask_.AtomInCharMask(svatom))
              mdist->solventAtoms.push_back( svatom );
      }
      // For solvent molecules that will be kept, record what the atom number
      // will be in the new stripped parm.
      if (nclosest < closestWaters_) {
        keptWaterAtomNum_[nclosest] = newnatom;
        stripMask_.AddUnit( Mol->MolUnit() );
        newnatom += Mol->NumAtoms();
        ++nclosest;
      }
      //SolventMols[solventMol].mask.PrintMaskAtoms("solvent");
      ++mdist;
    }
    ++molnum;
  }

  // Setup distance atom mask
  // NOTE: Should ensure that no solvent atoms are selected!
  if ( setup.Top().SetupIntegerMask(distanceMask_) ) return Action::ERR;
  if (distanceMask_.None()) {
    mprintf("Warning: Distance mask '%s' contains no atoms.\n",
            distanceMask_.MaskString());
    return Action::SKIP;
  }
  distanceMask_.MaskInfo();

  // Check the total number of solvent atoms to be kept.
  NsolventAtoms *= closestWaters_;
  mprintf("\tKeeping %i solvent atoms.\n",NsolventAtoms);
  if (NsolventAtoms < 1) {
    mprintf("Warning: # of solvent atoms to be kept is < 1.\n");
    return Action::SKIP;
  }
  // Store original # of molecules.
  // NOTE: This is stored so that it can be used in the OpenMP section
  //       of action. I dont think iterators are thread-safe.
  NsolventMolecules_ = nSolvent;
  // Create stripped Parm
  if (newParm_!=0) delete newParm_;
  newParm_ = setup.Top().modifyStateByMask(stripMask_);
  if (newParm_==0) {
    mprinterr("Error: Could not create new topology.\n");
    return Action::ERR;
  }
  setup.SetTopology( newParm_ );
  newParm_->Brief("Closest topology:");
  // Allocate space for new frame
  newFrame_.SetupFrameV( setup.Top().Atoms(), setup.CoordInfo() );

  // If prefix given then output stripped parm
  topWriter_.WriteTops( *newParm_ );

# ifdef CUDA
  // Allocate space for simple arrays that will be sent to GPU.
  if (GPU_MEM_ != 0) delete[] GPU_MEM_;
  size_t gpu_mem_size = SolventMols_.size() + // V_distances_
                        (NsolventAtoms * SolventMols_.size() * 3); // V_atom_coords_
  if (!useMaskCenter_) gpu_mem_size += (distanceMask_.Nselected() * 3); // U_atom_coords_
  GPU_MEM_ = new double[ gpu_mem_size ];
  V_atom_coords_ = GPU_MEM_;
  if (useMaskCenter_) {
    V_distances_   = V_atom_coords_ + (NsolventAtoms * SolventMols_.size() * 3);
  } else {
    U_atom_coords_ = V_atom_coords_ + (NsolventAtoms * SolventMols_.size() * 3);
    V_distances_   = U_atom_coords_ + (distanceMask_.Nselected() * 3);
  }
# else
  // Allocate temp space for selected solute atom coords
  if (useMaskCenter_)
    U_cell0_coords_.resize( 3 );
  else
    U_cell0_coords_.resize( distanceMask_.Nselected() * 3 );
# endif

  return Action::MODIFY_TOPOLOGY;
}

// Action_Closest::DoAction()
/** Find the minimum distance between atoms in distanceMask and each 
  * solvent Mask.
  */
Action::RetType Action_Closest::DoAction(int frameNum, ActionFrame& frm) {
  double maxD, Dist2;
  Iarray::const_iterator solvent_atom;

  if (imageOpt_.ImagingEnabled()) {
    // Calculate max possible imaged distance
    maxD = frm.Frm().BoxCrd().Param(Box::X) + frm.Frm().BoxCrd().Param(Box::Y) + 
           frm.Frm().BoxCrd().Param(Box::Z);
    maxD *= maxD;
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );
  } else {
    // If not imaging, set max distance to an arbitrarily large number
    maxD = DBL_MAX;
  }
#ifdef CUDA
// -----------------------------------------------------------------------------
  // Copy solvent atom coords to array
  int NAtoms = SolventMols_[0].solventAtoms.size(); // guaranteed to be same size due to setup
  for (int sMol = 0; sMol < NsolventMolecules_; ++sMol) {
    for(int sAtom = 0; sAtom < NAtoms; sAtom++) {
      int index =  (sAtom * 3 ) + (sMol * 3 * NAtoms);
      const double *a = frm.Frm().XYZ(SolventMols_[sMol].solventAtoms[sAtom]);
      V_atom_coords_[index + 0] = a[0];
      V_atom_coords_[index + 1] = a[1];
      V_atom_coords_[index + 2] = a[2];
    }
  }
  
  if (useMaskCenter_) {
    Vec3 maskCenter = frm.Frm().VGeometricCenter( distanceMask_ );
    Action_Closest_Center( V_atom_coords_, V_distances_, maskCenter.Dptr(),
                           maxD, NsolventMolecules_, NAtoms, imageOpt_.ImagingType(),
                           frm.Frm().BoxCrd().XyzPtr(),
                           frm.Frm().BoxCrd().UnitCell().Dptr(),
                           frm.Frm().BoxCrd().FracCell().Dptr() );
  } else {
    int NSAtoms = distanceMask_.Nselected();
    for (int nsAtom = 0; nsAtom < NSAtoms; ++nsAtom) {
      const double* a = frm.Frm().XYZ(distanceMask_[nsAtom]);
      U_atom_coords_[nsAtom*3 + 0] = a[0];
      U_atom_coords_[nsAtom*3 + 1] = a[1];
      U_atom_coords_[nsAtom*3 + 2] = a[2];
    }

    Action_Closest_NoCenter( V_atom_coords_, V_distances_, U_atom_coords_,
                             maxD, NsolventMolecules_, NAtoms, NSAtoms, imageOpt_.ImagingType(),
                             frm.Frm().BoxCrd().XyzPtr(),
                             frm.Frm().BoxCrd().UnitCell().Dptr(),
                             frm.Frm().BoxCrd().FracCell().Dptr() );
  }
  // Copy distances back into SolventMols_
  for (int sMol = 0; sMol < NsolventMolecules_; sMol++)
    SolventMols_[sMol].D = V_distances_[sMol];

#else /* Not CUDA */
// -----------------------------------------------------------------------------
  if (imageOpt_.ImagingType() == ImageOption::NONORTHO) {
    // ----- NON-ORTHORHOMBIC IMAGING ------------
    // Wrap all solute atoms back into primary cell and save coords
    if (useMaskCenter_) {
      double* uFrac = &U_cell0_coords_[0];
      //  Calc COM and convert to frac coords
      Vec3 center = frm.Frm().BoxCrd().FracCell() * frm.Frm().VGeometricCenter( distanceMask_ );
      // Wrap to primary unit cell
      center[0] = center[0] - floor(center[0]);
      center[1] = center[1] - floor(center[1]);
      center[2] = center[2] - floor(center[2]);
      // Convert back to Cartesian
      frm.Frm().BoxCrd().UnitCell().TransposeMult( uFrac, center.Dptr() );
    } else {
      Image::WrapToCell0( U_cell0_coords_, frm.Frm(), distanceMask_, frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell() );
    }
    // Calculate closest distance of every solvent image to solute
    int mnum;
#   ifdef _OPENMP
#   pragma omp parallel private(mnum, Dist2, solvent_atom)
    {
#   pragma omp for
#   endif
    for (mnum = 0; mnum < NsolventMolecules_; mnum++)
    {
      MolDist& Mol = SolventMols_[mnum];
      Mol.D = maxD;
      for (solvent_atom = Mol.solventAtoms.begin();
           solvent_atom != Mol.solventAtoms.end(); ++solvent_atom)
      {
        // Convert to frac coords
        Vec3 vFrac = frm.Frm().BoxCrd().FracCell() * Vec3( frm.Frm().XYZ( *solvent_atom ) );
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
              Vec3 vCart = frm.Frm().BoxCrd().UnitCell().TransposeMult( vFrac + Vec3(ix, iy, iz) );
              // Loop over all solute atoms
              for (unsigned int idx = 0; idx < U_cell0_coords_.size(); idx += 3)
              {
                double x = vCart[0] - U_cell0_coords_[idx  ];
                double y = vCart[1] - U_cell0_coords_[idx+1];
                double z = vCart[2] - U_cell0_coords_[idx+2];
                Dist2 = x*x + y*y + z*z;
                Mol.D = std::min(Dist2, Mol.D);
              } // END loop over solute atoms
            } // END loop over images (Z)
      } // END loop over solvent mol atoms
    } // END loop over solvent molecules
#   ifdef _OPENMP
    } /* END pragma omp parallel */
#   endif
  } else {
    // ----- ORTHORHOMBIC/NO IMAGING -------------
    if (useMaskCenter_) {
      //  Calc COM
      Vec3 center = frm.Frm().VGeometricCenter( distanceMask_ );
      U_cell0_coords_[0] = center[0];
      U_cell0_coords_[1] = center[1];
      U_cell0_coords_[2] = center[2];
    } else {
      // Store selected solute coordinates.
      int idx = 0;
      for (AtomMask::const_iterator atm = distanceMask_.begin();
                                    atm != distanceMask_.end(); ++atm, idx += 3)
      {
        const double* XYZ = frm.Frm().XYZ( *atm );
        U_cell0_coords_[idx  ] = XYZ[0];
        U_cell0_coords_[idx+1] = XYZ[1];
        U_cell0_coords_[idx+2] = XYZ[2];
      }
    }
    // Loop over all solvent molecules
    int mnum;
#   ifdef _OPENMP
#   pragma omp parallel private(mnum, Dist2, solvent_atom)
    {
#   pragma omp for
#   endif
    for (mnum = 0; mnum < NsolventMolecules_; mnum++)
    {
      MolDist& Mol = SolventMols_[mnum];
      Mol.D = maxD;
      for (solvent_atom = Mol.solventAtoms.begin();
           solvent_atom != Mol.solventAtoms.end(); ++solvent_atom)
      {
        Vec3 Vcoord( frm.Frm().XYZ( *solvent_atom ) );
        // Loop over all solute atoms
        for (unsigned int idx = 0; idx < U_cell0_coords_.size(); idx += 3)
        {
          Vec3 Ucoord( U_cell0_coords_[idx], U_cell0_coords_[idx+1], U_cell0_coords_[idx+2] );
          if (imageOpt_.ImagingType() == ImageOption::ORTHO)
            Dist2 = DIST2_ImageOrtho( Vcoord, Ucoord, frm.Frm().BoxCrd() );
          else
            Dist2 = DIST2_NoImage( Vcoord, Ucoord );
          Mol.D = std::min(Dist2, Mol.D);
        } // END loop over solute atoms
      } // END loop over solvent mol atoms
    } // END loop over solvent molecules
#   ifdef _OPENMP
    } /* END pragma omp parallel */
#   endif
  }
#endif /* CUDA */
// -----------------------------------------------------------------------------

  // Sort distances
  std::sort( SolventMols_.begin(), SolventMols_.end(), moldist_cmp() );
  // Add first closestWaters solvent atoms to stripMask
  std::vector<MolDist>::const_iterator solventend = SolventMols_.begin() + closestWaters_;
  Iarray::const_iterator katom = keptWaterAtomNum_.begin();
  for ( std::vector<MolDist>::const_iterator solvent = SolventMols_.begin();
                                             solvent != solventend;
                                           ++solvent, ++katom ) 
  {
    //mprintf("DEBUG:\tmol %i ",(*solvent).mol);
    //solvent->mask.PrintMaskAtoms("Mask");
    stripMask_.AddMaskAtPosition( solvent->mask, *katom );
    // Record which water molecules are closest if requested
    if (outFile_!=0) {
      int fnum = frm.TrajoutNum() + 1;
      framedata_->Add(Nclosest_, &fnum);
      moldata_->Add(Nclosest_, &(solvent->mol));
      Dist2 = sqrt( solvent->D );
      distdata_->Add(Nclosest_, &Dist2);
      solvent_atom = solvent->mask.begin();
      int solvent_first_atom = *solvent_atom + 1; 
      atomdata_->Add(Nclosest_, &solvent_first_atom);
      ++Nclosest_;
    }
    // DEBUG - print first closestWaters distances
    //mprintf("DEBUG: Mol %i   D2= %lf   Atom0= %i\n",(*it).mol, (*it).D, (*it).mask->Selected[0]);
  }

  // Modify and set frame
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("action_stripMask");
  newFrame_.SetFrame(frm.Frm(), stripMask_);
  frm.SetFrame( &newFrame_ );

  return Action::MODIFY_COORDS;
}

#ifdef MPI
/** Since datasets are actually # frames * closestWaters_ in length, need to
  * sync here.
  */
int Action_Closest::SyncAction() {
  if (outFile_ == 0) return 0;
  // Get total number of closest entries.
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &Nclosest_, 1, MPI_INT, &(rank_frames[0]) );
  for (int rank = 1; rank < trajComm_.Size(); rank++)
    Nclosest_ += rank_frames[ rank ];
  framedata_->Sync( Nclosest_, rank_frames, trajComm_ );
  framedata_->SetNeedsSync( false );
  moldata_->Sync( Nclosest_, rank_frames, trajComm_ );
  moldata_->SetNeedsSync( false );
  distdata_->Sync( Nclosest_, rank_frames, trajComm_ );
  distdata_->SetNeedsSync( false );
  atomdata_->Sync( Nclosest_, rank_frames, trajComm_ );
  atomdata_->SetNeedsSync( false );
  return 0;
}
#endif
