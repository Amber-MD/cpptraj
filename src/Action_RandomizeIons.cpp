#include <cmath> // sqrt
#include "Action_RandomizeIons.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_RandomizeIons::Action_RandomizeIons() :
  overlap_(0.0),
  min_(0.0),
  debug_(0)
{}

void Action_RandomizeIons::Help() const {
  mprintf("\t<mask> [around <mask> by <distance>] [overlap <value>]\n"
          "\t[noimage] [seed <value>]\n"
          "  Swap positions of ions in <mask> with randomly selected solvent molecules.\n");
}

// Action_RandomizeIons::Init()
Action::RetType Action_RandomizeIons::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get first mask
  std::string ionmask = actionArgs.GetMaskNext();
  if (ionmask.empty()) {
    mprinterr("Error: randomizeions: No mask for ions specified.\n");
    return Action::ERR;
  }
  if (ions_.SetMaskString( ionmask )) return Action::ERR;

  // Get Keywords
  image_.InitImaging( !actionArgs.hasKey("noimage") );
  int seed = actionArgs.getKeyInt("seed", -1);
  overlap_ = actionArgs.getKeyDouble("overlap", 3.5);
  min_ = actionArgs.getKeyDouble("by", 3.5);
  // Pre-square overlap and min
  overlap_ *= overlap_;
  min_ *= min_;
  // If no around mask specified, leave blank
  std::string aroundmask = actionArgs.GetStringKey("around");
  if (!aroundmask.empty()) {
    if (around_.SetMaskString( aroundmask )) return Action::ERR;
  }
  
  // INFO
  mprintf("    RANDOMIZEIONS: Swapping postions of ions in mask '%s' with solvent.\n",
          ions_.MaskString());
  mprintf("\tNo ion can get closer than %.2f angstroms to another ion.\n", sqrt( overlap_ ));
  if (around_.MaskStringSet())
    mprintf("\tNo ion can get closer than %.2f angstroms to atoms in mask '%s'\n",
            sqrt( min_ ), around_.MaskString());
  if (!image_.UseImage())
    mprintf("\tImaging of the coordinates will not be performed.\n");
  if (seed > 0)
    mprintf("\tRandom number generator seed is %i\n", seed);
  RN_.rn_set( seed );

  return Action::OK;
}

// Action_RandomizeIons::Setup()
Action::RetType Action_RandomizeIons::Setup(ActionSetup& setup) {
  if (setup.Top().Nsolvent() < 1) {
    mprinterr("Error: This command only works if solvent information has been specified.\n");
    return Action::ERR;
  }

  // Set up ion mask
  if (setup.Top().SetupIntegerMask( ions_ )) return Action::ERR;
  if ( ions_.None() ) {
    mprintf("Warning: Ion mask '%s' has no atoms.\n", ions_.MaskString());
    return Action::SKIP;
  }
  mprintf("\tIon mask is '%s' (%i atoms)\n", ions_.MaskString(), ions_.Nselected());

  // Set up the around mask if necessary
  if (around_.MaskStringSet()) {
    if (setup.Top().SetupIntegerMask( around_ )) return Action::ERR;
    if ( around_.None() ) {
      mprintf("Warning: Around mask '%s' has no atoms.\n", around_.MaskString());
    } else {
      mprintf("\tAround mask is '%s' (%i atoms)\n", around_.MaskString(),
              around_.Nselected());
    }
  }

  // Check that each ion is only a single atom residue.
  // NOTE: Should this be a molecule check instead? If so can then get rid of ResSize 
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    int res = setup.Top()[*ion].ResNum();
    if (debug_ > 0)
      mprintf("\tAtom %i is in residue %i which is %i atoms\n",
              *ion+1, res+1, setup.Top().Res( res ).NumAtoms() );
    if ( setup.Top().Res( res ).NumAtoms() > 1 ) {
      mprintf("Warning: randomizeions: Ion atom %i belongs to residue %i which\n",
              *ion + 1, res + 1);
      mprintf("Warning:                contains more than 1 atom (%i)!\n", 
              setup.Top().Res( res ).NumAtoms());
    }
  }

  // Save each solvent molecule.
  solvMols_.clear();
  solvMols_.reserve( setup.Top().Nsolvent() );
  for (Topology::mol_iterator Mol = setup.Top().MolStart();
                              Mol != setup.Top().MolEnd(); ++Mol)
  {
    if ( Mol->IsSolvent() ) {
      solvMols_.push_back( Mol->MolUnit() ); 
    }
  }
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  // Allocate solvent molecule considered for swap mask
  solvent_.resize( solvMols_.size() );

  return Action::OK;
}

// Action_RandomizeIons::DoAction()
Action::RetType Action_RandomizeIons::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 ucell, recip;

  if (image_.ImageType() == NONORTHO)
    frm.Frm().BoxCrd().ToRecip(ucell, recip);
  // Loop over all solvent molecules and mark those that are too close to the solute
  int n_active_solvent = 0;
  for (unsigned int idx = 0; idx != solvMols_.size(); idx++) {
    solvent_[idx] = true;
    if (around_.MaskStringSet()) {
      // Is solvent molecule too close to any atom in the around mask?
      // For speed, use only the first solvent atom.
      const double* solventXYZ = frm.Frm().XYZ( solvMols_[idx].Front() );
      for (AtomMask::const_iterator atom = around_.begin(); atom != around_.end(); ++atom)
      {
        double dist = DIST2( solventXYZ, frm.Frm().XYZ(*atom), image_.ImageType(),
                             frm.Frm().BoxCrd(), ucell, recip);
        if (dist < min_) {
          solvent_[idx] = false;
          //mprintf("RANDOMIZEIONS: water %i only %.2f ang from around @%i\n",
          //        smolnum, sqrt(dist), *atom+1);
          break;
        }
      }
    }
    if (solvent_[idx]) ++n_active_solvent;
  }
  if (n_active_solvent < ions_.Nselected()) {
    mprinterr("Error: Fewer active solvent molecules (%i) than ions (%i)\n",
              n_active_solvent, ions_.Nselected());
    return Action::ERR;
  }

  // DEBUG - print solvent molecule mask
  if (debug_ > 2) {
    mprintf("RANDOMIZEIONS: The following waters (first atoms) are ACTIVE so far:\n");
    unsigned int smoltot = 0;
    for (unsigned int idx = 0; idx != solvMols_.size(); idx++) {
      if (solvent_[idx]) {
        mprintf(" %6i ", solvMols_[idx].Front()+1 );
        ++smoltot;
        if (smoltot%10 == 0) mprintf("\n");
      }
    }
    mprintf("RANDOMIZEIONS: A total of %u waters (out of %zu) are active\n",
            smoltot, solvent_.size());
  }

  // Outer loop over all ions
  for (AtomMask::const_iterator ion1 = ions_.begin(); ion1 != ions_.end(); ++ion1)
  {
    //mprintf("RANDOMIZEIONS: Processing ion atom %i\n", *ion1+1);
    // Is a potential solvent molecule close to any of the ions (except this one)?
    for (unsigned int idx = 0; idx != solvMols_.size(); idx++) {
      if (solvent_[idx]) {
        // For speed, use only the first solvent atom.
        const double* solventXYZ = frm.Frm().XYZ( solvMols_[idx].Front() );
        // This solvent is active; check distance to all other ions
        for (AtomMask::const_iterator ion2 = ions_.begin(); ion2 != ions_.end(); ++ion2)
        {
          if (*ion1 != *ion2) {
            double dist = DIST2( solventXYZ, frm.Frm().XYZ(*ion2), image_.ImageType(),
                                 frm.Frm().BoxCrd(), ucell, recip);
            if (dist < overlap_) {
              // This solvent mol is too close to another ion.
              solvent_[idx] = false;
              //mprintf("RANDOMIZEIONS: water %i only %.2f ang from ion @%i\n",
              //        smolnum, sqrt(dist), *ion2+1);
              break;
            }
          }
        } // END inner loop over ions
      }
    } // END loop over solvent molecules
    // The solvent_ array should now be true for all solvent molecules eligible
    // to swap with ion. Select an elegible solvent molecule.
    int loop = 1;
    int swapMol = 0;
    while (loop > 0 && loop < (int)solvMols_.size()) {
      // Run the random number generator so that the same number is not produced
      // when the seed was set manually.
      swapMol = (int)(RN_.rn_gen() * (double)solvMols_.size());
      if ( solvent_[swapMol] ) 
        loop = -1;
      else 
        ++loop;
    }

    // If a suitable solvent molecule was found, swap it.
    if (loop > 0) {
      mprintf("Warning: Tried to swap ion @%i with %i random waters\n",*ion1+1,loop);
      mprintf("Warning: and couldn't meet criteria; skipping.\n");
    } else {
      if (debug_ > 2)
        mprintf("RANDOMIZEIONS: Swapping solvent mol %i for ion @%i\n", swapMol+1, *ion1+1);
      const double* ionXYZ = frm.Frm().XYZ( *ion1 );
      // Use only the first solvent atom
      const double* watXYZ = frm.Frm().XYZ( solvMols_[swapMol].Front() );
      Vec3 trans( ionXYZ[0] - watXYZ[0],
                  ionXYZ[1] - watXYZ[1],
                  ionXYZ[2] - watXYZ[2]);
      frm.ModifyFrm().Translate( trans, solvMols_[swapMol] );
      trans.Neg();
      frm.ModifyFrm().Translate( trans, *ion1 );
    }
  } // END outer loop over all ions

  return Action::MODIFY_COORDS;
}
