#include <cmath> // sqrt
#include <cstdlib> // random, srandom
#include "Action_RandomizeIons.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_RandomizeIons::Action_RandomizeIons() :
  overlap_(0),
  min_(0),
  aroundmask_(0),
  seed_(1)
{
  // Imaging on by default.
  useImage_ = true;
}

/** randomizeions <mask> [around <mask> by <distance>] [overlap <value>] 
  *               [noimage] [seed <value>]
  */
int Action_RandomizeIons::init() {
  // Get first mask
  char* ionmask = actionArgs.getNextMask();
  if (ionmask == NULL) {
    mprinterr("Error: randomizeions: No mask for ions specified.\n");
    return 1;
  }
  ions_.SetMaskString( ionmask );

  // Get Keywords
  if (actionArgs.hasKey("noimage"))
    useImage_ = false;
  seed_ = actionArgs.getKeyInt("seed", -1);
  overlap_ = actionArgs.getKeyDouble("overlap", 3.5);
  min_ = actionArgs.getKeyDouble("by", 3.5);
  // Pre-square overlap and min
  overlap_ *= overlap_;
  min_ *= min_;
  // If no around mask specified, leave blank
  aroundmask_ = actionArgs.getKeyString("around", 0);
  if (aroundmask_!=0)
    around_.SetMaskString( aroundmask_ );
  
  // INFO
  mprintf("    RANDOMIZEIONS: swapping the postions of the ions in mask [%s]\n", 
          ions_.MaskString());
  mprintf("                   with the solvent. No ions can get closer than %5.2f angstroms\n",
          sqrt( overlap_ ));
  mprintf("                   to another ion.\n");
  if (aroundmask_!=0) {
    mprintf("                   No ion can get closer than %5.2f angstroms to mask [%s]\n",
            sqrt( min_ ), around_.MaskString());
  }
  if (!useImage_)
    mprintf("                   Imaging of the coordinates will not be performed.\n");
  if (seed_ > 0) {
    mprintf("                   Random number generator seed is %i\n", seed_);
    srandom((unsigned) seed_);
  }

  return 0;
}

int Action_RandomizeIons::setup() {
  if (currentParm->Nsolvent()==0) {
    mprinterr("Warning: randomizeions: This command only works if solvent information\n");
    mprinterr("Warning:                has been specified.");
    //mprinterr(" See the \"solvent\" command.");
    mprinterr("\n");
    return 1;
  }

  // Set up ion mask
  if (currentParm->SetupIntegerMask( ions_ )) return 1;
  if ( ions_.None() ) {
    mprinterr("Warning: randomizeions: Mask [%s] has no atoms.\n", ions_.MaskString());
    return 1;
  }
  mprintf("\tIon mask is [%s] (%i atoms)\n", ions_.MaskString(), ions_.Nselected());

  // Check that each ion is only a single atom residue.
  // NOTE: Should this be a molecule check instead? If so can then get rid of ResSize 
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    int res = (*currentParm)[*ion].ResNum();
    if ( currentParm->ResSize( res ) > 1 ) {
      mprintf("Warning: randomizeions: Ion atom %i belongs to residue %i which\n",
              *ion + 1, res + 1);
      mprintf("Warning:                contains more than 1 atom (%i)!\n", 
              currentParm->ResSize( res ));
    }
  }

  // Check the solvent information to make sure that each solvent listed has the
  // same number of atoms in each molecule; otherwise a uniform trajectory is not
  // possible and therefore this command will be ignored.
  Topology::mol_iterator solvmol = currentParm->SolventStart();
  int NsolventAtoms = (*solvmol).NumAtoms();
  ++solvmol;
  for (; solvmol != currentParm->SolventEnd(); solvmol++) {
    if ( NsolventAtoms != (*solvmol).NumAtoms() ) {
      mprinterr("Warning: randomizeions: Solvent molecules in %s are not of uniform size.\n",
                currentParm->c_str());
      mprinterr("       First solvent mol = %i atoms, this solvent mol = %i atoms.\n",
                NsolventAtoms, (*solvmol).NumAtoms());
      return 1;
    }
  }

  // Allocate solvent molecule mask
  solvent_.resize( currentParm->Nsolvent() );

  return 0;
}

int Action_RandomizeIons::action() {
  double ucell[9], recip[9], trans[3], watXYZ[3], ionXYZ[3];

  if (imageType_==Frame::NONORTHO)
    currentFrame->BoxToRecip(ucell, recip);
  // loop over all solvent molecules and mark those that are too close to the solute
  std::vector<bool>::iterator smask = solvent_.begin();
  for (Topology::mol_iterator solvmol = currentParm->SolventStart();
                              solvmol != currentParm->SolventEnd(); solvmol++)
  {
    *smask = true;
    // is solvent molecule to near any atom in the around mask?
    if (aroundmask_ != 0) {
      for (AtomMask::const_iterator atom = around_.begin(); atom != around_.end(); ++atom)
      {
        double dist = currentFrame->DIST2( (*solvmol).BeginAtom(), *atom, imageType_, 
                                           ucell, recip);
        if (dist < min_) {
          *smask = false;
          mprintf("RANDOMIZEIONS: water @%i only %.2f ang from around @%i\n",
                  (*solvmol).BeginAtom()+1, sqrt(dist), *atom+1);
          break;
        }
      }
    }
    ++smask;
  }

  // Outer loop over all ions
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    mprintf("RANDOMIZEIONS: Processing ion atom %i\n", *ion+1);
    // is a potential solvent molecule close to any of the ions (except this one)?
    smask = solvent_.begin();
    for (Topology::mol_iterator solvmol = currentParm->SolventStart();
                                solvmol != currentParm->SolventEnd(); solvmol++)
    {
      if (*smask) {
        // if this solvent is active, check distance to all other ions
        for (AtomMask::const_iterator ion2 = ions_.begin(); ion2 != ions_.end(); ++ion2)
        {
          if (*ion != *ion2) {
            double dist = currentFrame->DIST2( (*solvmol).BeginAtom(), *ion2, imageType_,
                                               ucell, recip);
            if (dist < overlap_) {
              *smask = false;
              mprintf("RANDOMIZEIONS: water @%i only %.2f ang from ion @%i\n",
                      (*solvmol).BeginAtom()+1, sqrt(dist), *ion2+1);
              break;
            }
          }
        } // END inner loop over ions
      }
      ++smask;
    } // END loop over solvent molecules

    // solvent should now be true for all solvent molecules eligible to
    // swap with ion.
    int loop = 1;
    int swapMol = 0;
    while (loop > 0 && loop < 10000) {
      // Run the random number generator so that the same number is not produced
      // when the seed was set manually.
      swapMol = random() % currentParm->Nsolvent();
      if ( solvent_[swapMol] ) 
        loop = -1;
      else 
        ++loop;
    }

    // If a suitable solvent molecule was found, swap it.
    if (loop > 0) {
      mprintf("Warning: randomizeions: Tried to swap ion @%i with %i random waters\n",*ion+1,loop);
      mprintf("Warning:                and couldn't meet criteria; skipping.\n");
    } else {
      mprintf("RANDOMIZEIONS: Swapping solvent mol %i for ion @%i\n", swapMol+1, *ion+1);
      currentFrame->GetAtomXYZ( ionXYZ, *ion );
      Topology::mol_iterator solvmol = currentParm->SolventStart() + swapMol;
      currentFrame->GetAtomXYZ( watXYZ, (*solvmol).BeginAtom() );
      trans[0] = ionXYZ[0] - watXYZ[0];
      trans[1] = ionXYZ[1] - watXYZ[1];
      trans[2] = ionXYZ[2] - watXYZ[2];
      currentFrame->Translate( trans, (*solvmol).BeginAtom(), (*solvmol).EndAtom() );
      trans[0] = -trans[0];
      trans[1] = -trans[1];
      trans[2] = -trans[2];
      currentFrame->Translate( trans, *ion, *ion );
    }

  } // END outer loop over all ions

  return 0;
}
