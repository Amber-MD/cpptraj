#include "MaskArray.h"
#include "AtomMask.h"
#include "Topology.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
Cpptraj::MaskArray::MaskArray() :
  type_(BY_ATOM),
  maxAtomsPerMask_(0),
  sameNumAtomsPerMask_(false)
{}

/** Save max atoms per mask, check if # atoms is the same as last time. */
void Cpptraj::MaskArray::checkAtomsPerMask( int nselected ) {
  if (maxAtomsPerMask_ == 0) {
    maxAtomsPerMask_ = nselected;
    sameNumAtomsPerMask_ = true;
  } else {
    if (sameNumAtomsPerMask_) {
      if (nselected != maxAtomsPerMask_)
        sameNumAtomsPerMask_ = false;
    }
    if (nselected > maxAtomsPerMask_)
      maxAtomsPerMask_ = nselected;
  }
}

/** Set up masks. */
int Cpptraj::MaskArray::SetupMasks(AtomMask const& maskIn, Topology const& topIn)
{
  if (type_ == BY_MOLECULE && topIn.Nmol() < 1) {
    mprintf("Warning: '%s' has no molecule information, cannot setup by molecule.\n",
             topIn.c_str());
    return 1;
  }
  masks_.clear();
  if ( maskIn.None() ) {
    mprintf("Warning: Nothing selected by mask '%s'\n", maskIn.MaskString());
    return 0;
  }
  int last = -1;
  int current = 0;
  maxAtomsPerMask_ = 0;
  sameNumAtomsPerMask_ = true;
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm)
  {
    switch (type_) {
      case BY_ATOM     : current = *atm; break;
      case BY_RESIDUE  : current = topIn[*atm].ResNum(); break;
      case BY_MOLECULE : current = topIn[*atm].MolNum(); break;
    }
    if (current != last) {
      if (!masks_.empty())
        checkAtomsPerMask( masks_.back().Nselected() );
      masks_.push_back( AtomMask() );
      masks_.back().SetNatoms( topIn.Natom() );
    }
    masks_.back().AddSelectedAtom( *atm );
    last = current;
  }
  if (!masks_.empty())
    checkAtomsPerMask( masks_.back().Nselected() );

  return 0;
}

void Cpptraj::MaskArray::Debug() const {
  // DEBUG
  mprintf("DEBUG: %zu masks created (max atoms=%i, same=%i)\n", masks_.size(),
          maxAtomsPerMask_, (int)sameNumAtomsPerMask_);
  for (Marray::const_iterator it = masks_.begin(); it != masks_.end(); ++it)
  {
    mprintf("  %6u :", it-masks_.begin());
    for (AtomMask::const_iterator atm = it->begin(); atm != it->end(); ++atm)
      mprintf(" %i", *atm + 1);
    mprintf("\n");
  }
}
