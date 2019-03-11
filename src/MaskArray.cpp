#include "MaskArray.h"
#include "Topology.h"
#include "CpptrajStdio.h"

int Cpptraj::MaskArray::SetMaskExpression(std::string const& exprIn, SelectType typeIn)
{
  // Set the mask expression first.
  mainMask_.ResetMask();
  if (mainMask_.SetMaskString( exprIn )) {
    mprinterr("Error: '%s' is an invalid mask string.\n", exprIn.c_str());
    return 1;
  }
  type_ = typeIn;
  return 0;
}

/** Set up masks. */
int Cpptraj::MaskArray::SetupMasks(Topology const& topIn)
{
  if (type_ == BY_MOLECULE && topIn.Nmol() < 1) {
    mprintf("Warning: '%s' has no molecule information, cannot setup by molecule.\n",
             topIn.c_str());
    return 1;
  }
  if (topIn.SetupIntegerMask( mainMask_ )) return 1;
  masks_.clear();
  int last = -1;
  int current = 0;
  for (AtomMask::const_iterator atm = mainMask_.begin(); atm != mainMask_.end(); ++atm)
  {
    switch (type_) {
      case BY_ATOM     : current = *atm; break;
      case BY_RESIDUE  : current = topIn[*atm].ResNum(); break;
      case BY_MOLECULE : current = topIn[*atm].MolNum(); break;
    }
    if (current != last) {
      masks_.push_back( AtomMask() );
      masks_.back().SetNatoms( topIn.Natom() );
    }
    masks_.back().AddSelectedAtom( *atm );
    last = current;
  }

  // DEBUG
  mprintf("DEBUG: %zu masks created\n", masks_.size());
  for (Marray::const_iterator it = masks_.begin(); it != masks_.end(); ++it)
  {
    mprintf("  %6u :", it-masks_.begin());
    for (AtomMask::const_iterator atm = it->begin(); atm != it->end(); ++atm)
      mprintf(" %i", *atm + 1);
    mprintf("\n");
  }

  return 0;
}
