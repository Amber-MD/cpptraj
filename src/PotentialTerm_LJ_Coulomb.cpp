#include "PotentialTerm_LJ_Coulomb.h"
#include "Topology.h"
#include "CharMask.h"

/**  CONSTRUCTOR */
PotentialTerm_LJ_Coulomb::PotentialTerm_LJ_Coulomb() :
  PotentialTerm(SIMPLE_LJ_Q),
  nonbond_(0),
  Evdw_(0),
  Eelec_(0)
{}

/** Setup nonbonds. */
int PotentialTerm_LJ_Coulomb::SetupTerm(Topology const& topIn, CharMask const& maskIn,
                                  EnergyArray& Earray)
{
  selectedAtoms_.clear();
  nonSelectedAtoms_.clear();
  typeIndices_.clear();

  selectedAtoms_.reserve( maskIn.Nselected() );
  nonSelectedAtoms_.reserve( topIn.Natom() - maskIn.Nselected() );
  for (int i = 0; i < topIn.Natom(); i++) {
    if (maskIn.AtomInCharMask( i ))
      selectedAtoms_.push_back( i );
    else
      nonSelectedAtoms_.push_back( i );
    typeIndices_.push_back( topIn[i].TypeIndex() );
  }

  nonbond_ = &(topIn.Nonbond());

  return 0;
}
