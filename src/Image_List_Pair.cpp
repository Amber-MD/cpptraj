#include "Image_List_Pair.h"
#include "AtomMask.h"
#include "Topology.h"
#include "CpptrajStdio.h"

/** Set up list to contain all residues selected by mask expression */
int Image::List_Pair::SetupList(Topology const& Parm, std::string const& maskExpression)
{
  AtomMask mask;
  if (mask.SetMaskString(maskExpression)) return 1;
  if ( Parm.SetupIntegerMask( mask ) ) return 1;
  mask.MaskInfo();
  if (mask.None()) return -1;
  std::vector<int> resnums = Parm.ResnumsSelectedBy( mask );
  for (std::vector<int>::const_iterator it = resnums.begin(); it != resnums.end(); ++it)
  {
    begin_.push_back( Parm.Res(*it).FirstAtom() );
    end_.push_back( Parm.Res(*it).LastAtom() );
  }
  return 0;
}

/** \return Unit containing all residues. */
Unit Image::List_Pair::AllEntities() const {
  Unit unitOut;
  for (unsigned int jdx = 0; jdx != begin_.size(); ++jdx)
    for (int idx = begin_[jdx]; idx != end_[jdx]; ++idx)
      unitOut.AddIndex( idx );
  return unitOut;
}

/** Print residues to be imaged to STDOUT. */
void Image::List_Pair::PrintEntities() const {
  mprintf("\tImaging residues defined by atom ranges:\n");
  for (unsigned int jdx = 0; jdx != begin_.size(); ++jdx)
    mprintf("\t  %10i to %10i\n", begin_[jdx]+1, end_[jdx]);
}
