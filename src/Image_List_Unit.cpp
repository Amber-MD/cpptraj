#include "Image_List_Unit.h"
#include "AtomMask.h"
#include "Topology.h"

int Image::List_Unit::SetupList(Topology const& Parm, std::string const& maskExpression)
{
  AtomMask mask;
  if (mask.SetMaskString(maskExpression)) return 1;
  if ( Parm.SetupIntegerMask( mask ) ) return 1;
  mask.MaskInfo();
  if (mask.None()) return -1;
  std::vector<int> molnums = Parm.MolnumsSelectedBy( mask );
  for (std::vector<int>::const_iterator it = molnums.begin(); it != molnums.end(); ++it)
  {
    units_.push_back( Parm.Mol(*it).MolUnit() );
  }
  return 0;
}
