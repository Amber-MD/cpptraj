#include "Image_List_Pair.h"
#include "AtomMask.h"
#include "Topology.h"

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
