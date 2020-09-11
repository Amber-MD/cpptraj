#include "Image_List_Mask.h"
#include "Topology.h"

int Image::List_Mask::SetupList(Topology const& Parm, std::string const& maskExpression)
{
  mask_.ResetMask();
  if (mask_.SetMaskString(maskExpression)) return 1;
  if ( Parm.SetupIntegerMask( mask_ ) ) return 1;
  mask_.MaskInfo();
  if (mask_.None()) return -1;
  return 0;
}
