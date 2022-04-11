#include "Image_List_Mask.h"
#include "Topology.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

/** Set up list to include atoms selected by mask expression. */
int Image::List_Mask::SetupList(Topology const& Parm, std::string const& maskExpression)
{
  mask_.ResetMask();
  if (mask_.SetMaskString(maskExpression)) return 1;
  if ( Parm.SetupIntegerMask( mask_ ) ) return 1;
  mask_.MaskInfo();
  if (mask_.None()) return -1;
  return 0;
}

/** \return Unit containing all atoms. */
Unit Image::List_Mask::AllEntities() const {
  Unit unitOut;
  for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); ++at)
    unitOut.AddIndex( *at );
  return unitOut;
}

/** Print atoms to be imaged to STDOUT. */
void Image::List_Mask::PrintEntities() const {
  std::string range = ArrayToRangeExpression( mask_.Selected(), 1 );
  mprintf("\tImaging the following atoms: %s\n", range.c_str());
}
