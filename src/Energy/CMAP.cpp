#include "CMAP.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
CMAP::CMAP() {}

/** Calculate CMAP energy */
double CMAP::Ene_CMAP(CmapArray const& Cmaps, Frame const& frameIn) {
  double ene_cmap = 0.0;

  for (CmapArray::const_iterator cmap = Cmaps.begin();
                                 cmap != Cmaps.end(); ++cmap
  {
    const double* ixyz = frameIn.XYZ( cmap->A1() );
    const double* jxyz = frameIn.XYZ( cmap->A2() );
    const double* kxyz = frameIn.XYZ( cmap->A3() );
    const double* lxyz = frameIn.XYZ( cmap->A4() );
    const double* mxyz = frameIn.XYZ( cmap->A5() );

  }
  return ene_cmap;
}
