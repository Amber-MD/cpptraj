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
    
