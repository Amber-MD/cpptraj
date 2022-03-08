#ifndef INC_CENTROID_COORD_H
#define INC_CENTROID_COORD_H
#include "Centroid_Coord.h"
#include "../CpptrajFile.h"

void Cpptraj::Cluster::Centroid_Coord::Print(std::string const& fnameIn) const {
  // Write Amber coordinate format
  CpptrajFile cOut;
  FileName fname(fnameIn);
  if (cOut.OpenWrite(fname)) return;
  cOut.Printf("%-80s\n", fname.base());
  int col = 0;
  for (int ic = 0; ic != cframe_.size(); ic++) {
    cOut.Printf("%8.3f", cframe_[ic]);
    ++col;
    if (col == 10) {
      cOut.Printf("\n");
      col = 0;
    }
  }
  if (col > 0) cOut.Printf("\n");
}

#endif
