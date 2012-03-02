// CoordList
#include "CoordList.h"
#include <cstddef>

// CONSTRUCTOR
CoordList::CoordList() {
  ncoords = 0;
}

// DESTRUCTOR
CoordList::~CoordList() {
  for (int crd = 0; crd < ncoords; crd++)
    delete[] coordList[crd];
}

// CoordList::AddCoordsByMask()
/** Given an array of double coordinates and a corresponding atom mask,
  * store the selected atoms in a float array.
  * Currently only for use with Integer masks.
  */
int CoordList::AddCoordsByMask(double *Xin, AtomMask *MaskIn) {
  if (Xin==NULL || MaskIn==NULL) return 1;
  //if (MaskIn->Selected==NULL) return 1;
  if (MaskIn->None()) return 1;
  // Allocate space for coords
  int natom = MaskIn->Nselected;
  float *coord = new float[ natom * 3 ];
  int catom3 = 0;
  for (int maskidx = 0; maskidx < natom; maskidx++) {
    int xatom3 = MaskIn->Selected[maskidx] * 3;
    coord[catom3  ] = (float)Xin[xatom3  ];
    coord[catom3+1] = (float)Xin[xatom3+1];
    coord[catom3+2] = (float)Xin[xatom3+2];
    catom3+=3;
  }
  coordList.push_back( coord );
  natomList.push_back( natom );
  ncoords++;
  return 0;
}

// CoordList::Coord()
/** Return a pointer to the given coord. Set corresponding number of atoms in
  * natom.
  */
float *CoordList::Coord(int crd, int *natom) {
  if (crd<0 || crd >= ncoords) return NULL;
  *natom = natomList[crd];
  return coordList[crd];
}

// CoordList::operator[]()
/** Return a pointer to the given coord. For use when natom is already known.
  */
float *CoordList::operator[](int crd) {
  if (crd<0 || crd >= ncoords) return NULL;
  return coordList[crd];
}

