#include <cmath>
#include "PairList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

PairList::PairList() {}

/** This leads to cellNeighbor_ dimensions of 7x10 */
const int PairList::cellOffset_ = 3;

// PairList::InitPairList()
int PairList::InitPairList(double cutIn, double skinNBin) {
  std::fill(translateVec_, translateVec_+18, Vec3(0.0));
  if (Fill_CellNeighbor()) return 1;
  cutList_ = cutIn + skinNBin;
  mprintf("DEBUG: cutList= %12.5f\n", cutList_);
  nGridX_0_ = -1;
  nGridY_0_ = -1;
  nGridZ_0_ = -1;
  maxNptrs_ = ((2*cellOffset_ + 1) * (2*cellOffset_ + 1) + 1 ) / 2;
  mprintf("DEBUG: max number of pointers= %i\n", maxNptrs_);
  return 0;
}

/** Set up the cellNeighbor_ array.
  * The neighbor cells of a cell of interest (call it A) that are only
  * "forward" of that cell reside in the plane of the cell, and three
  * planes "forward" in the z direction.
  *
  *       A = cell for which we want to get neighbor list
  *       x = forward neighbor cell within 3 cells
  *       o = same as x except this cell has same x index as A
  *
  *            i3          i3+1            i3+2         i3+3
  *        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  *        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  *        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  * ---->       Axxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  *                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  *                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  *                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
  *
  * A cell and its neighbors span the x direction over 7 cells (3 on each 
  * side and the cell iself) cellNeighbor_ array contains a 0, 1, or -1 
  * for whether the neighbor cell has x index outside the unit cell and
  * needs to be translated along the x uc vector positive one cell, negative,
  * or not at all.
  * There are 10 cases of neighbor cells in x direction, cellNeighbor_(*,1)
  * (0000000). All have x indices within the unit cell. Cases 2,3,4 are 
  * special for the row of neighbors containing the cell A itself (see
  * arrow). This is the only case where neighbors with x index are not
  * included in the list search since those  cells are "before" the cell of
  * interest, and this is a "look forward only" method. So for this row of
  * cells, only 4 cellNeighbor_ values are needed: cell A, and the three
  * cells to right.
  * The cases represent whether the neighbors extend out of the unit cell by
  * one, two, or three cells. Entry 1 is for cell A and must be 0 since it
  * must be in the unit cell. (last 4 entries are ignored for this set).
  *          (*,2)  ( 0 0 0 1 0 0 0)
  *          (*,3)  ( 0 0 1 1 0 0 0)
  *          (*,4)  ( 0 1 1 1 0 0 0)
  * Cases 5,6,7 are for neighbors that extend to the left out of the UC.
  *          (*,5)  (-1 0 0 0 0 0 0)
  *          (*,6)  (-1-1 0 0 0 0 0)
  *          (*,7)  (-1-1-1 0 0 0 0)
  * Cases 8,9,10 are same as 2,3,4 except that there are 7 cells in all other
  * rows:
  *          (*,8)  ( 0 0 0 0 0 0 1)
  *          (*,9)  ( 0 0 0 0 0 1 1)
  *          (*,10) ( 0 0 0 0 1 1 1)
  */
int PairList::Fill_CellNeighbor() {
  // Sanity check. Currently wired for 3 cells in forward direction.
  if (cellOffset_ != 3) {
    mprinterr("Internal Error: PairList::Fill_CellNeighbor():\n"
              "Internal Error: Should be only 3 cells but currently %i\n",
              cellOffset_);
    return 1;
  }
  // Most cells will not be translated, so set xtran to 0 for all
  // possibilities, then fill in the 1 and -1 entries for neighbor
  // cells that are beyond UC edges.
  // CASE 1: Zero out array
  for (int i = 0; i < 2*cellOffset_+1; i++)
    for (int j = 0; j < cellOffset_*3+1; j++)
      cellNeighbor_[i][j] = 0;
  // CASES 2,3,4
  for (int j = 0; j < cellOffset_; j++)
    for (int i = cellOffset_-j; i < cellOffset_+1; i++)
      cellNeighbor_[i][j+1] = 1;
  // CASES 5,6,7
  for (int j = 0; j < cellOffset_; j++)
    for (int i = 0; i <= j; i++)
      cellNeighbor_[i][cellOffset_+j+1] = -1;
  // CASES 8,9,10
  for (int j = 0; j < cellOffset_; j++)
    for (int i = 0; i < j+1; i++) 
      cellNeighbor_[2*cellOffset_-i][2*cellOffset_+1+j] = 1;

//  for (int i = 0; i < 7; i++)
//    for (int j = 0; j < 10; j++)
//      mprintf("XTRAN %3i%3i%3i\n", i+1, j+1, cellNeighbor_[i][j]);
  for (int j = 0; j < 10; j++) {
    mprintf("XTRAN %3i", cellNeighbor_[0][j]);
    for (int i = 1; i < 7; i++)
      mprintf("%3i", cellNeighbor_[i][j]);
    mprintf("\n");
  }
  return 0;
}

/// Round floating point to nearest whole number.
static inline double ANINT(double xIn) {
  double fpart, ipart;
  fpart = modf(xIn, &ipart);
  if (fpart < 0.0) fpart = -fpart;
  if (fpart < 0.5)
    return ipart;
  if (xIn > 0.0)
    return ipart + 1.0;
  else
    return ipart - 1.0;
}

// PairList::MapCoords()
void PairList::MapCoords(Frame const& frmIn, Matrix_3x3 const& ucell,
                         Matrix_3x3 const& recip, AtomMask const& maskIn)
{
  t_map_.Start();
  Frac_.clear();
  Frac_.reserve( maskIn.Nselected() );
  Image_.clear();
  Image_.reserve( maskIn.Nselected() );

  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
  {
    Vec3 fc = recip * Vec3(frmIn.XYZ(*atom));
    // TODO: Use ANINT below to have frac coords/grid that matches sander
    //       for now. Should eventually just use the floor() call to bound
    //       between 0.0 and 1.0
    //Frac_.push_back( Vec3(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2])) );
    // Wrap back into primary cell between -.5 and .5
    Frac_.push_back( Vec3(fc[0]-ANINT(fc[0]), fc[1]-ANINT(fc[1]), fc[2]-ANINT(fc[2])) );
    //mprintf("FRAC %10.5f%10.5f%10.5f\n", Frac_.back()[0], Frac_.back()[1], Frac_.back()[2]);
    Image_.push_back( ucell.TransposeMult( Frac_.back() ) );
  }
  mprintf("DEBUG: Mapped coords for %zu atoms.\n", Frac_.size());
  // Allocate memory
  atomCell_.resize( Frac_.size() );
  atomGridIdx_.resize( Frac_.size() );
  t_map_.Stop();
}

// PairList::FillTranslateVec()
void PairList::FillTranslateVec(Matrix_3x3 const& ucell) {
  int iv = 0;
  for (int i3 = 0; i3 < 2; i3++)
    for (int i2 = -1; i2 < 2; i2++)
      for (int i1 = -1; i1 < 2; i1++)
        translateVec_[iv++] = ucell.TransposeMult( Vec3(i1, i2, i3) );
  for (int i = 0; i < 18; i++)
    mprintf("TRANVEC %3i%12.5f%12.5f%12.5f\n", i+1, translateVec_[i][0],
            translateVec_[i][1], translateVec_[i][2]);
}

// PairList::SetupGrids()
int PairList::SetupGrids(Vec3 const& recipLengths) {

  int nghb0 = cellOffset_;
  int nghb1 = nghb0;
  int nghb2 = nghb0;
  int nghb3 = nghb0;

  double dc1 = cutList_ / (double)nghb1;
  double dc2 = cutList_ / (double)nghb2;
  double dc3 = cutList_ / (double)nghb3;

  nGridX_ = std::max(1, (int)(recipLengths[0] / dc1)); // nucgrd1
  nGridY_ = std::max(1, (int)(recipLengths[1] / dc2));
  nGridZ_ = std::max(1, (int)(recipLengths[2] / dc3));

  // TODO Add non-periodic case
  // Check short range cutoff
  dc1 = recipLengths[0] / (double)nGridX_;
  dc2 = recipLengths[1] / (double)nGridY_;
  dc3 = recipLengths[2] / (double)nGridZ_;
  double cut = (double)nghb1 * dc1;
  if (nghb2*dc2 < cut)
    cut = (double)nghb2*dc2;
  if (nghb3*dc3 < cut)
    cut = (double)nghb3*dc3;
  //if(nogrdptrs)cut=cutlist
  // DEBUG
  mprintf("Number of grids per unit cell in each dimension: %i %i %i\n",
          nGridX_, nGridY_, nGridZ_);
  //mprintf("Unit cell edge lengths in each dimension: %9.3f %9.3f %9.3f\n",);
  mprintf("Distance between parallel faces of unit cell: %9.3f %9.3f %9.3f\n",
          recipLengths[0], recipLengths[1], recipLengths[2]);
  mprintf("Distance between faces of short range grid subcells: %9.3f %9.3f %9.3f\n",
          dc1, dc2, dc3);
  mprintf("Resulting cutoff from subcell neighborhoods is %f\n", cut);
  if (cut < cutList_) {
    mprinterr("Error: Resulting cutoff %f too small for lower limit %f\n", cut, cutList_);
    return 1;
  }
  // Allocation
  nGridMax_ = nGridX_ * nGridY_ * nGridZ_;
  nLoGrid_.resize( nGridMax_ );
  nHiGrid_.resize( nGridMax_ );
  nAtomsInGrid_.resize( nGridMax_ );
  idxOffset_.resize( nGridMax_ );
  myGrids_.resize( nGridMax_ );
  neighborPtr_.resize( nGridMax_ );
  neighborTrans_.resize( nGridMax_ );
  for (int i = 0; i != nGridMax_; i++) {
    neighborPtr_[i].resize( maxNptrs_ );
    neighborTrans_[i].resize( maxNptrs_ );
  }
  size_t nbrSize = (size_t)(nGridMax_ * maxNptrs_) * sizeof(int);
  mprintf("DEBUG: Grid memory total: %s\n", 
          ByteString(((nLoGrid_.size() + nHiGrid_.size() + myGrids_.size() +
                       nAtomsInGrid_.size() + idxOffset_.size()) * sizeof(int)) +
                     (2 * nbrSize) +
                     ((myGrids_.size()) * sizeof(bool)),
                     BYTE_DECIMAL).c_str());
  return 0;
}

/** Grid mapped atoms in the unit cell into grid subcells according
  * to the fractional coordinates.
  */
void PairList::GridUnitCell() {
  nAtomsInGrid_.assign( nGridMax_, 0 ); // TODO do in setup?
  //TODO for non-periodic
  double shift = 0.5;
  // Find out which grid subcell each atom is in.
  for (unsigned int i = 0; i != Frac_.size(); i++) {
    Vec3 const& frac = Frac_[i];
    int i1 = (int)((frac[0] + shift) * (double)nGridX_);
    int i2 = (int)((frac[1] + shift) * (double)nGridY_);
    int i3 = (int)((frac[2] + shift) * (double)nGridZ_);
    int idx = (i3*nGridX_*nGridY_)+(i2*nGridX_)+i1;
    atomCell_[i] = idx;
//    mprintf("GRID atom assigned to cell %6i%6i%10.5f%10.5f%10.5f\n", i+1, idx+1,
//            frac[0],frac[1],frac[2]);
    if (idx < 0 || idx >= nGridMax_) { // Sanity check
      mprinterr("Internal Error: Grid is out of range (>= %i || < 0)\n", nGridMax_);
      return;
    }
    nAtomsInGrid_[idx]++;
  }

  // Find the offset of the starting atoms for each grid subcell.
  idxOffset_[0] = 0;
  for (int i = 1; i < nGridMax_; i++) {
    idxOffset_[i] = idxOffset_[i-1] + nAtomsInGrid_[i-1];
//    mprintf("INDOFF %6i\n", idxOffset_[i]);
    // Reset atom count in each cell.
    nAtomsInGrid_[i-1] = 0;
  }
  nAtomsInGrid_[nGridMax_-1] = 0;

  // Get list of atoms sorted by grid cell so that atoms in first subcell
  // are first, atoms in second subcell come after that, etc.
  for (unsigned int i = 0; i != Frac_.size(); i++) {
    int idx = atomCell_[i];
    int j = nAtomsInGrid_[idx] + idxOffset_[idx];
    nAtomsInGrid_[idx]++;
    atomGridIdx_[j] = (int)i;
  }
  for (unsigned int j = 0; j != atomGridIdx_.size(); j++)
    mprintf("INDATG %6i\n", atomGridIdx_[j]+1);
}

/** Determine list of nearby cell line centers that need to be searched 
  * for pair list.
  * Example: 2d grid of cells, X is cell of interest. Here we are using 3
  *          cells to the cutoff, so need 3 cells distant in each direction.
  *          We only point to the center of each string of 7 cells in the
  *          x direction, and the pair list routine knows to check this cell
  *          and three on each side of it in the +x and -x directions.
  *          This routine will list the center of cell strings:
  *         A                          B
  *    ..................    ..................   and three more like layer B
  *    ..................    ..................
  *    .....ooo#ooo......    .....ooo#ooo......
  *    .....ooo#ooo......    .....ooo#ooo......
  *    .....ooo#ooo......    .....ooo#ooo......
  *    ........Xooo......    .....ooo#ooo......
  *    ..................    .....ooo#ooo......
  *    ..................    .....ooo#ooo......
  *    ..................    .....ooo#ooo......
  *    ..................    ..................
  *         cell X and its      Next layer ahead
  *         neighbors o
  *         center cell #
  *         (ahead only search)
  *
  * The pointer list contains the identity of all the # cells and the X cell.
  */
void PairList::GridPointers(int myindexlo, int myindexhi) {
  for (std::vector<Iarray>::iterator it = neighborPtr_.begin();
                                     it != neighborPtr_.end(); ++it)
    it->clear();
  int index0_min = 999999999; // TODO FIXME
  int index0_max = -1;
  int nGridXY = nGridX_ * nGridY_;
  int nghb1 = cellOffset_;
  int nghb2 = cellOffset_;
  int nghb3 = cellOffset_;
  for (int i3 = 0; i3 < nGridZ_; i3++)
  {
    int i3grd = nGridXY * i3;
    for (int i2 = 0; i2 < nGridY_; i2++)
    {
      int i2grd = i3grd + nGridX_ * i2;
      for (int i1 = 0; i1 < nGridX_; i1++)
      {
        int idx = i2grd + i1;
        if (idx >= myindexlo && idx < myindexhi)
        {
          myGrids_[idx] = 1;
          int index0 = idx - myindexlo;
          index0_max = std::max(index0_max, index0);
          index0_min = std::min(index0_min, index0);
          // Setup xtran index for this cell
          int xtindex_1, xtindex_2;
          if (i1 + nghb1 >= nGridX_) {
            xtindex_1 = 1 + i1 + nghb1 - nGridX_;
            xtindex_2 = xtindex_1 + 2*nghb1;
          } else if (i1 - nghb1 < 0) {
            xtindex_1 = 0;
            xtindex_2 = 2*nghb1 + 2 - i1;
          } else {
            xtindex_1 = 0;
            xtindex_2 = 0;
          }
          // First do the plane that this cell is in.
          // First row of neighbors starts with the cell itself (Cell0)
          int num = 0;
          neighborPtr_[index0][num] = idx;
          // No y or z translate, this cell has to be in unit cell.
          neighborTrans_[index0][num] = 5 + (xtindex_1 << 8) ;
          for (int jj2 = 0; jj2 < nghb1; jj2++)
            myGrids_[ idx + jj2 - nGridX_*cellNeighbor_[jj2][xtindex_1] ] = 1;
          // Next rows of neighbors are the nghb2 rows of 2*nghb1+1 cells
          // above. Use the center cell of the row as the pointer (with
          // the same x position as Cell0). This way the pointer cell
          // must be in the unit cell w.r.t. the x direction.
          num = 1;
          for (int j2 = 0; j2 < nghb2; j2++, num++) {
            int j2grd = idx + j2*nGridX_;
            if (j2 + i2 <= nGridY_)
              neighborTrans_[index0][num] = 5 + (xtindex_2 << 8);
            else {
              j2grd -= nGridXY;
              neighborTrans_[index0][num] = 8 + (xtindex_2 << 8);
            }
            neighborPtr_[index0][num] = j2grd;
            for (int jj2 = -nghb1; jj2 <= nghb1; jj2++)
              myGrids_[ j2grd + jj2 - nGridX_*cellNeighbor_[jj2+nghb1][xtindex_2] ] = 1;
          }
          // Now there are nghb3 planes of (2*nghb1+1)(2*nghb2+1) cells
          // ahead that are good neighbors.
          for (int j3 = 0; j3 < nghb3; j3++) {
            int j3grd = j3 * nGridXY;
            int k3off;
            if (j3 + i3 >= nGridZ_) {
              j3grd = j3grd - nGridXY * nGridZ_;
              k3off = 9;
            } else
              k3off = 0;
            for (int j2 = -nghb2; j2 <= nghb2; nghb2++, num++) {
              int j2grd = idx + j3grd + j2*nGridX_;
              if(j2+i2 >= nGridY_) {
                neighborTrans_[index0][num] = 8 + k3off + (xtindex_2 << 8);
                j2grd -= nGridXY;
              } else if (j2+i2 < 0) {
                neighborTrans_[index0][num] = 2 + k3off + (xtindex_2 << 8);
                j2grd = j2grd + nGridXY;
              } else
                neighborTrans_[index0][num] = 5 + k3off + (xtindex_2 << 8);
              neighborPtr_[index0][num] = j2grd;
              for (int jj2=-nghb1; jj2 <= nghb1; jj2++) {
                myGrids_[ j2grd + jj2 - nGridX_*cellNeighbor_[jj2+nghb1][xtindex_2] ] = 1;
              }
            }
          }
        } // end if grid index is mine
      } // end for i1
    } // end for i2
  } // end for i3
  for (int i1 = 0; i1 < nGridMax_; i1++)
    for (int i2 = 0; i2 < maxNptrs_; i2++)
      mprintf("NGHBPTR %6i%6i%6i\n", i1, i2, neighborPtr_[i1][i2]);
}

// PairList::CreatePairList()
int PairList::CreatePairList(Frame const& frmIn, AtomMask const& maskIn) {
  Matrix_3x3 ucell, recip;
  frmIn.BoxCrd().ToRecip(ucell, recip);
  MapCoords(frmIn, ucell, recip, maskIn);

  FillTranslateVec(ucell);

  if (SetupGrids(frmIn.BoxCrd().RecipLengths(recip))) return 1;

  // Maybe for parallelization later.
  int myindexlo = 0;
  int myindexhi = nGridMax_;

  GridUnitCell();

  // Only calculate grid pointers if number of grids has changed.
  if (nGridX_ != nGridX_0_ ||
      nGridY_ != nGridY_0_ ||
      nGridZ_ != nGridZ_0_)
  {
    nGridX_0_ = nGridX_;
    nGridY_0_ = nGridY_;
    nGridZ_0_ = nGridZ_;
    GridPointers(myindexlo, myindexhi);
  }

  return 0;
}
