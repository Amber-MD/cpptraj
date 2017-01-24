#include "PairList.h"
#include "CpptrajStdio.h"

PairList::PairList() {}

// PairList::InitPairList()
int PairList::InitPairList() {
  std::fill(translateVec_, translateVec_+18, Vec3(0.0));
  if (Fill_CellNeighbor()) return 1;
  return 0;
}

/** This leads to cellNeighbor_ dimensions of 7x10 */
const int PairList::cellOffset_ = 3;

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
    // Wrap back into primary cell
    //Frac_.push_back( Vec3(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2])) );
    Frac_.push_back( Vec3(fc[0]-(int)fc[0], fc[1]-(int)fc[1], fc[2]-(int)fc[2]) );
    Image_.push_back( ucell.TransposeMult( Frac_.back() ) );
  }
  mprintf("DEBUG: Mapped coords for %zu atoms.\n", Frac_.size());
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

// PairList::CreatePairList()
void PairList::CreatePairList(Frame const& frmIn, AtomMask const& maskIn) {
  Matrix_3x3 ucell, recip;
  frmIn.BoxCrd().ToRecip(ucell, recip);
  MapCoords(frmIn, ucell, recip, maskIn);

  FillTranslateVec(ucell);

}
