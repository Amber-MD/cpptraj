#include <cmath> // floor
#include <algorithm> // min and max
#include "PairList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString()

PairList::PairList() :
  cutList_(0.0),
  debug_(0),
  nGridX_(-1),
  nGridY_(-1),
  nGridZ_(-1),
  nGridX_0_(-1),
  nGridY_0_(-1),
  nGridZ_0_(-1)
{}

/** This leads to cellNeighbor_ dimensions of 7x10 */
const int PairList::cellOffset_ = 3;

// PairList::InitPairList()
int PairList::InitPairList(double cutIn, double skinNBin, int debugIn) {
  debug_ = debugIn;
  std::fill(translateVec_, translateVec_+18, Vec3(0.0));
  //if (Fill_CellNeighbor()) return 1;
  cutList_ = cutIn + skinNBin;
  nGridX_0_ = -1;
  nGridY_0_ = -1;
  nGridZ_0_ = -1;
  //maxNptrs_ = ((2*cellOffset_ + 1) * (2*cellOffset_ + 1) + 1 ) / 2;
  //mprintf("DEBUG: max number of pointers= %i\n", maxNptrs_);
  return 0;
}

int PairList::SetupPairList(Box const& boxIn) {
  Matrix_3x3 ucell, recip;
  boxIn.ToRecip(ucell, recip);
  return SetupPairList( boxIn.Type(), boxIn.RecipLengths(recip) );
}

// PairList::SetupPairList()
int PairList::SetupPairList(Box::BoxType typeIn, Vec3 const& recipLengthsIn) {
  Timer t_setup;
  t_setup.Start();
  if (typeIn == Box::NOBOX) {
    mprinterr("Error: Pair list code currently requires box coordinates.\n");
    return 1;
  }

  // Allocate/reallocate memory
  if (SetupGrids(recipLengthsIn)) return 1;
  t_setup.Stop();
  t_setup.WriteTiming(1, "Pair List Setup:");
  mprintf("\tGrid dimensions: %i %i %i (%zu total).\n", nGridX_, nGridY_, nGridZ_, cells_.size());
  return 0;
}

// PairList::FillTranslateVec()
/** Fill the translate vector array with offset values based on this
  * unit cell. Only need forward direction, so no -Z.
  */
void PairList::FillTranslateVec(Matrix_3x3 const& ucell) {
  int iv = 0;
  for (int i3 = 0; i3 < 2; i3++)
    for (int i2 = -1; i2 < 2; i2++)
      for (int i1 = -1; i1 < 2; i1++)
        translateVec_[iv++] = ucell.TransposeMult( Vec3(i1, i2, i3) );
  //for (int i = 0; i < 18; i++)
  //  mprintf("TRANVEC %3i%12.5f%12.5f%12.5f\n", i+1, translateVec_[i][0],
  //          translateVec_[i][1], translateVec_[i][2]);
}

// PairList::CreatePairList()
int PairList::CreatePairList(Frame const& frmIn, Matrix_3x3 const& ucell,
                             Matrix_3x3 const& recip, AtomMask const& maskIn)
{
  t_total_.Start();
  // Calculate translation vectors based on current unit cell.
  FillTranslateVec(ucell);
  // If box size has changed a lot this will reallocate grid
  t_gridpointers_.Start();
  if (SetupGrids(frmIn.BoxCrd().RecipLengths(recip))) return 1;
  t_gridpointers_.Stop();
  // Place atoms in grid cells
  t_map_.Start();
  GridUnitCell(frmIn, ucell, recip, maskIn);
  t_map_.Stop();
  t_total_.Stop();
  return 0;
}

// PairList::GridAtom()
void PairList::GridAtom(int atomIdx, Vec3 const& frac, Vec3 const& cart) {
  int i1 = (int)((frac[0]) * (double)nGridX_);
  int i2 = (int)((frac[1]) * (double)nGridY_);
  int i3 = (int)((frac[2]) * (double)nGridZ_);
  int idx = (i3*nGridX_*nGridY_)+(i2*nGridX_)+i1;
  //mprintf("GRID2 atom assigned to cell %6i%6i%10.5f%10.5f%10.5f\n",
  //        atomIdx+1, idx+1, frac[0], frac[1], frac[2]);
  if (idx < 0 || idx >= (int)cells_.size()) { // Sanity check
    mprinterr("Internal Error: Grid %i is out of range (>= %zu || < 0)\n",
              idx, cells_.size());
    return;
  }
  cells_[idx].AddAtom( AtmType(atomIdx, cart) );
  Frac_.push_back( frac );
}

/** Place selected atoms into grid cells. Convert to fractional coords, wrap
  * into primary cell, then determine grid cell.
  */
void PairList::GridUnitCell(Frame const& frmIn, Matrix_3x3 const& ucell,
                             Matrix_3x3 const& recip, AtomMask const& maskIn)
{
  // Clear any existing atoms in cells.
  for (Carray::iterator cell = cells_.begin(); cell != cells_.end(); ++cell)
    cell->ClearAtoms();
  Frac_.clear();
  Frac_.reserve( maskIn.Nselected() );
  if (frmIn.BoxCrd().Type() == Box::ORTHO) {
    // Orthogonal imaging
    for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
    {
      const double* XYZ = frmIn.XYZ(*atom);
      Vec3 fc( XYZ[0]*recip[0],    XYZ[1]*recip[4],    XYZ[2]*recip[8]   );
      Vec3 fcw(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2]));
      Vec3 ccw(fcw[0]*ucell[0],    fcw[1]*ucell[4],    fcw[2]*ucell[8]   );
      GridAtom( atom-maskIn.begin(), fcw, ccw );
    }
  } else {
    // Non-orthogonal imaging
    for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
    {
      Vec3 fc = recip * Vec3(frmIn.XYZ(*atom));
      Vec3 fcw(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2]));
      Vec3 ccw = ucell.TransposeMult( fcw );
      GridAtom( atom-maskIn.begin(), fcw, ccw );
    }
  }
}

// PairList::SetupGrids()
/** Determine grid sizes. If this is the first time this routine is called
  * or the grid sizes have changed (re)allocate the grid and set up the
  * grid pointers.
  */
int PairList::SetupGrids(Vec3 const& recipLengths) {
  int offsetX = cellOffset_;
  int offsetY = cellOffset_;
  int offsetZ = cellOffset_;

  double dc1 = cutList_ / (double)offsetX;
  double dc2 = cutList_ / (double)offsetY;
  double dc3 = cutList_ / (double)offsetZ;

  nGridX_ = std::max(1, (int)(recipLengths[0] / dc1)); // nucgrd1
  nGridY_ = std::max(1, (int)(recipLengths[1] / dc2));
  nGridZ_ = std::max(1, (int)(recipLengths[2] / dc3));

  // Check if grid (re)allocation needs to happen
  if ( nGridX_ == nGridX_0_ &&
       nGridY_ == nGridY_0_ &&
       nGridZ_ == nGridZ_0_ )
    return 0;
  if (nGridX_0_ != -1) // -1 is the initial allocation
    mprintf("Warning: Unit cell size has changed so much that grid must be recalculated.\n"
            "Warning: Old sizes= {%i, %i, %i}  New sizes= {%i, %i, %i}\n",
            nGridX_0_, nGridY_0_, nGridZ_0_, nGridX_, nGridY_, nGridZ_);
  nGridX_0_ = nGridX_;
  nGridY_0_ = nGridY_;
  nGridZ_0_ = nGridZ_;

  // TODO Add non-periodic case
  // Check short range cutoff
  dc1 = recipLengths[0] / (double)nGridX_;
  dc2 = recipLengths[1] / (double)nGridY_;
  dc3 = recipLengths[2] / (double)nGridZ_;
  double cut = (double)offsetX * dc1;
  if (offsetY*dc2 < cut)
    cut = (double)offsetY*dc2;
  if (offsetZ*dc3 < cut)
    cut = (double)offsetZ*dc3;
  //if(nogrdptrs)cut=cutlist
  // Allocation
  int nGridMax = nGridX_ * nGridY_ * nGridZ_;
  cells_.clear();
  cells_.resize( nGridMax );
  if (debug_ > 0) {
    mprintf("DEBUG: Number of grids per unit cell in each dimension: %i %i %i\n",
            nGridX_, nGridY_, nGridZ_);
    //mprintf("Unit cell edge lengths in each dimension: %9.3f %9.3f %9.3f\n",);
    mprintf("DEBUG: Distance between parallel faces of unit cell: %9.3f %9.3f %9.3f\n",
            recipLengths[0], recipLengths[1], recipLengths[2]);
    mprintf("DEBUG: Distance between faces of short range grid subcells: %9.3f %9.3f %9.3f\n",
            dc1, dc2, dc3);
    mprintf("DEBUG: Resulting cutoff from subcell neighborhoods is %f\n", cut);
    mprintf("%zu total grid cells\n", cells_.size());
  }
  if (cut < cutList_) {
    mprinterr("Error: Resulting cutoff %f too small for lower limit %f\n", cut, cutList_);
    return 1;
  }

  // NOTE: myindex are for parallelization later (maybe).
  int myindexlo = 0;
  int myindexhi = (int)cells_.size();
  CalcGridPointers(myindexlo, myindexhi);

  PrintMemory();

  return 0;
}

static inline void CheckOffset(int nGrid, int& offset, char dir) {
  if ((nGrid+1) / 2 <= offset) {
    offset = std::max((nGrid / 2) - 1, 1);
    mprintf("Warning: %c cell offset reset to %i\n", dir, offset);
  }
}

/** Calculate all of the forward neighbors for each grid cell and if
  * needed which translate vector is appropriate.
  * +X: need self and offsetX cells to the "right".
  * +Y: need 2*offsetX+1 cells in "above" offsetY rows.
  * +Z: need 2*offsetX+1 by 2*offsetY+1 cells in "above" offsetZ rows.
  * The translate vector is calculated based on offset indices. The x
  * and y offset indices are actually +1 so we can calculate an index
  * via idx = oz*3*3 + oy*3 + ox. oz is either 0 or 1 since we only
  * care about forward direction.
  */
void PairList::CalcGridPointers(int myindexlo, int myindexhi) {
  //Matrix<bool> PairCalcd;
  //PairCalcd.resize(nGridMax_, 0); // Half matrix
  //std::fill(PairCalcd.begin(), PairCalcd.end(), false);
  //  int idx = (i3*nGridX_*nGridY_)+(i2*nGridX_)+i1;
  int offsetX = cellOffset_;
  int offsetY = cellOffset_;
  int offsetZ = cellOffset_;
  // Check the cell offsets. If they are too big comared to the grid
  // size we will end up calculating too many values.
  CheckOffset(nGridX_, offsetX, 'X');
  CheckOffset(nGridY_, offsetY, 'Y');
  CheckOffset(nGridZ_, offsetZ, 'Z');

  int NP_ = 0;
  int nGridXY = nGridX_ * nGridY_;
  for (int nz = 0; nz != nGridZ_; nz++)
  {
    int idx3 = nz * nGridXY;
    for (int ny = 0; ny != nGridY_; ny++)
    {
      int idx2 = idx3 + (ny*nGridX_);
      for (int nx = 0; nx != nGridX_; nx++)
      {
        int idx = idx2 + nx; // Absolute grid cell index
        if (idx >= myindexlo && idx < myindexhi) {
          Iarray& Nbr = cells_[idx].neighborPtr_;
          Iarray& Ntr = cells_[idx].neighborTrans_;
          //int NP = 0;
//          mprintf("DBG: Cell %3i%3i%3i (%i):", nx,ny,nz, idx);
          // Get this cell and all cells ahead in the X direction.
          // This cell is always a "neighbor" of itself.
          int maxX = nx + offsetX + 1;
          for (int ix = nx; ix < maxX; ix++, NP_++) {
            // Wrap ix if necessary
            if (ix < nGridX_) {
//              mprintf(" %i+0", idx2+ix);
              Nbr.push_back( idx2 + ix );
              Ntr.push_back( 4 ); // No translation. 0 0 0
            } else {
//              mprintf(" %i+1", idx2+ix - nGridX_);
              Nbr.push_back( idx2 + ix - nGridX_ );
              Ntr.push_back( 5 ); // Translate by +1 in X. 1 0 0
            }
          }
          // Get all cells in the Y+ direction
          int minX = nx - offsetX;
          int minY = ny + 1;
          int maxY = ny + offsetY + 1;
          for (int iy = minY; iy < maxY; iy++) {
            // Wrap iy if necessary
            int wy, oy;
            if (iy < nGridY_) {
              wy = iy;
              oy = 1; // 0
            } else {
              wy = iy - nGridY_;
              oy = 2; // 1
            }
            int jdx2 = idx3 + (wy*nGridX_);
            for (int ix = minX; ix < maxX; ix++, NP_++) {
              // Wrap ix if necessary
              int wx, ox;
              if (ix < 0) {
                wx = ix + nGridX_;
                ox = 0; // -1
              } else if (ix < nGridX_) {
                wx = ix;
                ox = 1; // 0
              } else {
                wx = ix - nGridX_;
                ox = 2; // 1
              }
              // Calc new index
              int jdx = jdx2 + wx; // Absolute neighbor grid cell index
//              mprintf(" %i%+i%+i", jdx, ox-1, oy-1);
              Nbr.push_back( jdx );
              int tidx = oy*3 + ox;
              Ntr.push_back( tidx );
            }
          }
          // Get all cells in the +Z direction
          int minZ = nz + 1;
          int maxZ = nz + offsetZ + 1;
          for (int iz = minZ; iz < maxZ; iz++) {
            // Wrap iz if necessary
            int wz, oz;
            if (iz < nGridZ_) {
              wz = iz;
              oz = 0; // 0
            } else {
              wz = iz - nGridZ_;
              oz = 1; // 1
            }
            int jdx3 = wz * nGridXY;
            minY = ny - offsetY;
            for (int iy = minY; iy < maxY; iy++) {
              // Wrap iy if necessary
              int wy, oy;
              if (iy < 0) {
                wy = iy + nGridY_;
                oy = 0; // -1
              } else if (iy < nGridY_) {
                wy = iy;
                oy = 1; // 0
              } else {
                wy = iy - nGridY_;
                oy = 2; // 1
              }
              int jdx2 = jdx3 + (wy*nGridX_);
              for (int ix = minX; ix < maxX; ix++, NP_++) {
                // Wrap ix if necessary
                int wx, ox;
                if (ix < 0) {
                  wx = ix + nGridX_;
                  ox = 0; // -1
                } else if (ix < nGridX_) {
                  wx = ix;
                  ox = 1; // 0
                } else {
                  wx = ix - nGridX_;
                  ox = 2; // 1
                }
                // Calc new index
                int jdx = jdx2 + wx; // Absolute neighbor grid cell index
//                mprintf(" %i%+i%+i%+i", jdx, ox-1, oy-1, oz);
                Nbr.push_back( jdx );
                int tidx = oz*9 + oy*3 + ox;
                Ntr.push_back( tidx );
              }
            }
          }

//          mprintf("\n");
          //if (NP > maxNptrs_) {
          //  mprinterr("Error: You overflowed! (%i)\n", NP);
          //  return;
          //}
//          mprintf("Grid %3i%3i%3i (%i): %zu neighbors\n", nx,ny,nz, idx, Nbr.size());
//          for (unsigned int i = 0; i != Nbr.size(); i++)
//            mprintf("\t%i [%i]\n", Nbr[i], Ntr[i]);
          //mprintf("\n");
/*          for (unsigned int i = 0; i != Nbr.size(); i++) {
            if (PairCalcd.element(idx,Nbr[i]))
              mprintf("Warning: Interaction %i %i already calcd.\n", idx, Nbr[i]);
            else
              PairCalcd.setElement(idx,Nbr[i], true);
          }*/
        } // END my cell
      } // nx
    } // ny
  } // nz
}

void PairList::Timing(double total) const {
  t_total_.WriteTiming(2, "Pair List: ", total);
  t_map_.WriteTiming(3,          "Map Coords:      ", t_total_.Total());
  t_gridpointers_.WriteTiming( 3,"Recalc Grid Ptrs:", t_total_.Total());
}

void PairList::PrintMemory() const {
  size_t total = 0;
  for (Carray::const_iterator cell = cells_.begin(); cell != cells_.end(); ++cell)
    total += cell->MemSize();
  total += ((Frac_.size() * sizeof(Vec3)) + sizeof(Varray));
  mprintf("\tTotal Grid memory: %s\n", ByteString(total, BYTE_DECIMAL).c_str());
}
