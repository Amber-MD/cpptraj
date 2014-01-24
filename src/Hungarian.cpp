#include <cfloat> // DBL_MAX
#include <stdexcept>
#include "Hungarian.h"
#include "Constants.h" // SMALL
#include "CpptrajStdio.h"

// CONSTRUCTOR: FIXME: Could get expensive with exceptions, get rid of
Hungarian::Hungarian(DataSet_MatrixDbl const& mIn) :
  matrix_(mIn),
  lineThroughRow_(mIn.Nrows(), false),
  lineThroughCol_(mIn.Ncols(), false),
  assignRowToCol_(mIn.Ncols(), -1),
  assignColToRow_(mIn.Nrows(), -1),
  nrows_((int)mIn.Nrows()),
  ncols_((int)mIn.Ncols())
{
  if (matrix_.Kind() != DataSet_2D::FULL) {
    mprinterr("Internal Error: Hungarian: Requires full double-precision matrix.\n");
    throw(std::bad_alloc());
  }
  if (matrix_.Nrows() != matrix_.Ncols()) {
    mprinterr("Internal Error: Hungarian: Non-square matrix not yet supported.");
    throw(std::bad_alloc());
  }
}

/** Initialize matrix for Hungarian algorithm. **/
int Hungarian::Initialize(size_t Ncols) {
  if (matrix_.Allocate2D( Ncols, Ncols )) return 1;
  lineThroughRow_.assign(matrix_.Nrows(), false);
  lineThroughCol_.assign(matrix_.Ncols(), false);
  assignRowToCol_.assign(matrix_.Ncols(), -1);
  assignColToRow_.assign(matrix_.Nrows(), -1);
  nrows_ = (int)matrix_.Nrows();
  ncols_ = (int)matrix_.Ncols();
  return 0;
}

// Hungarian::AssignRowsToColumns()
/** Assign as many rows to columns as possible. This is done by iteratively
  * finding the row/col with the lowest number of zeros, then assigning that 
  * to the first unassigned col/row which has a zero. This is repeated until
  * no more assignments can be made.
  */
int Hungarian::AssignRowsToColumns() {
  int Nassigned = 0;
  // For this stage, lines will be drawn through ASSIGNED rows/cols.
  lineThroughRow_.assign(matrix_.Nrows(), false);
  lineThroughCol_.assign(matrix_.Ncols(), false);
  assignRowToCol_.assign(matrix_.Ncols(), -1);
  assignColToRow_.assign(matrix_.Ncols(), -1);
  int assigned = 1;
  while (assigned > 0) { 
    assigned = 0;
    // Find row with lowest number of zeros
    int minRow = -1;
    int minRowZeros = matrix_.Nrows() + 1;
    for (int row = 0; row < nrows_; ++row) {
      if (!lineThroughRow_[row]) {
        int elt = row * ncols_;
        int nzeros = 0;
        for (int col = 0; col < ncols_; ++col, ++elt) {
          if (!lineThroughCol_[col] && matrix_[elt] < Constants::SMALL) ++nzeros;
        }
        if (nzeros != 0 && nzeros < minRowZeros) {
          minRowZeros = nzeros;
          minRow = row;
        }
      }
    }
    // Find col with lowest number of zeros
    int minCol = -1;
    int minColZeros = ncols_ + 1;
    for (int col = 0; col < ncols_; ++col) {
      if (!lineThroughCol_[col]) {
        int elt = col;
        int nzeros = 0;
        for (int row = 0; row < nrows_; ++row, elt += ncols_) {
          if (!lineThroughRow_[row] && matrix_[elt] < Constants::SMALL) ++nzeros;
        }
        if (nzeros != 0 && nzeros < minColZeros) {
          minColZeros = nzeros;
          minCol = col;
        }
      }
    }
#   ifdef DEBUG_HUNGARIAN
    mprintf("Min Row %i (%i), min Column %i (%i)\n", minRow, minRowZeros,
            minCol, minColZeros);
#   endif
    if (minRow == -1 && minCol == -1) // No zeros left
      break;
    if (minColZeros < minRowZeros) {
#     ifdef DEBUG_HUNGARIAN 
      mprintf("\tColumn %i has min # of zeros (%i)\n", minCol, minColZeros);
#     endif
      // Assign column to the first unassigned row whose elt is zero.
      int elt = minCol; 
      for (int row = 0; row < nrows_; ++row, elt += ncols_) {
        if (matrix_[elt] < Constants::SMALL && !lineThroughRow_[row]) {
          lineThroughRow_[row] = true;
          lineThroughCol_[minCol] = true;
          assignRowToCol_[minCol] = row;
          assignColToRow_[row] = minCol;
#         ifdef DEBUG_HUNGARIAN
          mprintf("\tAssigned row %i to column %i\n", row, minCol);
#         endif
          assigned = 1;
          break;
        }
      }
    } else if (minRowZeros >= minColZeros) { // Preference given to rows
#     ifdef DEBUG_HUNGARIAN
      mprintf("Row %i has min # of zeros (%i)\n", minRow, minRowZeros);
#     endif
      // Assign row to the first unassigned col whose elt is zero.
      int elt = minRow * ncols_;
      for (int col = 0; col < ncols_; ++col, ++elt) {
        if (matrix_[elt] < Constants::SMALL && !lineThroughCol_[col]) {
          lineThroughRow_[minRow] = true;
          lineThroughCol_[col] = true;
          assignRowToCol_[col] = minRow;
          assignColToRow_[minRow] = col;
#         ifdef DEBUG_HUNGARIAN
          mprintf("\tAssigned row %i to column %i\n", minRow, col);
#         endif
          assigned = 1;
          break;
        }
      }
    }
    Nassigned += assigned;
  } // END while loop
# ifdef DEBUG_HUNGARIAN
  // DEBUG - Print assignments
  mprintf("%i ASSIGNMENTS:\n", Nassigned);
  for (int col = 0; col < ncols_; ++col)
    mprintf("\tAssigned row %i to column %i\n", assignRowToCol_[col], col);
# endif
  return Nassigned;
} 

/** \return Array containing which row idx matches which column idx. */
std::vector<int> Hungarian::Optimize() {
# ifdef DEBUG_HUNGARIAN
  PrintMatrix("INITIAL MATRIX");
# endif
  // Reduce elements in each row by the minimum in each row
  for (int row = 0; row < nrows_; ++row) {
    double minval = DBL_MAX;
    int elt0 = row * ncols_;
    int elt = elt0;
    for (int col = 0; col < ncols_; ++col, ++elt)
      if (matrix_[elt] < minval) minval = matrix_[elt];
    for (int col = 0; col < ncols_; ++col, ++elt0)
      matrix_[elt0] -= minval;
  }
# ifdef DEBUG_HUNGARIAN
  PrintMatrix("AFTER ROW REDUCTION");
# endif
  // Reduce elements in each col by the minimum in each col
  for (int col = 0; col < ncols_; ++col) {
    double minval = DBL_MAX;
    int elt0 = col;
    int elt = elt0;
    for (int row = 0; row < nrows_; ++row, elt += ncols_)
      if (matrix_[elt] < minval) minval = matrix_[elt];
    for (int row = 0; row < nrows_; ++row, elt0 += ncols_)
      matrix_[elt0] -= minval;
  }
# ifdef DEBUG_HUNGARIAN
  PrintMatrix("AFTER COL REDUCTION");
# endif
  // Loop until every row is assigned a column
  int max_iterations = nrows_ * ncols_;
  int iterations = 0;
  while (iterations < max_iterations) {
    // Attempt to assign each row to one column
    if (AssignRowsToColumns() == nrows_) break; // Assignments successful
    // Draw minimum number of lines required to cross out all zeros
    CoverZeroElements();
    // Update matrix according to lines
    UpdateMatrix();
    ++iterations;
  }
  return assignRowToCol_;
}  

/** Update matrix according to how lines have been drawn. The minimum
  * uncovered element is added to an element if covered twice, and
  * subtracted if uncovered. Elements which are covered once are left
  * alone.
  */
void Hungarian::UpdateMatrix() {
  // Find the minimum uncovered element
  double min_uncovered = DBL_MAX;
  for (int row = 0; row < nrows_; ++row) {
    if (!lineThroughRow_[row]) {
      for (int col = 0; col < ncols_; ++col) {
        if (!lineThroughCol_[col]) {
          double matrix_elt = matrix_.GetElement(row, col);
          if (matrix_elt < min_uncovered) {
#           ifdef DEBUG_HUNGARIAN
            mprintf("\tNew min %f at %i, %i\n",matrix_elt,row,col);
#           endif
            min_uncovered = matrix_elt;
          }
        }
      }
    }
  }
  // Update matrix elements
# ifdef DEBUG_HUNGARIAN
  mprintf("MIN UNCOVERED ELEMENT= %f\n", min_uncovered);
# endif
  int elt = 0;
  for (int row = 0; row < nrows_; ++row) {
    for (int col = 0; col < ncols_; ++col, ++elt) {
      if (lineThroughRow_[row] && lineThroughCol_[col]) // Marked twice: add
        matrix_[elt] += min_uncovered;
      else if (!lineThroughRow_[row] && !lineThroughCol_[col]) // Unmarked: subtract
        matrix_[elt] -= min_uncovered;
    }
  }
# ifdef DEBUG_HUNGARIAN
  PrintMatrix("AFTER ADDITION/SUBTRACTION OF MIN UNCOVERED");
# endif
}

#ifdef DEBUG_HUNGARIAN
// Hungarian::PrintLines()
void Hungarian::PrintLines(const char* title) {
  mprintf("%s\n", title);
  mprintf("\tLines through rows:");
  for (int row = 0; row < nrows_; ++row)
    if (lineThroughRow_[row]) mprintf(" %i", row);
  mprintf("\n\tLines through cols:");
  for (int col = 0; col < ncols_; ++col)
    if (lineThroughCol_[col]) mprintf(" %i", col);
  mprintf("\n");
}

// Hungarian::PrintMatrix()
void Hungarian::PrintMatrix(const char* Title) {
  mprintf("    %s\n",Title);
  int elt = 0;
  for (int row = 0; row < nrows_; ++row) {
    for (int col = 0; col < ncols_; ++col)
      mprintf(" %8.4f", matrix_[elt++]);
    mprintf("\n");
  }
}
#endif

// Hungarian::CoverZeroElements()
/** Cover all zeros in the matrix by marking as few rows and/or columns as 
  * possible. 
  */
int Hungarian::CoverZeroElements() {
  lineThroughRow_.assign(matrix_.Nrows(), false);
  lineThroughCol_.assign(matrix_.Ncols(), false);
  // Mark unassigned rows
  for (int row = 0; row < nrows_; ++row) {
    if (assignColToRow_[row] == -1)
      lineThroughRow_[row] = true;
  }
# ifdef DEBUG_HUNGARIAN
  PrintLines("MARK UNASSIGNED ROWS:");
# endif
  bool closed = false;
  while (!closed) {
    closed = true;
    // Mark columns with zeros in marked rows
    for (int row = 0; row < nrows_; ++row) {
      if (lineThroughRow_[row]) {
        int elt = row * ncols_;
        for (int col = 0; col < ncols_; ++col, ++elt) {
          if (matrix_.GetElement(row, col) < Constants::SMALL)
            lineThroughCol_[col] = true;
        }
      }
    }
#   ifdef DEBUG_HUNGARIAN
    PrintLines("MARK COLS WITH ZEROS IN MARKED ROWS:");
#   endif
    // Mark all rows having assignments in marked columns
    for (int col = 0; col < ncols_; ++col) {
      if (lineThroughCol_[col]) {
        for (int row = 0; row < nrows_; ++row) {
          if (!lineThroughRow_[row] && assignColToRow_[row] == col) {
            lineThroughRow_[row] = true;
            closed = false;
          }
        }
      }
    }
#   ifdef DEBUG_HUNGARIAN
    PrintLines("MARK ROWS HAVING ASSIGNMENTS IN MARKED COLS:");
#   endif
  }
  // Lines are drawn through marked columns and unmarked rows.
  for (std::vector<bool>::iterator mrow = lineThroughRow_.begin();
                                   mrow != lineThroughRow_.end(); ++mrow)
    *mrow = !(*mrow);
# ifdef DEBUG_HUNGARIAN
  PrintLines("FINAL LINES:");
# endif
  return 0;
} 
