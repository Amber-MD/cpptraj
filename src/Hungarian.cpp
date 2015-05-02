#include <cfloat> // DBL_MAX
#include "Hungarian.h"
#include "Constants.h" // SMALL
#ifdef DEBUG_HUNGARIAN
#  include "CpptrajStdio.h"
#endif

/** Initialize matrix for Hungarian algorithm. **/
int Hungarian::Initialize(size_t Ncols) {
  if (matrix_.resize( Ncols, Ncols )) return 1;
  lineThroughRow_.assign(matrix_.Nrows(), false);
  lineThroughCol_.assign(matrix_.Ncols(), false);
  assignRowToCol_.assign(matrix_.Ncols(), -1);
  assignColToRow_.assign(matrix_.Nrows(), -1);
  nrows_ = (int)matrix_.Nrows();
  ncols_ = (int)matrix_.Ncols();
  return 0;
}

// Hungarian::Optimize()
/** \return Array containing which row idx matches which column idx. */
std::vector<int> Hungarian::Optimize() {
# ifdef DEBUG_HUNGARIAN
  mprintf("----- HUNGARIAN ALGORITHM START ------------------------------------------------\n");
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
  // Loop until every row can be assigned a column
  int max_iterations = nrows_ * ncols_;
  int iterations = 0;
  while (iterations < max_iterations) {
    // Assign as many columns to rows as possible
    // If number of assignments == number of rows, done.
    if ( AssignRowsToColumns() == nrows_ ) break;
    // Cover zero elts with min # lines possible
    CoverZeroElements();
    // Update matrix according to min elements
    UpdateMatrix();
    ++iterations;
  }
# ifdef DEBUG_HUNGARIAN
  mprintf("--------------------------------------------------------------------------------\n");
# endif
  return assignRowToCol_;
}

// Hungarian::AssignRowsToColumns()
/** Try to assign columns to rows such that each row corresponds to one column
  * and the penalty incurred in each case is zero. This is done by iteratively
  * finding the row/col with the lowest number of zeros, then assigning that 
  * to the first unassigned col/row which has a zero. This is repeated until
  * no more assignments can be made.
  */
int Hungarian::AssignRowsToColumns() {
# ifdef DEBUG_HUNGARIAN
  mprintf("Assigning rows to columns:\n");
# endif
  int Nassigned = 0;
  assignRowToCol_.assign(matrix_.Ncols(), -1);
  assignColToRow_.assign(matrix_.Ncols(), -1);
  int assigned = 1;
  while (assigned > 0) { 
    assigned = 0;
    // Find row with lowest number of zeros
    int minRow = -1;
    int minRowZeros = nrows_ + 1;
    for (int row = 0; row < nrows_; ++row) {
      if (assignColToRow_[row] == -1) {
        int elt = row * ncols_;
        int nzeros = 0;
        for (int col = 0; col < ncols_; ++col, ++elt) {
          if (assignRowToCol_[col] == -1 && matrix_[elt] < Constants::SMALL) ++nzeros;
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
      if (assignRowToCol_[col] == -1) {
        int elt = col;
        int nzeros = 0;
        for (int row = 0; row < nrows_; ++row, elt += ncols_) {
          if (assignColToRow_[row] == -1 && matrix_[elt] < Constants::SMALL) ++nzeros;
        }
        if (nzeros != 0 && nzeros < minColZeros) {
          minColZeros = nzeros;
          minCol = col;
        }
      }
    }
    if (minRow == -1 && minCol == -1) break; // No zeros left
#   ifdef DEBUG_HUNGARIAN
    mprintf("  Min Row %i (%i zeros), min Column %i (%i zeros)\n", minRow, minRowZeros,
            minCol, minColZeros);
#   endif
    if (minColZeros < minRowZeros) {
#     ifdef DEBUG_HUNGARIAN 
      mprintf("\tColumn %i has min # of zeros (%i)\n", minCol, minColZeros);
#     endif
      // Assign column to the first unassigned row whose elt is zero.
      int elt = minCol; 
      for (int row = 0; row < nrows_; ++row, elt += ncols_) {
        if (matrix_[elt] < Constants::SMALL && assignColToRow_[row] == -1) {
          assignRowToCol_[minCol] = row;
          assignColToRow_[row] = minCol;
#         ifdef DEBUG_HUNGARIAN
          mprintf("\tAssigned row %i to column %i\n", row, minCol);
#         endif
          assigned = 1;
          break;
        }
      }
    } else if (minRowZeros <= minColZeros) { // Preference given to rows
#     ifdef DEBUG_HUNGARIAN
      mprintf("\tRow %i has min # of zeros (%i)\n", minRow, minRowZeros);
#     endif
      // Assign row to the first unassigned col whose elt is zero.
      int elt = minRow * ncols_;
      for (int col = 0; col < ncols_; ++col, ++elt) {
        if (matrix_[elt] < Constants::SMALL && assignRowToCol_[col] == -1) {
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
  mprintf("  %i Assignments:\n", Nassigned);
  for (int col = 0; col < ncols_; ++col)
    mprintf("\tAssigned row %i to column %i\n", assignRowToCol_[col], col);
# endif
  return Nassigned;
}

// Hungarian::UpdateMatrix()
/** Update matrix according to how lines have been drawn. The minimum
  * uncovered element is added to an element if covered twice, and
  * subtracted if uncovered. Elements which are covered once are left
  * alone.
  */
void Hungarian::UpdateMatrix() {
# ifdef DEBUG_HUNGARIAN
  mprintf("Updating matrix.\n");
# endif
  // Find the minimum uncovered element
  double min_uncovered = DBL_MAX;
  for (int row = 0; row < nrows_; ++row) {
    if (!lineThroughRow_[row]) {
      for (int col = 0; col < ncols_; ++col) {
        if (!lineThroughCol_[col]) {
          double matrix_elt = matrix_.element(col, row);
          if (matrix_elt < min_uncovered) {
#           ifdef DEBUG_HUNGARIAN
            mprintf("\tNew min %f at col=%i, row=%i\n",matrix_elt,col,row);
#           endif
            min_uncovered = matrix_elt;
          }
        }
      }
    }
  }
  // Update matrix elements
# ifdef DEBUG_HUNGARIAN
  mprintf("  Min Uncovered Element= %f\n", min_uncovered);
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
  PrintMatrix("After Addition/Subtraction of Min Uncovered Element");
# endif
}

#ifdef DEBUG_HUNGARIAN
// Hungarian::PrintLines()
void Hungarian::PrintLines(const char* title) {
  mprintf("  %s\n", title);
  mprintf("\tLines through rows:");
  for (int row = 0; row < nrows_; ++row)
    if (lineThroughRow_[row]) mprintf(" %i", row);
  mprintf("\n\tLines through cols:");
  for (int col = 0; col < ncols_; ++col)
    if (lineThroughCol_[col]) mprintf(" %i", col);
  mprintf("\n");
  int elt = 0;
  for (int row = 0; row < nrows_; ++row) {
    for (int col = 0; col < ncols_; ++col) {
      if (lineThroughRow_[row] && lineThroughCol_[col])
        mprintf("-%6s", "--|---");
      else if (lineThroughRow_[row])
        mprintf("-%6s", "------");
      else if (lineThroughCol_[col])
        mprintf(" %6s", "  |   ");
      else
        mprintf(" %6.2f", matrix_[elt]);
      elt++;
    }
    mprintf("\n");
  }
}

// Hungarian::PrintMatrix()
void Hungarian::PrintMatrix(const char* Title) {
  mprintf("  %s\n",Title);
  int elt = 0;
  for (int row = 0; row < nrows_; ++row) {
    for (int col = 0; col < ncols_; ++col)
      mprintf(" %6.2f", matrix_[elt++]);
    mprintf("\n");
  }
}
#endif

// Hungarian::CoverZeroElements()
/** Cover all zeros in the matrix by drawing lines through as few rows 
  * and/or columns as possible.
  */
void Hungarian::CoverZeroElements() {
  /* New algorithm:
   *   0: Calc # zeros in each row/col.
   *   1: On first iteration, pick row/col with max # zeros and mark, tie goes to row.
   *      On subsequent, pick row/col with max (unmarked zeros - unmarked nonzeros).
   *   2: Update zero/nonzero counts, also count unmarked non-zero cells.
   *   3: If no more zeros exit, otherwise go to 1.
   */
  lineThroughRow_.assign(matrix_.Nrows(), false);
  lineThroughCol_.assign(matrix_.Ncols(), false);
  typedef std::vector<int> Iarray;
  Iarray rowZeroCount(nrows_, 0);
  Iarray colZeroCount(ncols_, 0);
  Iarray rowNonzeroCount(nrows_, 0);
  Iarray colNonzeroCount(ncols_, 0);
  // Step 0: Get initial zero counts.
  int TotalZeros = 0;
  Matrix<double>::iterator elt = matrix_.begin(); // TODO: const_iterator
  for (int row = 0; row != nrows_; ++row)
  {
    for (int col = 0; col != ncols_; ++col)
    {
      if (*elt < Constants::SMALL) {
        rowZeroCount[row]++;
        colZeroCount[col]++;
        ++TotalZeros;
      } else {
        rowNonzeroCount[row]++;
        colNonzeroCount[col]++;
      }
      ++elt;
    }
  }
# ifdef DEBUG_HUNGARIAN
  mprintf("Initial zero counts (%i total):\n", TotalZeros);
# endif
  // Loop while there are still unmarked zeros.
  unsigned int iter = 0;
  while (TotalZeros > 0) {
    int maxDiff = -1;
    int maxZeroCount = -1;
    int maxIdx = -1;
    bool maxIsRow = true;
    // Find max # of zeros in rows.
    for (int row = 0; row != nrows_; ++row) {
      if (!lineThroughRow_[row]) {
#       ifdef DEBUG_HUNGARIAN
        if (iter == 0)
        mprintf("\tRow %i has %i zeros.\n",row,rowZeroCount[row]);
        //mprintf("\tRow %i has %i zeros, %i nonzeros.\n",row,rowZeroCount[row],rowNonzeroCount[row]);
#       endif
        int diff;
        if (iter > 0)
          diff = rowZeroCount[row] - rowNonzeroCount[row];
        else
          diff = rowZeroCount[row];
        if (maxZeroCount == -1 || diff > maxDiff)
        {
          maxDiff = diff;
          maxZeroCount = rowZeroCount[row];
          maxIdx = row;
        }
      }
    }
    // Find max # of zeros in cols.
    for (int col =0; col != ncols_; ++col) {
      if (!lineThroughCol_[col]) {
#       ifdef DEBUG_HUNGARIAN
        if (iter == 0)
        mprintf("\tCol %i has %i zeros.\n",col,colZeroCount[col]);
        //mprintf("\tCol %i has %i zeros, %i nonzeros.\n",col,colZeroCount[col],colNonzeroCount[col]);
#       endif
        int diff;
        if (iter > 0)
          diff = colZeroCount[col] - colNonzeroCount[col];
        else
          diff = colZeroCount[col];
        if (diff > maxDiff)
        {
          maxDiff = diff;
          maxZeroCount = colZeroCount[col];
          maxIdx = col;
          maxIsRow = false;
        }
      }
    }
    // Cross out row or column with max.
#   ifdef DEBUG_HUNGARIAN
    if (maxIsRow)
      mprintf("Max zero count %i in row %i\n", maxDiff, maxIdx);
      //mprintf("Max zero count %i (%i) in row %i\n", maxDiff, maxZeroCount, maxIdx);
    else
      mprintf("Max zero count %i in col %i\n", maxDiff, maxIdx);
      //mprintf("Max zero count %i (%i) in col %i\n", maxDiff, maxZeroCount, maxIdx);
#   endif
    if (maxIsRow) {
      // Cross out row, update each column
      lineThroughRow_[maxIdx] = true;
      elt = matrix_.begin() + (maxIdx * ncols_);
      //mprintf("DEBUG: Going across row %i:", maxIdx);
      for (int col = 0; col != ncols_; ++col, ++elt) {
        //mprintf(" %g", *elt);
        if (*elt < Constants::SMALL)
          colZeroCount[col]--;
        else
          colNonzeroCount[col]--;
      }
    } else {
      // Cross out column, update each row.
      lineThroughCol_[maxIdx] = true;
      elt = matrix_.begin() + maxIdx;
      //mprintf("DEBUG: Going down col %i:", maxIdx);
      for (int row = 0; row != nrows_; ++row, elt += ncols_) {
        //mprintf(" %g", *elt);
        if (*elt < Constants::SMALL)
          rowZeroCount[row]--;
        else
          rowNonzeroCount[row]--;
      }
    }
    // Update total # zeros
    TotalZeros -= maxZeroCount;
#   ifdef DEBUG_HUNGARIAN
    mprintf("Current zero counts (%i total):\n", TotalZeros);
#   endif
    ++iter;
  }
/* // OLD ALGORITHM
# ifdef DEBUG_HUNGARIAN
  mprintf("Drawing lines through rows/cols with zero elements\n");
# endif
  lineThroughRow_.assign(matrix_.Nrows(), false);
  lineThroughCol_.assign(matrix_.Ncols(), false);
  typedef std::vector<bool> Barray; // FIXME: In header?
  Barray markedRow(matrix_.Nrows(), false);
  Barray markedCol(matrix_.Ncols(), false);
  // Put minimal number of lines through rows/cols so all zeros covered.
  int Nmarks = 1;
  while (Nmarks > 0) {
    Nmarks = 0;
    // Mark all rows having no assignments.
    for (int row = 0; row < nrows_; row++) {
      if (!markedRow[row]) {
        if (assignColToRow_[row] == -1) {
          markedRow[row] = true;
#         ifdef DEBUG_HUNGARIAN
          mprintf("\tMarking row %i (has no assignment)\n", row);
#         endif
          Nmarks++;
        }
      }
    }
    // Mark all columns having zeros in marked rows
    for (int col = 0; col < ncols_; col++) {
      if (!markedCol[col]) {
        int elt = col;
        for (int row = 0; row < nrows_; ++row, elt += ncols_) {
          if (markedRow[row] && matrix_[elt] < Constants::SMALL) {
            markedCol[col] = true;
#           ifdef DEBUG_HUNGARIAN
            mprintf("\tMarking col %i (has zero in marked row %i)\n", col, row);
#           endif
            Nmarks++;
            break;
          }
        }
      }
    }
    // Mark all rows with zero in marked columns
    for (int row = 0; row < nrows_; row++) {
      if (!markedRow[row]) {
        int elt = row * ncols_;
        for (int col = 0; col < ncols_; ++col, ++elt) {
          if (markedCol[col] && matrix_[elt] < Constants::SMALL) {
            markedRow[row] = true;
#           ifdef DEBUG_HUNGARIAN
            mprintf("\tMarking row %i (has zero in marked column %i)\n", row, col);
#           endif
            Nmarks++;
            break;
          }
        }
      }
    }
  }
  // Draw lines through all marked columns and unmarked rows
  for (int col = 0; col < ncols_; col++)
    if (markedCol[col]) lineThroughCol_[col] = true;
  for (int row = 0; row < nrows_; row++)
    if (!markedRow[row]) lineThroughRow_[row] = true;
*/
# ifdef DEBUG_HUNGARIAN
  //mprintf("  Assigned %i lines\n", Nlines);
  PrintLines("Matrix With Lines:");
# endif
}
