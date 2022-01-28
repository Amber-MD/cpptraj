#ifndef INC_CLUSTER_DYNAMICMATRIX_H
#define INC_CLUSTER_DYNAMICMATRIX_H
#include <vector>
#include "../Matrix.h"
#include "../CpptrajStdio.h" // DEBUG
namespace Cpptraj {
namespace Cluster {

/// Hold between-cluster distances and whether cluster is still present.
class DynamicMatrix {
  public:
    DynamicMatrix() {}
    /// Indicate given row/col should be ignored.
    //void Ignore(int row)            { ignore_[row] = true;   }
    inline void Ignore(int);
    /// \return true if given row/col has been ignored.
    bool IgnoringRow(int row) const { return ignore_[row];   }
    /// \return Original number of rows in matrix
    size_t Nrows()            const { return Mat_.Nrows();   }
    /// Set the row and column of the smallest element not being ignored.
#   ifdef _OPENMP
    double FindMin(int&, int&);
#   else
    double FindMin(int&, int&) const;
#   endif
    /// \return an element.
    inline double GetCdist(int c, int r) const { return Mat_.element(c,r); }
    /// Print all matrix elements to STDOUT
    void PrintElements() const;
    /// Print matrix elements to STDOUT as a square
    void PrintElementsSquare() const;
    /// Print closest index to each col to STDOUT
    void PrintClosest() const;
    /// Add given element to matrix.
    //int AddCdist(double d)        { return Mat_.addElement((float)d); }
    /// Set element at column/row to given value
    //void SetCdist(int col, int row, double val) { Mat_.setElement(col, row, val); }
    inline void SetCdist(int, int, float);
    /// Set up matrix for given number of rows
    int SetupMatrix(size_t);
  private:
    inline void updateClosestIdx(unsigned int);

    Matrix<float> Mat_;        ///< Upper-triangle matrix holding cluster distances.
    std::vector<bool> ignore_; ///< If true, ignore the row/col when printing/searching etc.
    std::vector<int> closestIdx_; ///< For each cluster, index of the cluster closest to it
#   ifdef _OPENMP
    std::vector<int> minRow_;
    std::vector<int> minCol_;
    std::vector<float> minVal_;
#   endif
};

/** Update closest to idx */
void DynamicMatrix::updateClosestIdx(unsigned int idx) {
      closestIdx_[idx] = -1;
      float current_closest = 0;
      for (unsigned int jdx = 0; jdx != closestIdx_.size(); jdx++) {
        if (!ignore_[jdx] && idx != jdx) {
          if (closestIdx_[idx] == -1) {
            closestIdx_[idx] = jdx;
            current_closest = Mat_.element(idx, jdx);
          } else {
            float fdist = Mat_.element(idx, jdx);
            if (fdist < current_closest) {
              closestIdx_[idx] = jdx;
              current_closest = fdist;
            }
          }
        }
      } // END loop jdx
      mprintf("DEBUG: New closest to %u is %i\n", idx, closestIdx_[idx]);
}

/** Set cluster distance. Record closest distance for each col. */
void DynamicMatrix::SetCdist(int col, int row, float val) {
  mprintf("DEBUG: Setting cluster %i to cluster %i distance= %f\n", col, row, val);

  int update_col = -1;
  int update_row = -1;

  // Check col->row distance
  if (closestIdx_[col] < 0) {
    mprintf("DEBUG: This is the initial distance for col.\n");
    closestIdx_[col] = row;
  } else {
    float closest_dist_to_col = Mat_.element(col, closestIdx_[col]);
    mprintf("DEBUG: Current closest to %i is %i (%f)\n", col, closestIdx_[col], closest_dist_to_col); 
    if (val < closest_dist_to_col) {
      mprintf("DEBUG: Updated col.\n");
      closestIdx_[col] = row;
    } else if (row == closestIdx_[col] && val > closest_dist_to_col) {
      // We are increasing what was previously the closest distance. May need to find new closest distance to col.
      update_col = col;
    }
  }

  // Check row->col distance
  if (closestIdx_[row] < 0) {
    mprintf("DEBUG: This is the initial distance for row.\n");
    closestIdx_[row] = col;
  } else {
    float closest_dist_to_row = Mat_.element(row, closestIdx_[row]);
    mprintf("DEBUG: Current closest to %i is %i (%f)\n", row, closestIdx_[row], closest_dist_to_row);
    if (val < closest_dist_to_row) {
      mprintf("DEBUG: Updated row.\n");
      closestIdx_[row] = col;
    } else if (col == closestIdx_[row] && val > closest_dist_to_row) {
      // We are increasing what was previously the closest distance. May need to find new closest distance to row.
      update_row = row;
    }
  }

  // Actually update the element
  Mat_.setElement(col, row, val);

  // Update closest if necessary
  if (update_col != -1) {
    mprintf("DEBUG: Closest to column %i increased. Need to update.\n", update_col);
    updateClosestIdx(update_col);
  }
  if (update_row != -1) {
    mprintf("DEBUG: Closest to row %i increased. Need to update.\n", update_row);
    updateClosestIdx(update_row);
  }
}

/** Specify given row should be ignored. Update closest distances if necessary. */
void DynamicMatrix::Ignore(int row) {
  ignore_[row] = true;
  // If any remaining indices have the ignored row that needs to be updated
  for (unsigned int idx = 0; idx != closestIdx_.size(); idx++) {
    if (!ignore_[idx] && closestIdx_[idx] == row) {
      mprintf("DEBUG: closest to %u was just ignored; must be updated.\n", idx);
      updateClosestIdx(idx);
/*      closestIdx_[idx] = -1;
      float current_closest = 0;
      for (unsigned int jdx = 0; jdx != closestIdx_.size(); jdx++) {
        if (!ignore_[jdx] && idx != jdx) {
          if (closestIdx_[idx] == -1) {
            closestIdx_[idx] = jdx;
            current_closest = Mat_.element(idx, jdx);
          } else {
            float fdist = Mat_.element(idx, jdx);
            if (fdist < current_closest) {
              closestIdx_[idx] = jdx;
              current_closest = fdist;
            }
          }
        }
      } // END loop jdx
      mprintf("DEBUG: New closest to %u is %i\n", idx, closestIdx_[idx]);*/
    }
  } // END loop idx
}

}
}
#endif
