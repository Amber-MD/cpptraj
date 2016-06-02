#ifndef INC_CLUSTERMATRIX_H
#define INC_CLUSTERMATRIX_H
#include <vector>
#include "Matrix.h"
/// Used to hold distances between clusters.
class ClusterMatrix {
  public:
    ClusterMatrix() {}
    /// Indicate given row/col should be ignored.
    void Ignore(int row)            { ignore_[row] = true;   }
    /// \return true if given row/col has been ignored.
    bool IgnoringRow(int row) const { return ignore_[row];   }
    /// \return Original number of rows
    size_t Nframes()          const { return ignore_.size(); }
    /// Set the row and column of the smallest element not being ignored.
    double FindMin(int&, int&) const;
  private:
    Matrix<float> Mat_;        ///< Upper-triangle matrix holding cluster distances.
    std::vector<bool> ignore_; ///< If true, ignore the row/col when printing/searching etc.
};
#endif
