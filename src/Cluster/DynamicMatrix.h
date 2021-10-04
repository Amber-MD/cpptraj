#ifndef INC_CLUSTER_DYNAMICMATRIX_H
#define INC_CLUSTER_DYNAMICMATRIX_H
#include <vector>
#include "../Matrix.h"
namespace Cpptraj {
namespace Cluster {

/// Hold between-cluster distances and whether cluster is still present.
class DynamicMatrix {
  public:
    DynamicMatrix() {}
    /// Indicate given row/col should be ignored.
    void Ignore(int row)            { ignore_[row] = true;   }
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
    /// Add given element to matrix.
    int AddCdist(double d)        { return Mat_.addElement((float)d); }
    /// Set element at column/row to given value
    void SetCdist(int col, int row, double val) { Mat_.setElement(col, row, val); }
    /// Set up matrix for given number of rows
    int SetupMatrix(size_t);
  private:
    Matrix<float> Mat_;        ///< Upper-triangle matrix holding cluster distances.
    std::vector<bool> ignore_; ///< If true, ignore the row/col when printing/searching etc.
#   ifdef _OPENMP
    std::vector<int> minRow_;
    std::vector<int> minCol_;
    std::vector<float> minVal_;
#   endif
};

}
}
#endif
