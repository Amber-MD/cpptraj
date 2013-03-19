#ifndef INC_CLUSTERMATRIX_H
#define INC_CLUSTERMATRIX_H
#include "TriangleMatrix.h"
class ClusterMatrix : private TriangleMatrix {
  public:
    ClusterMatrix() : sieve_(1) {}
    /// Set up matrix with sieve value.
    ClusterMatrix(size_t, size_t);
    ClusterMatrix(const ClusterMatrix&);
    ClusterMatrix& operator=(const ClusterMatrix&);

    int SaveFile(std::string const&) const;
    int LoadFile(std::string const&, int);
    int SetupMatrix(size_t);
    /// Indicate given row/col should be ignored.
    void Ignore(int row)            { ignore_[row] = true;   }
    /// \return true if given row/col has been ignored.
    bool IgnoringRow(int row) const { return ignore_[row];   }
    /// \return Number of frames (original nrows)
    size_t Nframes()          const { return ignore_.size(); }
    /// Set the row and column of the smallest element.
    double FindMin(int&, int&) const;
    void PrintElements() const;
    inline double GetElement(int, int) const;
    inline void SetElement(int, int, double);
    using TriangleMatrix::Nelements;
    using TriangleMatrix::AddElement;
  private:
    static const unsigned char Magic_[];
    /// For reading/writing 8 byte integers
    typedef unsigned long long int uint_8;
    /// If true, ignore the row/col when printing/searching etc
    std::vector<bool> ignore_;
    size_t sieve_;
};
// Inline functions
double ClusterMatrix::GetElement(int row, int col) const {
  // row and col are based on original; convert to reduced
  // FIXME: This assumes GetElement will never be called for a sieved frame.
  return TriangleMatrix::GetElement(row / sieve_, col / sieve_);
}

void ClusterMatrix::SetElement(int row, int col, double val) {
  return TriangleMatrix::SetElement(row / sieve_, col / sieve_, val);
}
#endif
