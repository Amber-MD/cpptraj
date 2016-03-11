#ifndef INC_DATASET_CMATRIX_H
#define INC_DATASET_CMATRIX_H
#include "DataSet_2D.h"
#include "Matrix.h"
#include "ClusterSieve.h"
/// Used to hold pairwise distance matrix for clustering.
class DataSet_Cmatrix : public DataSet_2D {
  public:
    DataSet_Cmatrix() : DataSet_2D(CMATRIX, TextFormat(TextFormat::DOUBLE, 12, 4)) {}
    DataSet_Cmatrix(const DataSet_Cmatrix&);
    DataSet_Cmatrix& operator=(const DataSet_Cmatrix&);
    /// Access internal matrix pointer to interface with file IO
    float*       Ptr()                         { return Mat_.Ptr();         }
    float const* Ptr()                   const { return Mat_.Ptr();         }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return Mat_.size();        }
    void Info()                          const { return;                    }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    // ----- DataSet_2D functions ----------------
    int Allocate2D(size_t x,size_t y)          { return 1;                         }
    int AllocateHalf(size_t x)                 { return 1;                         }
    int AllocateTriangle(size_t x)             { return Mat_.resize(0,x);          }
    double GetElement(size_t x,size_t y) const { return (double)Mat_.element(x,y); }
    double GetElement(size_t i)          const { return (double)Mat_[i];           }
    /// \return Actual number of rows in matrix.
    size_t Nrows()                       const { return Mat_.Nrows();              }
    size_t Ncols()                       const { return Mat_.Ncols();              }
    double* MatrixArray()                const;
    MatrixKindType MatrixKind()          const { return TRI;                       }
    // -------------------------------------------
    /// Set up matrix with sieve value.
    int SetupWithSieve(size_t, size_t, int);
    /// Allocate ignore array for given # of original rows and sieve.
    int SetupIgnore(size_t, std::vector<char> const&, int);

    //int SaveFile(std::string const&) const;
    //int LoadFile(std::string const&, int);
    int SetupMatrix(size_t);
    /// Indicate given row/col should be ignored.
    void Ignore(int row)            { ignore_[row] = true;   }
    /// \return true if given row/col has been ignored.
    bool IgnoringRow(int row) const { return ignore_[row];   }
    /// \return Number of frames (original nrows)
    size_t Nframes()          const { return ignore_.size(); }
    /// \return An array containing sieved frame numbers.
    ClusterSieve::SievedFrames Sieved() const { return sievedFrames_.Frames(); }
    /// \return Sieve value
    int SieveValue()                    const { return sievedFrames_.Sieve();  }
    /// \return Sieve type
    ClusterSieve::SieveType SieveType() const { return sievedFrames_.Type();   }
    /// Set the row and column of the smallest element.
    double FindMin(int&, int&) const;
    void PrintElements() const;
    /// \return an element indexed by sievedFrames.
    inline double GetFdist(int, int) const;
    /// \return an element.
    inline double GetCdist(int c, int r) const { return Mat_.element(c,r); }
    inline void SetElement(int, int, double);
    size_t Nelements()        const { return Mat_.size();               }
    int AddElement(double d)        { return Mat_.addElement((float)d); }
    size_t DataSize() const;
    typedef Matrix<float>::iterator const_iterator;
    const_iterator begin() const { return Mat_.begin(); }
    const_iterator end()   const { return Mat_.end();   }
  private:
    static const unsigned char Magic_[];
    /// For reading/writing 8 byte unsigned integers
    typedef unsigned long long int uint_8;
    /// For reading/writing 8 byte signed integers
    typedef long long int sint_8;
    /// If true, ignore the row/col when printing/searching etc
    std::vector<bool> ignore_;
    Matrix<float> Mat_;
    ClusterSieve sievedFrames_;
};
// Inline functions
double DataSet_Cmatrix::GetFdist(int col, int row) const {
  // row and col are based on original; convert to reduced
  // FIXME: This assumes GetElement will never be called for a sieved frame.
  return Mat_.element(sievedFrames_.FrameToIdx(col), 
                      sievedFrames_.FrameToIdx(row));
}

void DataSet_Cmatrix::SetElement(int col, int row, double val) {
  Mat_.setElement(col, row, val);
}
#endif
