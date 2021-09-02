#ifndef INC_CLUSTER_CMATRIX_BINARY_H
#define INC_CLUSTER_CMATRIX_BINARY_H
#include <cstddef> // size_t
#include "../CpptrajFile.h" 
class FileName;
namespace Cpptraj {
namespace Cluster {
class Cframes;
/// Used to read pairwise distance cache files in cpptraj binary format.
class Cmatrix_Binary {
  public:
    Cmatrix_Binary();

    /// \return true if file is binary cpptraj cluster matrix file.
    static bool ID_Cmatrix(CpptrajFile&);

    /// \return Sieve value.
    int Sieve()          const { return sieve_; }
    /// \return Actual number of rows in matrix.
    size_t ActualNrows() const { return actual_nrows_; }
    /// \return Number of total frames originally.
    size_t Ntotal()      const { return ntotal_;} 

    /// Open cluster matrix file for reading. Set sieve and actual # rows. 
    int OpenCmatrixRead(FileName const&);
    /// Get cluster matrix element (col, row)
    double GetCmatrixElement(unsigned int, unsigned int);
    /// Get cluster matrix element (raw index)
    double GetCmatrixElement(unsigned int);
    /// Read cmatrix into given pointer TODO should this take Matrix/StatusArray?
    int GetCmatrix(float*, char*);

    /// Write cluster matrix TODO add a setup routine
    static int WriteCmatrix(FileName const&, const float*, Cframes const&, size_t, int);
  private:
    static const unsigned char Magic_[];
    /// For reading/writing 8 byte unsigned integers
    typedef unsigned long long int uint_8;
    /// For reading/writing 8 byte signed integers
    typedef long long int sint_8;
    
    CpptrajFile file_;
    int sieve_;
    size_t actual_nrows_;
    size_t ntotal_;       ///< Total number of frames in original data.
    off_t headerOffset_;
};

}
}
#endif
