#ifndef INC_DATASET_MODES_H
#define INC_DATASET_MODES_H
#include "DataSet_MatrixDbl.h"
/// Hold eigenvalues/eigenvectors and optionally averaged coords.
class DataSet_Modes : public DataSet {
  public:
    DataSet_Modes();
    ~DataSet_Modes();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Modes();}
    static const char* DeprecateFileMsg;
    // ----- DataSet functions -------------------
    size_t Size()                       const { return nmodes_; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    void Info()                         const { return;         }
    void Add( size_t, const void* ) {}
    void WriteBuffer(CpptrajFile&, SizeArray const&) const {} // TODO implement?
    int Allocate(SizeArray const&) { return 0; } // TODO implement?
    int Append(DataSet*) { return 1; }
    // -------------------------------------------
    typedef std::vector<double> Darray;
    typedef Darray::const_iterator AvgIt;
    AvgIt AvgBegin()                 const { return avgcrd_.begin(); } // TODO : Get rid of?
    Darray const& AvgCrd()           const { return avgcrd_; }
    Darray const& Mass()             const { return mass_;   }
    int NavgCrd()                    const { return (int)avgcrd_.size();  } // Project
    /// For reading directly into avgcrd buffer
    double* AvgFramePtr()                  { return &avgcrd_[0];          }
    const double* AvgFramePtr()      const { return &avgcrd_[0];          }
    void AllocateAvgCoords(int n)          { avgcrd_.resize(n, 0.0);      }

    int SetAvgCoords(DataSet_2D const&);
    int SetModes(bool, int, int, const double*, const double*);
    int CalcEigen(DataSet_2D const&,int);
    void PrintModes();
    int EigvalToFreq(double);
    int MassWtEigvect();
    int ReduceVectors();
    int Thermo(CpptrajFile&, int, double, double) const;

    double Eigenvalue(int i)         const { return evalues_[i];                } // IRED
    const double* Eigenvectors()     const { return evectors_;                  } // IRED
    const double* Eigenvector(int i) const { return evectors_ + (i * vecsize_); }
    int Nmodes()                     const { return nmodes_;                    } // Project
    int VectorSize()                 const { return vecsize_;                   } // Project
    bool IsReduced()                 const { return reduced_;                   }
  private:
    int ReduceCovar();
    int ReduceDistCovar();

    Darray avgcrd_;     ///< Average coordinates
    Darray mass_;       ///< Masses
    double* evalues_;   ///< Array of eigenvalues
    double* evectors_;  ///< Array of eigenvectors
    int nmodes_;        ///< Number of eigenmodes
    int vecsize_;       ///< Size of each eigenvector
    bool reduced_;      ///< True if modes have been reduced
    bool massWtEvecs_;  ///< True if eigenvectors have been mass-weighted
    bool evalsAreFreq_; ///< True if eigenvalues are in units of cm^-1
};
#endif
