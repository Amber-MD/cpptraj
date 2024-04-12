#ifndef INC_DATASET_MODES_H
#define INC_DATASET_MODES_H
#include "DataSet.h"
class DataSet_2D;
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
    size_t MemUsageInBytes() const;
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
    /// For reading directly into mass
    double* MassPtr()                      { return &mass_[0]; }
    /// For reading directly into eigenvlues array
    double* EvalPtr()                      { return evalues_; }
    /// For reading directly into eigenvectors array
    double* EvectPtr()                     { return evectors_; }

    int SetAvgCoords(DataSet_2D const&);
    int SetAvgCoords(std::vector<double> const&); // TODO deprecate above version for this one
    int SetModes(bool, int, int, const double*, const double*);
    /// Allocate memory for modes data
    int AllocateModes(unsigned int, unsigned int, unsigned int, unsigned int);
    /// Calculate all eigenvalues/eigenvectors for given matrix
    int CalcEigen_General(DataSet_2D const&);
    /// Calculate specific # of eigenvalues/eigenvectors for given symmetric matrix
    int CalcEigen(DataSet_2D const&,int);
    /// Multiply all elements of the specified eigenvector by given value
    void MultiplyEvecByFac(int, double);
    /// Descending sort by absolute value of eigenvalues
    void SortByAbsEigenvalue();
    void PrintModes();
    int EigvalToFreq(double);
    int MassWtEigvect();
    int ReduceVectors();
    /// Resize so that only the first N modes are saved
    int ResizeModes(int);
    int Thermo(CpptrajFile&, int, double, double) const;

    double Eigenvalue(int i)         const { return evalues_[i];                } // IRED
    const double* EigenvaluePtr()    const { return evalues_; }
    const double* Eigenvectors()     const { return evectors_;                  } // IRED
    const double* Eigenvector(int i) const { return evectors_ + (i * vecsize_); }
    int Nmodes()                     const { return nmodes_;                    } // Project
    int VectorSize()                 const { return vecsize_;                   } // Project
    bool IsReduced()                 const { return reduced_;                   }
    bool EvecsAreMassWtd()           const { return evecsAreMassWtd_;           }
    bool EvalsAreFreq()              const { return evalsAreFreq_;              }
  private:
    int ReduceCovar();
    int ReduceDistCovar();

    Darray avgcrd_;        ///< Average coordinates
    Darray mass_;          ///< Masses
    double* evalues_;      ///< Array of eigenvalues
    double* evectors_;     ///< Array of eigenvectors
    int nmodes_;           ///< Number of eigenmodes
    int vecsize_;          ///< Size of each eigenvector
    bool reduced_;         ///< True if modes have been reduced
    bool evecsAreMassWtd_; ///< True if eigenvectors have been mass-weighted
    bool evalsAreFreq_;    ///< True if eigenvalues are in units of cm^-1
};
#endif
