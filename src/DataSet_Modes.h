#ifndef INC_DATASET_MODES_H
#define INC_DATASET_MODES_H
#include "DataSet_MatrixDbl.h"
#include "Frame.h"
#include "Analysis.h"
/// Hold eigenvalues/eigenvectors and optionally averaged coords.
class DataSet_Modes : public DataSet {
  public:
    DataSet_Modes();
    ~DataSet_Modes();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Modes();}
    static const char* DeprecateFileMsg;
    // ----- DataSet functions -------------------
    size_t Size() const { return nmodes_; }
    int Sync()          { return 1;       }
    void Info()   const { return;         }
    void Add( size_t, const void* ) {} // TODO: Get rid of, only needed by DataSet_1D
    // -------------------------------------------
    void SetAvgCoords(DataSet_2D const&);
    void AllocateAvgCoords(int n)          { avgcrd_.SetupFrame( n / 3 ); }
    const Frame& AvgFrame()          const { return avgcrd_;              } // Project, Modes
    int NavgCrd()                    const { return (int)avgcrd_.size();  } // Project
    /// For reading directly into avgcrd buffer
    double* AvgFramePtr()                  { return avgcrd_.xAddress();   }

    int SetModes(bool, int, int, const double*, const double*);
    int CalcEigen(DataSet_2D const&,int);
    void PrintModes();
    int EigvalToFreq(double);
    int MassWtEigvect( DataSet_MatrixDbl::Darray const& );
    int ReduceVectors();
    int Thermo(CpptrajFile&, int, double, double) const;
    int NMWiz(CpptrajFile&, int, std::string, Topology const& );

    void SetType( DataSet_2D::MatrixType typeIn ) { type_ = typeIn; }

    double Eigenvalue(int i)         const { return evalues_[i];                } // IRED
    const double* Eigenvectors()     const { return evectors_;                  } // IRED
    const double* Eigenvector(int i) const { return evectors_ + (i * vecsize_); }
    int Nmodes()                     const { return nmodes_;                    } // Project
    int VectorSize()                 const { return vecsize_;                   } // Project
    DataSet_2D::MatrixType Type()    const { return type_;                      } // Project
    bool IsReduced()                 const { return reduced_;                   }
  private:
    int ReduceCovar();
    int ReduceDistCovar();

    Frame avgcrd_;                ///< Average coordinates
    double* evalues_;             ///< Array of eigenvalues
    double* evectors_;            ///< Array of eigenvectors
    int nmodes_;                  ///< Number of eigenmodes
    int vecsize_;                 ///< Size of each eigenvector
    DataSet_2D::MatrixType type_; ///< Type of matrix modes were generated from
    bool reduced_;                ///< True if modes have been reduced
};
#endif
