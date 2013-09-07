#ifndef INC_DATASET_MODES_H
#define INC_DATASET_MODES_H
#include "DataSet_MatrixDbl.h"
#include "Frame.h"
/// Hold eigenvalues/eigenvectors and optionally averaged coords.
class DataSet_Modes : public DataSet {
  public:
    DataSet_Modes();
    ~DataSet_Modes();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Modes();}
    // ----- DataSet functions -------------------
    size_t Size() const { return nmodes_; }
    int Sync()          { return 1;       }
    void Info()   const { return;         }
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
    void SetAvgCoords(DataSet_2D const&);
    int CalcEigen(DataSet_2D const&,int);
    void PrintModes();
    int WriteToFile(std::string const&);
    int ReadEvecFile(std::string const&, int, int);
    int EigvalToFreq();
    int MassWtEigvect( DataSet_MatrixDbl::Darray const& );
    int Reduce();
    int Thermo(CpptrajFile&, int, double, double) const;

    void SetType( DataSet_2D::MatrixType typeIn ) { type_ = typeIn; }

    const double* AvgCrd()           const { return (const double*)avgcrd_.xAddress(); } // Project
    const double* Eigenvalues()      const { return evalues_;                   } // IRED 
    double Eigenvalue(int i)         const { return evalues_[i];                } // IRED
    const double* Eigenvectors()     const { return evectors_;                  } // IRED
    const double* Eigenvector(int i) const { return evectors_ + (i * vecsize_); }
    int Nmodes()                     const { return nmodes_;                    } // Project
    int VectorSize()                 const { return vecsize_;                   } // Project
    int NavgCrd()                    const { return (int)avgcrd_.size();        } // Project
    DataSet_2D::MatrixType Type()    const { return type_;                      } // Project
  private:
    int ReduceCovar();
    int ReduceDistCovar();

    Frame avgcrd_;
    double* evalues_;
    double* evectors_;
    int nmodes_;
    int vecsize_;
    DataSet_2D::MatrixType type_;
    bool reduced_;
};
#endif
