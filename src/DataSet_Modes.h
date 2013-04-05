#ifndef INC_DATASET_MODES_H
#define INC_DATASET_MODES_H
#include "DataSet_MatrixDbl.h"
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
    void SetAvgCoords(int, const double*);
    int CalcEigen(DataSet_2D const&,int);
    void PrintModes();
    int WriteToFile(std::string const&);
    int ReadEvecFile(std::string const&, int, int);
    int EigvalToFreq();
    int MassWtEigvect( const double* );
    int ReduceCovar();
    int ReduceDistCovar(int);

    void SetType( DataSet_MatrixDbl::MatrixType typeIn ) { type_ = typeIn; }

    const double* AvgCrd()            { return avgcrd_;                    }
    const double* Eigenvalues()       { return evalues_;                   } 
    double Eigenvalue(int i)          { return evalues_[i];                }
    const double* Eigenvectors()      { return evectors_;                  }
    const double* Eigenvector(int i)  { return evectors_ + (i * vecsize_); }
    int Nmodes()                      { return nmodes_;                    }
    int VectorSize()                  { return vecsize_;                   }
    int NavgCrd()                     { return navgcrd_;                   }
    DataSet_MatrixDbl::MatrixType Type() { return type_;                      }
  private:
    double* avgcrd_;
    double* evalues_;
    double* evectors_;
    int nmodes_;
    int vecsize_;
    int navgcrd_;
    DataSet_MatrixDbl::MatrixType type_;
    bool reduced_;
};
#endif
