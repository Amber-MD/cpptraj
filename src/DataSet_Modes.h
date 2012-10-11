#ifndef INC_DATASET_MODES_H
#define INC_DATASET_MODES_H
#include "DataSet.h"
#include "Frame.h"
#include "DataSet_Matrix.h"
/// Hold eigenvalues/eigenvectors and optionally averaged coords.
class DataSet_Modes : public DataSet {
  public:
    DataSet_Modes();
    ~DataSet_Modes();

    // ---------- DataSet routines
    int Size() { return nmodes_; }
    // ---------------------------

    int CalcEigen(DataSet_Matrix&,int); // TODO: Make const ref
    void PrintModes();
    int EigvalToFreq();
    int MassWtEigvect( const double* );

    int Nmodes() { return nmodes_; }
    double Eigenvalue(int i) { return evalues_[i]; }
    const double* Eigenvector(int i) {
      return evectors_ + (i * vecsize_);
    }
  private:
    Frame avg_;
    double* evalues_;
    double* evectors_;
    int nmodes_;
    int vecsize_;
};
#endif
