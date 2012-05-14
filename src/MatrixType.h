#ifndef INC_MATRIXTYPE_H
#define INC_MATRIXTYPE_H
#include "Action.h"
#include "VectorType.h"
class MatrixType : public DataSet, public Action {
  public:
    enum matrixMode {
      MATRIX_NULL=0, MATRIX_DIST,      MATRIX_COVAR, MATRIX_MWCOVAR,
      MATRIX_CORREL, MATRIX_DISTCOVAR, MATRIX_IDEA,  MATRIX_IRED
    };
    static const char MatrixOutput[][10];

    MatrixType();
    ~MatrixType();

    void print();

    matrixMode Type()   { return type_; }
    int Mask1Tot()      { return mask1_.Nselected(); }
    int Snap()          { return snap_; }
    /// For interfacing with fortran routines
    // NOTE: Should MatrixPtr and VectPtr return copies?
    double* MatrixPtr() { return mat_; }
    double* VectPtr();
    const AtomMask& Mask1() { return mask1_; }
    Topology* Parm()    { return matrixParm_;}

  private:
    matrixMode type_;
    double* vect_;
    double* vect2_;
    double* mat_;
    int vectsize_;
    int matsize_;
    AtomMask mask1_;
    char* mask2expr_;
    AtomMask mask2_;
    int mask1tot_;
    int mask2tot_;
    Topology* matrixParm_;
    int Nelt_;
    int snap_;
    // IRED only
    std::vector<VectorType*> IredVectors_;

    // Only needed by action
    static const char MatrixModeString[][27];
    std::string filename_;
    int start_;
    int stop_;
    int offset_;
    int order_;
    enum OutputType { BYATOM=0, BYRESIDUE, BYMASK };
    OutputType outtype_;

    int init();
    int setup();
    int action();

    void Info();
    int HalfMatrixIndex(int, int);
};
#endif
