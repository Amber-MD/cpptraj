#ifndef INC_MATRIX_2D_H
#define INC_MATRIX_2D_H
#include "DataSet.h"
class Matrix_2D : public DataSet {
  public:
    Matrix_2D();
    ~Matrix_2D();
    Matrix_2D( const Matrix_2D& );
    Matrix_2D& operator=( const Matrix_2D& );
    double& operator[](int idx) { return elements_[idx]; }

    int Setup(int,int);
    int AddElement(double);
    void SetElement(int,int,double);
    double GetElement(int,int);
    int Nrows() const { return nrows_; }
    int Ncols() const { return ncols_; }
    int Nelements() const { return (int)nelements_; }

    int Xmax() { return nrows_ - 1; }
    int Size() { return (int)nelements_; }
    void Write2D( CpptrajFile&, int, int);
    void GetDimensions(std::vector<int>&);
  private:
    double* elements_;
    int ncols_;
    int nrows_;
    size_t nelements_;
    size_t currentElement_;

    int calcIndex(int i, int j) { return ( (i*ncols_)+j ); }
};
#endif
