#ifndef INC_MATRIX_3X3_H
#define INC_MATRIX_3X3_H
class Matrix_3x3 {
  public:
    Matrix_3x3();
    Matrix_3x3(double*);

    int Diagonalize( double*, double* );
    int Diagonalize_Sort(double (&)[3][3], double *);
    void Print(const char*);
  private:
    static const int MAX_ITERATIONS;
    double M_[9];
};
#endif
