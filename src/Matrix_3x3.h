#ifndef INC_MATRIX_3X3_H
#define INC_MATRIX_3X3_H
#include "Vec3.h"
#ifdef MPI
#include "Parallel.h"
#endif
class Matrix_3x3 {
  public:
    Matrix_3x3() {}
    Matrix_3x3(const Matrix_3x3&);
    Matrix_3x3(const double*);
    Matrix_3x3(double);
    Matrix_3x3(double,double,double);
    Matrix_3x3(double m0, double m1, double m2, double m3, double m4,
               double m5, double m6, double m7, double m8)
    {
      M_[0] = m0; M_[1] = m1; M_[2] = m2; M_[3] = m3; M_[4] = m4;
      M_[5] = m5; M_[6] = m6; M_[7] = m7; M_[8] = m8;
    }
    Matrix_3x3& operator=(const Matrix_3x3&);
 
    // NOTE: No bounds check!
    // TODO: Make const ref only?
    double  operator[](int idx) const { return M_[idx]; }
    double& operator[](int idx)       { return M_[idx]; }
    Vec3 Row1() const { return Vec3(M_);   }
    Vec3 Row2() const { return Vec3(M_+3); }
    Vec3 Row3() const { return Vec3(M_+6); }
    Vec3 Col1() const { return Vec3(M_[0], M_[3], M_[6]); }
    Vec3 Col2() const { return Vec3(M_[1], M_[4], M_[7]); }
    Vec3 Col3() const { return Vec3(M_[2], M_[5], M_[8]); }
    void Zero();
    void Print(const char*) const;

    int Diagonalize( Vec3& );
    int Diagonalize_Sort( Vec3& );
    int Diagonalize_Sort_Chirality(Vec3&,int);

    void Transpose();
    /// \return Matrix with rows and columns transposed.
    inline Matrix_3x3 Transposed() const;
    /// \return Result of multiplying this matrix times given 3x3 matrix TODO split into a void and const version
    Matrix_3x3& operator*=(const Matrix_3x3&);
    /// Multiply all elements of this matrix by scalar
    Matrix_3x3& operator*=(double);
    /// \return Result of multiplying this matrix times given scalar
    Matrix_3x3 operator*(double) const;
    
    void RotationAroundZ(double, double);
    void RotationAroundY(double, double);
    void CalcRotationMatrix(Vec3 const&, double);
    void CalcRotationMatrix(double, double, double);
    double RotationAngle();
    Vec3 AxisOfRotation(double);
    /// Multiply 3x3 matrix times double[3]
    void TimesVec(double* result, const double* rhs) const {
      double x = rhs[0];
      double y = rhs[1];
      double z = rhs[2];
      result[0] = ((M_[0]*x) + (M_[1]*y) + (M_[2]*z));
      result[1] = ((M_[3]*x) + (M_[4]*y) + (M_[5]*z));
      result[2] = ((M_[6]*x) + (M_[7]*y) + (M_[8]*z));
    }
    /// Multiply 3x3 matrix times 1x3 vector
    Vec3 operator*(Vec3 const& rhs) const {
      double x = rhs[0]; 
      double y = rhs[1]; 
      double z = rhs[2];
      return Vec3( ((M_[0]*x) + (M_[1]*y) + (M_[2]*z)),
                   ((M_[3]*x) + (M_[4]*y) + (M_[5]*z)),
                   ((M_[6]*x) + (M_[7]*y) + (M_[8]*z))  );
    }
    /// Multiply transpose of 3x3 matrix times double[3]
    void TransposeMult(double* result, const double* rhs) const {
      double x = rhs[0];
      double y = rhs[1];
      double z = rhs[2];
      result[0] = ((M_[0]*x) + (M_[3]*y) + (M_[6]*z));
      result[1] = ((M_[1]*x) + (M_[4]*y) + (M_[7]*z));
      result[2] = ((M_[2]*x) + (M_[5]*y) + (M_[8]*z));
    }
    /// Multiply transpose of 3x3 matrix times 1x3 vector
    Vec3 TransposeMult(Vec3 const& rhs) const {
      double x = rhs[0];
      double y = rhs[1];
      double z = rhs[2];
      return Vec3( ((M_[0]*x) + (M_[3]*y) + (M_[6]*z)),
                   ((M_[1]*x) + (M_[4]*y) + (M_[7]*z)),
                   ((M_[2]*x) + (M_[5]*y) + (M_[8]*z))  );
    }
    /// Multiply this times 3x3 matrix
    Matrix_3x3 operator*(Matrix_3x3 const&) const;
    /// Multiply this times transpose of 3x3 matrix
    Matrix_3x3 TransposeMult(Matrix_3x3 const&) const;
    // TODO: Get rid of this
    const double* Dptr() const { return M_; }
    double* Dptr() { return M_; }
#   ifdef MPI
    void BroadcastMatrix(Parallel::Comm const&);
    int SendMatrix(int, Parallel::Comm const&) const;
    int RecvMatrix(int, Parallel::Comm const&);
#   endif
  private:
    double M_[9];
    // The following three variables are set during Diagonalize_Sort. They
    // indicate the original ordering of the eigenvalues/eigenvectors. This
    // information can be used to prevent reflections when e.g. aligning
    // coordinates along principal axes (see e.g. Action_Principal).
    int i1_;
    int i2_;
    int i3_;
    static const int MAX_ITERATIONS;

    int jacobiCheckChirality();
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
Matrix_3x3 Matrix_3x3::Transposed() const {
  return Matrix_3x3( M_[0], M_[3], M_[6],
                     M_[1], M_[4], M_[7],
                     M_[2], M_[5], M_[8] );
}
#endif
