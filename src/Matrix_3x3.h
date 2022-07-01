#ifndef INC_MATRIX_3X3_H
#define INC_MATRIX_3X3_H
#include "Vec3.h"
#ifdef MPI
#include "Parallel.h"
#endif
class Matrix_3x3 {
  public:
    /// CONSTRUCTOR
    Matrix_3x3() {}
    /// COPY CONSTRUCTOR
    Matrix_3x3(const Matrix_3x3&);
    /// CONSTRUCTOR - pointer to flattened matrix array of size 9, assumes row-major
    Matrix_3x3(const double*);
    /// CONSTRUCTOR - Set all elements to a single number
    Matrix_3x3(double);
    /// CONSTRUCTOR - Set diagonal, all others to zero
    Matrix_3x3(double,double,double);
    /// CONSTRUCTOR - Set all elements individually (row major)
    Matrix_3x3(double m0, double m1, double m2, double m3, double m4,
               double m5, double m6, double m7, double m8)
    {
      M_[0] = m0; M_[1] = m1; M_[2] = m2; M_[3] = m3; M_[4] = m4;
      M_[5] = m5; M_[6] = m6; M_[7] = m7; M_[8] = m8;
    }
    /// ASSIGNMENT
    Matrix_3x3& operator=(const Matrix_3x3&);
 
    // NOTE: No bounds check!
    // TODO: Make const ref only?
    /// \return Element at specified index
    double  operator[](int idx) const { return M_[idx]; }
    /// \return Element at specified index
    double& operator[](int idx)       { return M_[idx]; }
    /// \return First row
    Vec3 Row1() const { return Vec3(M_);   }
    /// \return Second row
    Vec3 Row2() const { return Vec3(M_+3); }
    /// \return Third row
    Vec3 Row3() const { return Vec3(M_+6); }
    /// \return First column
    Vec3 Col1() const { return Vec3(M_[0], M_[3], M_[6]); }
    /// \return Second column
    Vec3 Col2() const { return Vec3(M_[1], M_[4], M_[7]); }
    /// \return Third column
    Vec3 Col3() const { return Vec3(M_[2], M_[5], M_[8]); }
    /// \return const pointer to internal matrix data
    const double* Dptr() const { return M_; }
    /// \return pointer to internal matrix data
    double* Dptr() { return M_; }

    /// Zero all elements of the matrix
    void Zero();
    /// Print matrix to stdout
    void Print(const char*) const;

    /// Diagonalize matrix, store eigenvectors in columns, set eigenvalues
    int Diagonalize( Vec3& );
    /// Diagonalize matrix and sort eigenvectors (in cols) by eigenvalue
    int Diagonalize_Sort( Vec3& );
    /// Diagonalize matrix, sort eignvectors, attempt to correct for eigenvector sign flips
    int Diagonalize_Sort_Chirality(Vec3&,int);

    /// Transpose the matrix
    void Transpose();
    /// \return Matrix with rows and columns transposed.
    inline Matrix_3x3 Transposed() const;

    /// \return Result of multiplying this matrix times given 3x3 matrix TODO split into a void and const version
    Matrix_3x3& operator*=(const Matrix_3x3&);
    /// Multiply all elements of this matrix by scalar
    Matrix_3x3& operator*=(double);
    /// \return Result of multiplying this matrix times given scalar
    Matrix_3x3 operator*(double) const;
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

    // ----- Rotation functions ------------------
    /// Calculate rotation matrix around given axis of given magnitude
    void CalcRotationMatrix(Vec3 const&, double);
    /// Calculate rotation matrix around X, Y, and Z axes
    void CalcRotationMatrix(double, double, double);
    /// \return Angle of rotation from matrix
    double RotationAngle() const;
    /// Decompose rotation matrix into Euler angles around each axis
    int RotationAngles(double&, double&, double&) const;
    /// Given theta, extract axis of rotation from rotation matrix
    Vec3 AxisOfRotation(double) const;

#   ifdef MPI
    /// Broadcast from master to other ranks
    void BroadcastMatrix(Parallel::Comm const&);
    /// Send matrix to rank
    int SendMatrix(int, Parallel::Comm const&) const;
    /// Receive matrix from rank
    int RecvMatrix(int, Parallel::Comm const&);
#   endif
  private:
    /// Calculate rotation matrix around Z axis.
    void RotationAroundZ(double, double);
    /// Calculate rotation matrix around Y axis.
    void RotationAroundY(double, double);
    /// Try to fix eigenvector sign flips (Diagonalize_Sort_Chirality())
    int jacobiCheckChirality();

    double M_[9];
    // The following three variables are set during Diagonalize_Sort. They
    // indicate the original ordering of the eigenvalues/eigenvectors. This
    // information can be used to prevent reflections when e.g. aligning
    // coordinates along principal axes (see e.g. Action_Principal).
    int i1_;
    int i2_;
    int i3_;
    static const int MAX_ITERATIONS_; ///< Max iterations to use in Jacobi
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
Matrix_3x3 Matrix_3x3::Transposed() const {
  return Matrix_3x3( M_[0], M_[3], M_[6],
                     M_[1], M_[4], M_[7],
                     M_[2], M_[5], M_[8] );
}
#endif
