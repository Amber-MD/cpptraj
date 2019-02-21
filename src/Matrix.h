#ifndef INC_MATRIX_H
#define INC_MATRIX_H
#include "ArrayIterator.h"
#include <cmath> // For linear least squares
/// Two-dimensional matrix template.
template <class T> class Matrix {
// TODO: Type may not be necessary here if in DataSet_2D
    enum MType { FULL = 0, HALF, TRIANGLE }; 
  public:
    /// CONSTRUCTOR
    Matrix() :
      elements_(0), ncols_(0L), nrows_(0L), nelements_(0L),
      maxElements_(0L), currentElement_(0L), type_(FULL),
      indexFxn_(calcFullIndex) {}
    /// DESTRUCTOR
    ~Matrix() { if (elements_!=0) delete[] elements_; }
    /// COPY CONSTRUCTOR
    Matrix( const Matrix& );
    /// ASSIGNMENT
    Matrix& operator=( const Matrix& );
    /// \return Element at absolute index idx
    T& operator[](size_t idx)                  { return elements_[idx];  }
    /// \return Const reference to element at absolute index idx
    const T& operator[](size_t idx)      const { return elements_[idx];  }
    /// \return total number of elements in the matrix.
    size_t size()                        const { return nelements_;      }
    /// \return current matrix type.
    MType Type()                         const { return type_;           }
    /// \return estimated size in bytes.
    static size_t sizeInBytes(size_t,size_t);
    /// Set up matrix for given number of cols and rows.
    int resize(size_t,size_t);
    /// \return element at specified col and row.
    const T& element(size_t,size_t) const;
    T&       element(size_t,size_t);
    /// \return number of rows (Y).
    size_t Nrows()     const { return nrows_;     }
    /// \return number of cols (X).
    size_t Ncols()     const { return ncols_;     }
    /// Add an element to the matrix in order.
    int addElement( const T& );
    /// Set element at col and row.
    void setElement(size_t, size_t, const T&);
    /// Add given value to element at col and row.
    void updateElement(size_t, size_t, const T&);
    /// \return pointer to internal array of elements.
    T const* Ptr()     const { return elements_;  }
    T* Ptr()                 { return elements_;  }
    /// Convert X and Y to index. 
    size_t CalcIndex(size_t x, size_t y) const { return indexFxn_(ncols_, x, y); }
    /// Iterator over matrix elements
    typedef ArrayIterator<T> iterator;
    /// Iterator to beginning of matrix elements.
    iterator begin() { return elements_;              }
    iterator begin() const { return elements_; }
    /// Iterator to end of matrix elements.
    iterator end()   { return elements_ + nelements_; }
    iterator end()   const { return elements_ + nelements_; }
    /// \return memory used by matrix in bytes.
    size_t DataSize() const {
      return (nelements_*sizeof(T)) + sizeof(T) +
             (5 * sizeof(size_t) + sizeof(MType) +
             sizeof(long int(*)()));
    }
    /// Clear matrix
    void clear() {
      if (elements_ != 0) delete[] elements_;
      ncols_ = 0;
      nrows_ = 0;
      nelements_ = 0;
      maxElements_ = 0;
      currentElement_ = 0;
    }
    /// Solve linear least-squares problem
    void LinearLeastSquares(T* b);
  private:
    T* elements_;           ///< Array of elements
    T diagElt_;             ///< For TRIANGLE, the value of the diagonal.
    size_t ncols_;          ///< Number of columns (X)
    size_t nrows_;          ///< Number of rows (Y)
    size_t nelements_;      ///< Total number of elements.
    size_t maxElements_;    ///< Max number of elements currently allocated for.
    size_t currentElement_; ///< Current element (for AddElement())
    MType type_;            ///< Current matrix type.
    /// Pointer to index calculator for current matrix type
    long int (*indexFxn_)(size_t,size_t,size_t);
    /// Full 2D matrix. 
    static long int calcFullIndex(size_t nX, size_t x, size_t y) {
      return (long int)( (y*nX)+x );
    }
    /// Upper-half matrix + diagonal.
    static long int calcHalfIndex(size_t nX, size_t xIn, size_t yIn) {
      size_t i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else {
        i = yIn;
        j = xIn;
      }
      return (long int)((i * nX) - (i * (i-1UL) / 2UL) + (j - i));
    }
    /// Upper-half matrix - diagonal.
    static long int calcTriIndex(size_t nX, size_t xIn, size_t yIn) {
      size_t i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else if (xIn > yIn) {
        i = yIn;
        j = xIn;
      } else // iIn == jIn, triangle matrix diagonal is indicated by -1 
        return -1L;
      size_t i1 = i + 1UL;
      return (long int)(( (nX * i) - ((i1 * i) / 2UL) ) + j - i1);
    }
};
// COPY CONSTRUCTOR
template<class T> Matrix<T>::Matrix(const Matrix& rhs) :
  elements_(0),
  diagElt_( rhs.diagElt_ ),
  ncols_( rhs.ncols_ ),
  nrows_( rhs.nrows_ ),
  nelements_( rhs.nelements_ ),
  maxElements_( rhs.maxElements_ ),
  currentElement_( rhs.currentElement_ ),
  type_( rhs.type_ ),
  indexFxn_( rhs.indexFxn_ )
{
  if (maxElements_ > 0L) {
    elements_ = new T[ maxElements_ ];
    std::copy( rhs.elements_, rhs.elements_ + nelements_, elements_ );
  }
}
// ASSIGNMENT
template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix& rhs) {
  if (this == &rhs) return *this;
  if (elements_!=0) {
    delete[] elements_;
    elements_ = 0;
  }
  ncols_ = rhs.ncols_;
  nrows_ = rhs.nrows_;
  nelements_ = rhs.nelements_;
  maxElements_ = rhs.maxElements_;
  diagElt_ = rhs.diagElt_;
  currentElement_ = rhs.currentElement_;
  type_  = rhs.type_;
  indexFxn_ = rhs.indexFxn_;
  if (maxElements_ > 0L) {
    elements_ = new T[ maxElements_ ];
    std::copy( rhs.elements_, rhs.elements_ + nelements_, elements_ );
  }
  return *this;
}
// Matrix::sizeInBytes()
/** \return Total size of Matrix in bytes. */
template<class T> size_t Matrix<T>::sizeInBytes(size_t nX, size_t nY) {
  size_t elSize;
  if (nX > 0L && nY > 0L)       // FULL
    elSize = nX * nY * sizeof( T );
  else if (nX > 0L && nY == 0L) // HALF
    elSize = (nX * (nX + 1UL) / 2UL) * sizeof( T );
  else if (nX == 0L && nY > 0L) // TRIANGLE
    elSize = (nY * (nY - 1UL) / 2UL) * sizeof ( T );
  else elSize = 0UL;
  return (elSize + sizeof(T*) + sizeof(T) + // elements, pointer, diagEl
          (5 * sizeof(size_t)) + sizeof(MType) + sizeof(long int(*)()));
}
// Matrix::resize()
/** Reallocate matrix with given number of columns and rows.
  * \param nX number of columns.
  * \param nY number of rows.
  * Three types of allocation:
  *   - cols > 0,  rows > 0  : Full matrix
  *   - cols > 0,  rows == 0 : Half matrix
  *   - cols == 0, rows > 0  : Triangle matrix
  * \return 0 if allocation successful, 1 otherwise.
  */
template<class T> int Matrix<T>::resize(size_t nX, size_t nY) {
  diagElt_ = T(); // Diagonal element default to zero.
  if (nX > 0L && nY > 0L) { // FULL
    ncols_ = nX;
    nrows_ = nY;
    nelements_ = ncols_ * nrows_;
    indexFxn_ = calcFullIndex;
    type_ = FULL;
  } else if (nX > 0L && nY == 0L) { // HALF
    ncols_ = nX;
    nrows_ = nX;
    nelements_ = ncols_ * (ncols_ + 1L) / 2L;
    indexFxn_ = calcHalfIndex;
    type_ = HALF;
  } else if (nX == 0L && nY > 0L) { // TRIANGLE
    ncols_ = nY;
    nrows_ = nY;
    nelements_ = (ncols_ * (ncols_ - 1L) / 2L);
    indexFxn_ = calcTriIndex;
    type_ = TRIANGLE;
  } else { // Both Zero, EMPTY
    ncols_ = 0L;
    nrows_ = 0L;
    nelements_ = 0L;
    return 1;
  }
  currentElement_ = 0L;
  if (nelements_ > 0L) {
    if (nelements_ > maxElements_) {
      if (elements_ != 0) delete[] elements_;
      elements_ = new T[ nelements_ ];
      if (elements_ == 0) return 1;
      maxElements_ = nelements_;
    }
    std::fill(elements_, elements_ + nelements_, T());
  }
  return 0;
}
// Matrix::addElement()
/** Add the input T to the element array and increment currentElement.
  * \return 1 on success, 0 if no more elements can be added.
  */
template<class T> int Matrix<T>::addElement(const T& elementIn) {
  if (currentElement_>=nelements_) return 0;
  elements_[currentElement_] = elementIn;
  ++currentElement_;
  return 1;
}
// Matrix::setElement()
/** Set element at specified column xIn and row yIn. */
template<class T> void Matrix<T>::setElement(size_t xIn, size_t yIn, const T& eltIn) {
  long int idx = indexFxn_(ncols_, xIn, yIn);
  elements_[idx] = eltIn;
}
// Matrix::updateElement()
/** Update specified element at column xIn and row yIn via += */
template<class T> void Matrix<T>::updateElement(size_t xIn, size_t yIn, const T& eltIn) {
  long int idx = indexFxn_(ncols_, xIn, yIn);
  elements_[idx] += eltIn;
}
// Matrix::element()
/** \return Constant reference to specified element at column xIn and row yIn. */
template<class T> const T& Matrix<T>::element(size_t xIn, size_t yIn) const {
  long int idx = indexFxn_(ncols_, xIn, yIn);
  if (idx < 0) return diagElt_; // In case of xIn == yIn for TRIANGLE
  return elements_[idx];
}
/** \return Reference to specified element at column xIn and row yIn. */
template<class T> T& Matrix<T>::element(size_t xIn, size_t yIn) {
  long int idx = indexFxn_(ncols_, xIn, yIn);
  if (idx < 0) return diagElt_; // In case of xIn == yIn for TRIANGLE
  return elements_[idx];
}
// Matrix::LinearLeastSquares()
/** Produce a vector containing the solution to an over-determined system of linear equations
  * represented by the matrix itself and an input vector.  This will work for float and
  * double-precision matrices.  The input vector will be overwritten with the solution, and
  * the matrix itself will be turned into an upper-triangular "R" result of the QR
  * decomposition. */
template<class T> void Matrix<T>::LinearLeastSquares(T* b)
{
  unsigned int i, j, k;
  double tnm_v, tnm_v2, tempval, sign_v;
  double multval, temp_b, pivot;

  // Eliminate zero rows
  k = 0;
  for (i = 0; i < nrows_; i++) {
    double msum = 0.0;
    for (j = 0; j < ncols_; j++) {
      msum += fabs(elements_[i*ncols_ + j]);
    }
    if (msum >= 1.0e-6) {
      for (j = 0; j < ncols_; j++) {
        elements_[k*ncols_ + j] = elements_[i*ncols_ + j];
      }
      b[k] = b[i];
      k++;
    }
  }
  if (k == 0 || k < ncols_) {
    for (i = 0; i < nrows_; i++) {
      b[i] = 0.0;
    }
    return;
  }
  int nEq = k;

  // Eliminate zero columns, but remember the original order of variables
  // so that variables correspnding to zero columns will get set to zero
  k = 0;
  int* solnID;
  solnID = new int[ncols_];
  for (i = 0; i < ncols_; i++) {
    double msum = 0.0;
    for (j = 0; j < nrows_; j++) {
      msum += fabs(elements_[j*ncols_ + i]);
    }
    if (msum >= 1.0e-6) {
      for (j = 0; j < nrows_; j++) {
	elements_[j*ncols_ + k] = elements_[j*ncols_ + i];
      }
      solnID[i] = k;
      k++;
    }
    else {
      solnID[i] = -1;
    }
  }
  unsigned int nVar = k;
  
  // Scratch arrays
  double* v;
  double* vprime;
  v = new double[nEq];
  vprime = new double[nEq];
  double *tmp;

  // Loop over columns
  for (k = 0; k < nVar; k++) {

    // Compute the kth column of Q*
    tnm_v2 = 0.0;
    for (i = 0; i < nEq - k; i++) {
      v[i] = elements_[(i + k)*ncols_ + k];
      tnm_v2 += v[i] * v[i];
    }
    sign_v = (v[0] >= 0.0) ? 1.0 : -1.0;
    tnm_v = sqrt(tnm_v2);
    tnm_v2 -= v[0] * v[0];
    v[0] += sign_v * tnm_v;
    tnm_v = 1.0/sqrt(tnm_v2 + v[0]*v[0]);
    for (i = 0; i < nEq - k; i++) {
      v[i] = v[i] * tnm_v;
    }

    // Update A as R evolves
    for (i = 0; i < nVar - k; i++) {
      vprime[i] = 0.0;
    }
    for (i = 0; i < nEq - k; i++) {
      tmp = &elements_[(i + k)*ncols_ + k];
      tempval = v[i];
      for (j = 0; j < nVar - k; j++) {
	vprime[j] += tempval * tmp[j];
      }
    }
    for (i = 0; i < nEq - k; i++) {
      tmp = &elements_[(i + k)*ncols_ + k];
      tempval = 2.0 * v[i];
      for (j = 0; j < nVar - k; j++) {
	tmp[j] -= tempval * vprime[j];
      }
    }

    // Update b as Q* evolves
    tmp = &b[k];
    tempval = 0.0;
    for (i = 0; i < nEq - k; i++) {
      tempval += v[i] * tmp[i];
    }
    tempval *= 2.0;
    for (i = 0; i < nEq - k; i++) {
      tmp[i] -= tempval * v[i];
    }
  }

  // Back-substitution for the solution
  for (i = nVar - 1; i > 0; i--) {
    pivot = 1.0 / elements_[i*ncols_ + i];
    temp_b = b[i];
    for (j = i-1; j >= 0; j--) {
      multval = elements_[j*ncols_ + i] * pivot;
      b[j] -= multval * temp_b;
    }
    b[i] *= pivot;
  }
  b[0] /= elements_[0];

  // Re-order the solved variables, placing zeroes if any columns were skipped
  for (i = ncols_ - 1; i >= 0; i--) {
    b[i] = (solnID[i] >= 0) ? b[solnID[i]] : 0.0;
  }
  
  // Free allocated memory
  delete[] solnID;
  delete[] v;
  delete[] vprime;
}
#endif
