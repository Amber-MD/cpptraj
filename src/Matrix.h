#ifndef INC_MATRIX_H
#define INC_MATRIX_H
#include <cstring> // memset, memcpy
/// Two-dimensional matrix template.
template <class T> class Matrix {
  public:
    enum MType { FULL = 0, HALF, TRI }; 
    Matrix() :
      elements_(0), ncols_(0L), nrows_(0L), nelements_(0L), 
      currentElement_(0L), type_(FULL), calcIndex(calcFullIndex) {}
    ~Matrix() { if (elements_!=0) delete[] elements_; }
    Matrix( const Matrix& );
    Matrix& operator=( const Matrix& );
    T& operator[](size_t idx)                  { return elements_[idx];  }
    const T& operator[](size_t idx)      const { return elements_[idx];  }
    /// \return total number of elements in the matrix.
    size_t size()                        const { return nelements_;      }
    /// \return current matrix type.
    MType Type()                         const { return type_;           }
    /// Set up matrix for given number of cols and rows.
    int resize(size_t,size_t);
    /// \return element at specified col and row.
    const T& element(size_t,size_t) const;
    /// \return number of rows (Y).
    size_t Nrows()     const { return nrows_;     }
    /// \return number of cols (X).
    size_t Ncols()     const { return ncols_;     }
    /// Add an element to the matrix in order.
    int addElement( const T& );
    /// Set element at col and row.
    void setElement(size_t,size_t, const T&);
    /// \return pointer to internal array of elements.
    T const* Ptr()     const { return elements_;  }
    T* Ptr()                 { return elements_;  }
    // DEBUG
    size_t CalcIndex(size_t x, size_t y) { return calcIndex(ncols_, x, y); }
  private:
    T* elements_;           ///< Array of elements
    size_t ncols_;          ///< Number of columns (X)
    size_t nrows_;          ///< Number of rows (Y)
    size_t nelements_;      ///< Total number of elements.
    size_t currentElement_; ///< Current element (for AddElement())
    MType type_;            ///< Current matrix type.
    /// Pointer to index calculator for current matrix type
    size_t (*calcIndex)(size_t,size_t,size_t);
    /// Full 2D matrix. 
    static size_t calcFullIndex(size_t nX,size_t x,size_t y) { return ( (y*nX)+x ); }
    /// Upper-half matrix + diagonal.
    static size_t calcHalfIndex(size_t nX, size_t xIn, size_t yIn) {
      size_t i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else {
        i = yIn;
        j = xIn;
      }
      return (i * nX - (i * (i-1L) / 2L) + (j - i));
    }
    /// Upper-half matrix - diagonal.
    static size_t calcTriIndex(size_t nX, size_t xIn, size_t yIn) {
      size_t i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else if (xIn > yIn) {
        i = yIn;
        j = xIn;
      } else // iIn == jIn, invalid for triangle matrix
        return 0L;
      size_t i1 = i + 1L;
      return ( ( (nX * i) - ((i1 * i) / 2L) ) + j - i1 );
    }
};
// COPY CONSTRUCTOR
template<class T> Matrix<T>::Matrix(const Matrix& rhs) :
  elements_(0),
  ncols_( rhs.ncols_ ),
  nrows_( rhs.nrows_ ),
  nelements_( rhs.nelements_ ),
  currentElement_( rhs.currentElement_ )
{
  if (nelements_ > 0L) {
    elements_ = new T[ nelements_ ];
    memcpy( elements_, rhs.elements_, nelements_*sizeof(T) );
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
  currentElement_ = rhs.currentElement_;
  if (nelements_ > 0L) {
    elements_ = new T[ nelements_ ];
    memcpy( elements_, rhs.elements_, nelements_*sizeof(T) );
  }
  return *this;
}
// Matrix::resize()
/** \param nX number of columns.
  * \param nY number of rows.
  * Three types of allocation:
  *   - cols > 0, rows > 0: Full matrix
  *   - cols > 0, rows == 0: Half matrix
  *   - cols == 0, rows > 0: Triangle matrix
  */
template<class T> int Matrix<T>::resize(size_t nX, size_t nY) {
  if (elements_!=0) {
    delete[] elements_;
    elements_ = 0;
  }
  if (nX > 0L && nY > 0L) { // FULL
    ncols_ = nX;
    nrows_ = nY;
    nelements_ = ncols_ * nrows_;
    calcIndex = calcFullIndex;
    type_ = FULL;
  } else if (nX > 0L && nY == 0L) { // HALF
    ncols_ = nX;
    nrows_ = nX;
    nelements_ = ncols_ * (ncols_ + 1L) / 2L;
    calcIndex = calcHalfIndex;
    type_ = HALF;
  } else if (nX == 0L && nY > 0L) { // TRIANGLE
    ncols_ = nY;
    nrows_ = nY;
    nelements_ = ncols_ * (ncols_ - 1L) / 2L;
    calcIndex = calcTriIndex;
    type_ = TRI;
  } else { // Both Zero, EMPTY
    return 1;
  }
  currentElement_ = 0L;
  elements_ = new T[ nelements_ ];
  memset(elements_, 0, nelements_*sizeof(T) );
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
template<class T> void Matrix<T>::setElement(size_t xIn, size_t yIn, const T& eltIn) {
  size_t idx = calcIndex(ncols_, xIn, yIn);
  elements_[idx] = eltIn;
}
// Matrix::element()
template<class T> const T& Matrix<T>::element(size_t xIn, size_t yIn) const {
  size_t idx = calcIndex(ncols_, xIn, yIn);
  return elements_[idx];
}
#endif
