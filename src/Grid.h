#ifndef INC_GRID_H
#define INC_GRID_H
#include "ArrayIterator.h"
#include <algorithm>
/// Three-dimensional grid template.
template <class T> class Grid {
  public:
    Grid() : nx_(0), ny_(0), nz_(0), nelements_(0), grid_(0) {}
    ~Grid() { if (grid_ != 0) delete[] grid_; }
    Grid( const Grid& );
    Grid& operator=( const Grid& );
    T& operator[](size_t idx)             { return grid_[idx];  }
    const T& operator[](size_t idx) const { return grid_[idx];  }
    /// \ return total number of grid points.
    size_t size()                   const { return nelements_;  }
    /// \return estimated size in bytes for given dimensions.
    static size_t sizeInBytes(size_t, size_t, size_t);
    /// Set up grid for given X, Y, and Z dimensions.
    int resize(size_t,size_t,size_t);
    /// \return Size of X dimension.
    size_t NX() const { return nx_; }
    /// \return Size of Y dimension.
    size_t NY() const { return ny_; }
    /// \return Size of Z dimension.
    size_t NZ() const { return nz_; }
    /// Increment grid point by given value.
    long int incrementBy(size_t, size_t, size_t, const T&);
    /// Set grid point to value.
    void setGrid(size_t, size_t, size_t, const T&);
    /// \return element at a specified grid point.
    const T& element(size_t, size_t, size_t) const;
    /// Convert X, Y, and Z bin #s to index.
    long int CalcIndex(size_t x, size_t y, size_t z) const { 
      return (long int)(x*(ny_*nz_))+(y*nz_)+z;
    }
    // NOTE: This way of calculating overall index is consistent with
    //       how Matrix index is calcd but not used for backwards compat.
    //long int CalcIndex(size_t x, size_t y, size_t z) const { return (z*nx_*ny_)+(y*nx_)+x; }
    /// Convert linear index to X, Y, and Z bin indices.
    void ReverseIndex(long int idx, size_t& x, size_t& y, size_t& z) const {
      //TODO check negative index?
      x = (size_t)idx / (ny_*nz_);
      y = ((size_t)idx / nz_) % ny_;
      z = (size_t)idx % nz_;
    }
    /// Iterator over grid elements.
    typedef ArrayIterator<T> iterator;
    iterator begin() { return grid_;              }
    iterator end()   { return grid_ + nelements_; }
    /// \return memory usage in bytes
    size_t DataSize() const {
      return ( (4*sizeof(size_t)) + (nelements_ * sizeof(T)) );
    }
  private:
    size_t nx_;        ///< Grid X dimension.
    size_t ny_;        ///< Grid Y dimension.
    size_t nz_;        ///< Grid Z dimension.
    size_t nelements_; ///< Total number of grid points.
    T* grid_;          ///< Array of grid points.
};
/** \return estimated size in bytes for given dimensions. */
template <class T> size_t Grid<T>::sizeInBytes(size_t nX, size_t nY, size_t nZ) {
  return ((4*sizeof(size_t)) + sizeof(T*) + (nX*nY*nZ*sizeof(T)));
}
// COPY CONSTRUCTOR
template <class T> Grid<T>::Grid(const Grid& rhs) :
  nx_(rhs.nx_),
  ny_(rhs.ny_),
  nz_(rhs.nz_),
  nelements_(rhs.nelements_),
  grid_(0)
{
  if (nelements_ > 0L) {
    grid_ = new T[ nelements_ ];
    std::copy( rhs.grid_, rhs.grid_ + nelements_, grid_ );
  }
}
// ASSIGNMENT
template <class T> Grid<T>& Grid<T>::operator=(const Grid& rhs) {
  if (this == &rhs) return *this;
  if (grid_!=0) {
    delete[] grid_;
    grid_ = 0;
  }
  nx_ = rhs.nx_;
  ny_ = rhs.ny_;
  nz_ = rhs.nz_;
  nelements_ = rhs.nelements_;
  if (nelements_ > 0L) {
    grid_ = new T[ nelements_ ];
    std::copy( rhs.grid_, rhs.grid_ + nelements_, grid_ );
  }
  return *this;
}
// Grid::resize()
template <class T> int Grid<T>::resize(size_t x, size_t y, size_t z) {
  if (grid_!=0) {
    delete[] grid_;
    grid_ = 0;
  }
  nx_ = x;
  ny_ = y;
  nz_ = z;
  nelements_ = nx_ * ny_ * nz_;
  if (nelements_ > 0L) {
    grid_ = new T[ nelements_ ];
    if (grid_ == 0) return 1;
    std::fill( grid_, grid_ + nelements_, T() );
  }
  return 0;
}
// Grid::incrementBy()
/** \return Index of underlying bin incremented. */
template <class T> long int Grid<T>::incrementBy(size_t x, size_t y, size_t z, const T& eltIn) {
  long int idx = CalcIndex(x,y,z);
  grid_[idx] += eltIn;
  return idx;
}
// Grid::setGrid()
template <class T> void Grid<T>::setGrid(size_t x, size_t y, size_t z, const T& eltIn) {
  grid_[CalcIndex(x,y,z)] = eltIn;
}
// Grid::element()
template <class T> const T& Grid<T>::element(size_t x, size_t y, size_t z) const {
  return grid_[CalcIndex(x,y,z)];
}
#endif
