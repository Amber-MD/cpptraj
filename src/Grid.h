#ifndef INC_GRID_H
#define INC_GRID_H
#include <cstring> // memset, memcpy
#include "ArrayIterator.h"
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
    /// Set up grid for given X, Y, and Z dimensions.
    int resize(size_t,size_t,size_t);
    /// \return element at a specified grid point.
    const T& element(size_t,size_t,size_t) const;
    /// \return Size of X dimension.
    size_t NX() const { return nx_; }
    /// \return Size of Y dimension.
    size_t NY() const { return ny_; }
    /// \return Size of Z dimension.
    size_t NZ() const { return nz_; }
    /// Increment grid point.
    void increment(size_t,size_t,size_t);
    /// Increment grid point by given value.
    void incrementBy(size_t,size_t,size_t, const T&);
    /// Set grid point.
    void setGrid(size_t,size_t,size_t, const T&);
    /// Convert X, Y, and Z to index.
    size_t CalcIndex(size_t x, size_t y, size_t z) const { return (x*ny_*nz_)+(y*nz_)+z; }
    /// Iterator over grid elements.
    typedef ArrayIterator<T> iterator;
    iterator begin() { return grid_;              }
    iterator end()   { return grid_ + nelements_; }
  private:
    size_t nx_;        ///< Grid X dimension.
    size_t ny_;        ///< Grid Y dimension.
    size_t nz_;        ///< Grid Z dimension.
    size_t nelements_; ///< Total number of grid points.
    T* grid_;          ///< Array of grid points.
};
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
    memcpy( grid_, rhs.grid_, nelements_ * sizeof(T) );
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
    memcpy( grid_, rhs.grid_, nelements_ * sizeof(T) );
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
    memset(grid_, 0, nelements_ * sizeof(T));
  }
  return 0;
}
// Grid::incrementBy()
template <class T> void Grid<T>::incrementBy(size_t x, size_t y, size_t z, const T& eltIn) {
  size_t idx = CalcIndex(x,y,z);
  grid_[idx] += eltIn;
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
