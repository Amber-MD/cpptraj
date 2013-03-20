#ifndef INC_GRID_H
#define INC_GRID_H
#include <string>
#include "ArgList.h"
#include "Topology.h"
/// Hold a 3-dimensional grid.
class Grid {
  public:
    Grid();
    ~Grid();
    Grid(const Grid&);
    Grid& operator=(const Grid&);
    /// Divide all elements by a constant -- no protection for div. by zero
    void operator/=(const double fac) {
      for (int i = 0; i < gridsize_; i++)
        grid_[i] /= fac;
    }
    /// Indicate which kind of gridding to perform
    enum GridModeType { ORIGIN = 0, BOX, CENTER };
    static const char* HelpText;
    GridModeType GridMode()      { return mode_;       } // TODO: Obsolete
    AtomMask const& CenterMask() { return centerMask_; } // TODO: Obsolete
    /// Initialize grid from argument list.
    int GridInit(const char*, ArgList&);
    /// Initialize the grid from a given size and resolution
    int GridInitSizeRes(const char*, double[3], double[3], std::string const&);
    /// Initialize grid from Density file
    int InitFromFile(std::string const&, std::string const&);
    /// Print information about the grid, allocate memory.
    void GridInfo();
    /// Setup grid based on given topology.
    int GridSetup(Topology const&);
    /// Grid the given XYZ point.
    int GridPoint(Vec3 const&);
    /// Add the given value to the density
    void AddDensity(int i, int j, int k, double val) {
      if (i > 0 && j > 0 && k > 0) {
        int idx = i * ny_ * nz_ + j * nz_ + k;
        if (idx < gridsize_) grid_[i*ny_*nz_+j*nz_+k] += val;
      }
    }
    /// Grid the given frame
    void GridFrame(Frame& currentFrame, AtomMask const& mask);
    /// Grid point (for backwards compat. with Action_Dipole)
    int BinPoint(double, double, double); 
    /// Print Xplor format grid density
    void PrintXplor(std::string const&, const char*, std::string);
    void PrintPDB(std::string const&, double, double);
    /// Print DX format grid density w/ either a customizable or default origin
    void PrintDX(std::string const&);
    void PrintDX(std::string const&, double, double, double);
    // DEBUG
    void PrintEntireGrid();
    /// Return number of bins in the X dimension. 
    int NX() { return nx_; }
    /// Return number of bins in the Y dimension.
    int NY() { return ny_; }
    /// Return number of bins in the Z dimension.
    int NZ() { return nz_; }
    /// Return size of the grid
    int GridSize() { return gridsize_; }
    /// Real X coordinate of grid center.
    double SX() { return sx_; }
    /// Real Y coordinate of grid center.
    double SY() { return sy_; }
    /// Real Z coordinate of grid center.
    double SZ() { return sz_; }
    /// Return X coordinate of bin center
    double Xcrd(int i) { return (double)i*dx_ - sx_ + 0.5*dx_; }
    /// Return Y coordinate of bin center
    double Ycrd(int j) { return (double)j*dy_ - sy_ + 0.5*dy_; }
    /// Return Z coordinate of bin center
    double Zcrd(int k) { return (double)k*dz_ - sz_ + 0.5*dz_; }
    /// Return X coordinate of bin corner
    double Xbin(int i) { return (double)i*dx_ - sx_; }
    /// Return Y coordinate of bin corner
    double Ybin(int j) { return (double)j*dy_ - sy_; }
    /// Return Z coordinate of bin corner
    double Zbin(int k) { return (double)k*dz_ - sz_; }
    /// Return population of bin at i, j, k 
    double GridVal(int i, int j, int k) {
      return (double)grid_[i*ny_*nz_ + j*nz_ + k];
    }
    /// Allows you to assign a specific value to a particular grid point
    void SetGridVal(int i, int j, int k, float val) {
      grid_[i*ny_*nz_ + j*nz_ + k] = val;
    }
    void SetGridVal(int i, int j, int k, double val) {
      grid_[i*ny_*nz_ + j*nz_ + k] = (float) val;
    }
    void SetGridVal(int i, float val) {
      grid_[i] = val;
    }
    void SetGridVal(int i, double val) {
      grid_[i] = (float) val;
    }
    /** Returns a grid with 0 everywhere except local maxima in density,
      * filtering out all points whose total density is less than the filter
      * \returns Grid instance with 0s everywhere except local maxima, which are unchanged
      */
    Grid& ExtractPeaks(double min_filter);

    // Iterator over grid
    class iterator : public std::iterator<std::forward_iterator_tag, float> 
    {
      public:
        iterator() : ptr_(0) {}
        iterator(const iterator& rhs) : ptr_(rhs.ptr_) {}
        iterator(float* pin) : ptr_(pin) {}
        iterator& operator=(const iterator& rhs) {
          if (this == &rhs) return *this;
          ptr_ = rhs.ptr_;
          return *this;
        }
        // Relations
        bool operator==(const iterator& rhs) { return (ptr_==rhs.ptr_);}
        bool operator!=(const iterator& rhs) { return (ptr_!=rhs.ptr_);}
        // Increment
        iterator& operator++() { 
          ++ptr_; 
          return *this;
        }
        iterator operator++(int) {
          iterator tmp(*this);
          ++(*this);
          return tmp;
        }
        // Value
        float& operator*() { return *ptr_; }
        // Address
        float* operator->() { return ptr_; }
      private:
        float* ptr_;
    };
    iterator begin() { return grid_; }
    iterator end() { return (grid_ + gridsize_); }
  private:
    int Allocate();

    float increment_; ///< Set to -1 if negative, 1 if not.
    GridModeType mode_;
    AtomMask centerMask_;
    double dx_;
    double dy_;
    double dz_;
    double sx_;
    double sy_;
    double sz_;
    int gridsize_;
    int nx_;
    int ny_;
    int nz_;
    float* grid_;
    std::string callingRoutine_;
};
#endif
