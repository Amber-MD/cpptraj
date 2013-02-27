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
    /// Indicate which kind of gridding to perform
    enum GridModeType { ORIGIN = 0, BOX, CENTER };
    static void Help();
    GridModeType GridMode()      { return mode_;       } // TODO: Obsolete
    AtomMask const& CenterMask() { return centerMask_; } // TODO: Obsolete
    /// Initialize grid from argument list.
    int GridInit(const char*, ArgList&);
    /// Initialize grid from Density file
    int GridInit(std::string const&, const char*);
    /// Print information about the grid, allocate memory.
    void GridInfo();
    /// Setup grid based on given topology.
    int GridSetup(Topology*);
    /// Grid the given XYZ point.
    int GridPoint(Vec3 const&);
    /// Grid the given frame
    void GridFrame(Frame& currentFrame, AtomMask const& mask);
    /// Grid point (for backwards compat. with Action_Dipole) 
    int BinPoint(double, double, double); 
    /// Print Xplor format grid density
    void PrintXplor(std::string const&, const char*, std::string);
    void PrintPDB(std::string const&, double, double);
    /// Print DX format grid density
    void PrintDX(std::string const&);
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
