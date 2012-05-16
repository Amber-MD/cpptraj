#ifndef INC_GRID_H
#define INC_GRID_H
#include <string>
#include "ArgList.h"
//#include "AtomMask.h"
#include "Topology.h"
#include "CpptrajFile.h"
class Grid {
  public:
    Grid();
    ~Grid();

    int GridInit(const char*, ArgList&);
    void GridInfo();
    int GridAllocate();
    int GridSetup(Topology*);
    //void GridPrint(bool, double, double);
    void GridPrintHeader(CpptrajFile&);

    int SX() { return sx_; }
    int SY() { return sy_; }
    int SZ() { return sz_; }
    int NX() { return nx_; }
    int NY() { return ny_; }
    int NZ() { return nz_; }
    bool GridBox() { return box_; }

    double GridVal(int i, int j, int k) {
      return grid_[i*ny_*nz_ + j*nz_ + k];
    }

    /** Main grid routine. */
    // NOTE: Placed in header for speed - does it actually inline though?
    void GridPoint( double xx, double yy, double zz) {
      int i = (int) (xx / dx_) - 1;
      int j = (int) (xx / dx_) - 1;
      int k = (int) (xx / dx_) - 1;
      // TODO: Benchmark this against multiple comparisons to individual bounds
      int idx = i*ny_*nz_ + j*nz_ + k;
      if (idx >= 0 && idx < gridsize_)
        grid_[idx] += increment_;
    }
  private:
    //double max_;
    int increment_; ///< Set to -1 if negative, 1 if not.
    bool box_;
    //bool invert_;
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
    //int frames_;
    float* grid_;
    //float* dipolex_;
    //float* dipoley_;
    //float* dipolez_;
    //std::string filename_;
    std::string callingRoutine_;
    //AtomMask gridmask_;
};
#endif
