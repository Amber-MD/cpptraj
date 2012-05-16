#ifndef INC_GRID_H
#define INC_GRID_H
#include <string>
#include "ArgList.h"
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
    void GridPrintHeader(CpptrajFile&);
    // DEBUG
    void PrintEntireGrid();

    int NX() { return nx_; }
    int NY() { return ny_; }
    int NZ() { return nz_; }
    bool GridBox() { return box_; }

    double GridVal(int i, int j, int k) {
      return grid_[i*ny_*nz_ + j*nz_ + k];
    }

    /** Main grid routine. */
    // NOTE: Placed in header for speed - does it actually inline though?
    void GridPoint( double xIn, double yIn, double zIn) {
      double xx = xIn + sx_;
      int i = (int) (xx / dx_) - 1;
      if (i >= 0 && i < nx_) {
        double yy = yIn + sy_;
        int j = (int) (yy / dy_) - 1;
        if (j >= 0 && j < ny_) {
          double zz = zIn + sz_;
          int k = (int) (zz / dz_) - 1;
          if (k >= 0 && k < nz_) {
            int idx = i*ny_*nz_ + j*nz_ + k;

            grid_[idx] += increment_;
          }
        }
      }
    }
  private:
    int increment_; ///< Set to -1 if negative, 1 if not.
    bool box_;
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
