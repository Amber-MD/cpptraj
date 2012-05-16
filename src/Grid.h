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

    /// Return population of bin at i, j, k 
    double GridVal(int i, int j, int k) {
      return (double)grid_[i*ny_*nz_ + j*nz_ + k];
    }

    /// Return X coordinate of bin
    double Xcrd(int i) { return (double)i*dx_ - nx_*dx_/2.0 + 0.5 * dx_; }
    /// Return Y coordinate of bin
    double Ycrd(int j) { return (double)j*dy_ - ny_*dy_/2.0 + 0.5 * dy_; }
    /// Return Z coordinate of bin
    double Zcrd(int k) { return (double)k*dz_ - nz_*dz_/2.0 + 0.5 * dz_; }
    

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
            grid_[idx] += increment_; // NOTE: Cast to float?
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
