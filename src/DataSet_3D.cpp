#include <cmath> // ceil
#include "DataSet_3D.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
DataSet_3D::DataSet_3D() {}

/** DESTRUCTOR */
DataSet_3D::~DataSet_3D() {}

/** COPY CONSTRUCTOR */
DataSet_3D::DataSet_3D(DataSet_3D const& rhs) :
  DataSet(rhs),
  gridBin_(rhs.gridBin_)
{}

/** ASSIGNMENT */
DataSet_3D& DataSet_3D::operator=(DataSet_3D const& rhs)
{
  if (this == &rhs) return *this;
  DataSet::operator=(rhs);
  gridBin_ = rhs.gridBin_;
  return *this;
}

// DataSet_3D::Allocate_N_O_Box()
int DataSet_3D::Allocate_N_O_Box(size_t nx, size_t ny, size_t nz, 
                                 Vec3 const& oxyz, Box const& boxIn)
{
  if (nx == 0 || ny == 0 || nz == 0) {
    mprinterr("Error: One or more grid sizes are 0: %zu %zu %zu\n", nx, ny, nz);
    return 1;
  }
  // Set origin and unit cell params.
  gridBin_.Setup_Sizes_Origin_Box(nx, ny, nz, oxyz, boxIn);
# ifdef DEBUG_DATASET_3D
  gridBin_.PrintDebug("Allocate_N_O_Box");
# endif
  return Allocate3D(nx, ny, nz);
}

// DataSet_3D::Allocate_N_O_D()
int DataSet_3D::Allocate_N_O_D(size_t nx, size_t ny, size_t nz,
                               Vec3 const& oxyz, Vec3 const& dxyz)
{
  if (nx == 0 || ny == 0 || nz == 0) {
    mprinterr("Error: One or more grid sizes are 0: %zu %zu %zu\n", nx, ny, nz);
    return 1;
  }
  // Set origin and spacing, calculate maximum (for binning).
  gridBin_.Setup_Sizes_Origin_Spacing(nx, ny, nz, oxyz, dxyz);
# ifdef DEBUG_DATASET_3D
  gridBin_.PrintDebug("Allocate_N_O_D");
# endif
  return Allocate3D(nx, ny, nz);
}

// Calc_Origin()
/** For even-spaced grids, origin is center - (N/2)*spacing.
  * For odd-spaced grids, origin is center - ((N-1/2)*spacing)+half_spacing
  */
/*
static double Calc_Origin(size_t N, double D) {
  size_t odd = N % 2;
  size_t half = (N - odd) / 2;
  if (odd)
    return -(((double)half * D) + (D * 0.5));
  else
    return -((double)half * D);
}
*/

/** \return Origin coords calculated from given center coords, spacings, and # of bins. */
/*
Vec3 DataSet_3D::calcOriginFromCenter(Vec3 const& cxyz,
                                      double dx, double dy, double dz,
                                      size_t nx, size_t ny, size_t nz)
{
  return Vec3( cxyz[0] + Calc_Origin(nx, dx),
               cxyz[1] + Calc_Origin(ny, dy),
               cxyz[2] + Calc_Origin(nz, dz) );
}*/

// DataSet_3D::Allocate_N_C_D()
int DataSet_3D::Allocate_N_C_D(size_t nx, size_t ny, size_t nz,
                               Vec3 const& cxyz, Vec3 const& dxyz)
{
  // Calculate origin from center coordinates.
  gridBin_.Setup_Sizes_Center_Spacing(nx, ny, nz, cxyz, dxyz);
# ifdef DEBUG_DATASET_3D
  gridBin_.PrintDebug("Allocate_N_C_D");
# endif
  return Allocate3D(nx, ny, nz);
/*  return Allocate_N_O_D(nx, ny, nz,
                        calcOriginFromCenter(cxyz, dxyz[0], dxyz[1], dxyz[2], nx, ny, nz),
                        dxyz);*/
/*
  Vec3 oxyz( cxyz[0] + Calc_Origin(nx, dxyz[0]),
             cxyz[1] + Calc_Origin(ny, dxyz[1]),
             cxyz[2] + Calc_Origin(nz, dxyz[2]) );
  return Allocate_N_O_D(nx,ny,nz,oxyz,dxyz);*/
}

// DataSet_3D::Allocate_X_C_D()
int DataSet_3D::Allocate_X_C_D(Vec3 const& sizes, Vec3 const& center, Vec3 const& dxyz)
{
  // Calculate bin counts
  /*size_t nx = (size_t)ceil(sizes[0] / dxyz[0]);
  size_t ny = (size_t)ceil(sizes[1] / dxyz[1]);
  size_t nz = (size_t)ceil(sizes[2] / dxyz[2]);
  return Allocate_N_C_D( nx, ny, nz, center, dxyz );*/
  GridBin::SizeArray gridSizes = gridBin_.Setup_Lengths_Center_Spacing(sizes, center, dxyz);
# ifdef DEBUG_DATASET_3D
  gridBin_.PrintDebug("Allocate_X_C_D");
# endif
  return Allocate3D(gridSizes[0], gridSizes[1], gridSizes[2]);
}

/** Set the center of the grid to given coords via moving the grid origin. */
/*void DataSet_3D::SetGridCenter(Vec3 const& cxyz) {
  gridBin_.SetOrigin( calcOriginFromCenter(cxyz,
                                           gridBin_.DX(), gridBin_.DY(), gridBin_.DZ(),
                                           NX(), NY(), NZ())
                    );
}*/

// DataSet_3D::GridInfo()
void DataSet_3D::GridInfo() const {
  Vec3 const& oxyz = gridBin_.GridOrigin();
  Vec3 cxyz = gridBin_.GridCenter();
  mprintf("\t\t-=Grid Dims=- %8s %8s %8s\n", "X", "Y", "Z");
  mprintf("\t\t        Bins: %8zu %8zu %8zu\n", NX(), NY(), NZ());
  mprintf("\t\t      Origin: %8g %8g %8g\n", oxyz[0], oxyz[1], oxyz[2]);
  //if (gridBin_.IsOrthoGrid()) {
    mprintf("\t\t     Spacing: %8g %8g %8g\n", gridBin_.DX(), gridBin_.DY(), gridBin_.DZ());
    mprintf("\t\t      Center: %8g %8g %8g\n", cxyz[0], cxyz[1], cxyz[2]);
    //mprintf("\tGrid max    : %8.3f %8.3f %8.3f\n", gridBin_.MX(), gridBin_.MY(), gridBin_.MZ());
  //} else {
    Box const& box = gridBin_.GridBox();
    mprintf("\t\tBox: %s ABC={%g %g %g} abg={%g %g %g}\n", box.CellShapeName(),
            box.Param(Box::X), box.Param(Box::Y), box.Param(Box::Z),
            box.Param(Box::ALPHA), box.Param(Box::BETA), box.Param(Box::GAMMA));
  //}
}

#ifdef MPI
/** Sum grid across ranks to master, ensure master has same orientation as final rank. */
int DataSet_3D::Sync(size_t total, std::vector<int> const& rank_frames, Parallel::Comm const& commIn)
{
  SyncGrid(total, rank_frames, commIn);
  gridBin_.Sync( commIn );
  return 0;
}
#endif
