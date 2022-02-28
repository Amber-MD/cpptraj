#ifndef GIST_ENTROPY_UTILS_H
#define GIST_ENTROPY_UTILS_H

#include <cstdlib>
#include <vector>
#include "Vec3.h"
#include <cmath>

namespace GistEntropyUtils {
  typedef std::vector<float> Farray;

  const double GIST_HUGE = 10000.0;

  /** Compute translational and 6D distances to elements in V_XYZ_vec and V_Q_vec, store the smallest in NNd and NNs;
    *
    * For each solvent molecule defined by three elements of V_XYZ_vec (its
    * position) and four elements of V_Q_vec (its quaternion), computes the 3D
    * distance to X, Y, Z, as well as the 6D distance to X, Y, Z, W4, X4, Y4, Z4
    * (also position + quaternion). If the squared 3D distance is < NNd, replaces
    * the current value. If the squared 6D distance is < NNs, replaces the
    * current value.
    *
    * \param center water X, Y, Z
    * \param W4 water quaternion W
    * \param X4 water quaternion X
    * \param Y4 water quaternion Y
    * \param Z4 water quaternion Z
    * \param V_XYZ_vec positions of solvent molecules
    * \param V_Q_vec quaternions of solvent molecules
    * \param omit Index of molecule that should be omitted. Set to negative to include all.
    * \param NNd (IN- and OUTPUT) stores the nearest translational distance found so far.
    * \param NNs (IN- and OUTPUT) stores the nearest 6D distance found so far.
    */
  void searchVectorsForNearestNeighbors6D(Vec3 center,
                                          float W4, float X4, float Y4, float Z4,
                                          const Farray& V_XYZ_vec, const Farray& V_Q_vec, int omit, double& NNd, double& NNs);

  /**
   * Search a 3D grid for the nearest neighbor of a molecule in 3D (translation) and 6D (translation+rotation)
   *
   * \param center coordinates of the reference molecule
   * \param W4 molecule orientation (quaternion)
   * \param X4 molecule orientation (quaternion)
   * \param Y4 molecule orientation (quaternion)
   * \param Z4 molecule orientation (quaternion)
   * \param V_XYZ vector (length N_voxels) of vectors (length N_other*3) of other molecules' positions
   * \param V_Q vector (length N_voxels) of vectors (length N_other*4) of other molecules' quaternions
   * \param grid_Nx number of grid voxels in x dimension.
   * \param grid_Ny number of grid voxels in y dimension.
   * \param grid_Nz number of grid voxels in z dimension.
   * \param grid_origin origin (xyz pos of the center of the first voxel) of the grid.
   * \param grid_spacing spacing of the grid, same in all directions.
   * \param n_layers maximum number of layers of voxels to search for nearest neighbors
   * \param omit_in_central_vox index of molecule in the central voxel to omit.
   * \return std::pair<double, double> SQUARED distances in 3D and 6D.
   */
  std::pair<double, double> searchGridNearestNeighbors6D(
      Vec3 center,
      float W4, float X4, float Y4, float Z4,
      const std::vector<Farray>& V_XYZ, const std::vector<Farray>& V_Q,
      int grid_Nx, int grid_Ny, int grid_Nz,
      Vec3 grid_origin,
      double grid_spacing,
      int n_layers, int omit_in_central_vox);

  const double SIX_CORR_SPACING = 0.01;

  // Index of voxel (x, y, z) in a grid with size (Nx, Ny, Nz).
  inline int voxel_num(int x, int y, int z, int Nx, int Ny, int Nz) {
      return (x * Ny + y) * Nz + z;
  }

  /**
   * Linear interpolation with extrapolation for out-of-bounds values.
   * 
   * Given discrete function f(x) with *values* at positions starting from zero with given *spacing*,
   * return the linearly interpolated value at *x*.
   * If x is out of bounds, extrapolates linearly using the two first (or last) values.
   * 
   * @param  x       position to interpolate values at
   * @param  values  discrete function values
   * @param  spacing spacing of x positions
   */
  double interpolate(double x, const double* values, size_t array_len, double spacing);

  // Interpolate SIX_CORR at a given position
  double sixVolumeCorrFactor(double NNs);
}
#endif