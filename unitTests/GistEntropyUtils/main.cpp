// Unit test for NameType class
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>

#include <string>
#include <algorithm>
#include "GistEntropyUtils.h"
#include "Vec3.h"

#define PI_SQR 9.869604401089358

using namespace GistEntropyUtils;

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

template<typename T>
bool close(T a, T b)
{
  T atol = 1e-7;
  T rtol = 1e-7;
  return fabs(a - b) < (atol + fabs(b) * rtol);
}

#define assertClose(a, b, msg) \
  do { \
    if (!close(a, b)) { \
      std::cerr << "Not close: " << a << " != " << b << std::endl; \
      return Err(msg); \
      } \
    } while(0)

#define assertEqual(a, b, msg) \
  do { \
    if (a != b) { \
      std::cerr << "Not equal: " << a << " != " << b << std::endl; \
      return Err(msg); \
      } \
    } while(0)

/* Helper function to test the translation part of searchVectorsForNearestNeighbors6D */
double searchTrans(Vec3 center,
                   float W4, float X4, float Y4, float Z4,
                   const Farray& V_XYZ, const Farray& V_Q, int omit)
{
  double NNd = GIST_HUGE;
  double NNs = GIST_HUGE;
  searchVectorsForNearestNeighbors6D(center, W4, X4, Y4, Z4, V_XYZ, V_Q, omit, NNd, NNs); 
  return NNd;
}

/* Helper function to test the 6D part of searchVectorsForNearestNeighbors6D */
double searchSix(Vec3 center,
                 float W4, float X4, float Y4, float Z4,
                 const Farray& V_XYZ, const Farray& V_Q, int omit)
{
  double NNd = GIST_HUGE;
  double NNs = GIST_HUGE;
  searchVectorsForNearestNeighbors6D(center, W4, X4, Y4, Z4, V_XYZ, V_Q, omit, NNd, NNs);
  return NNs;
}

int main() {
  const double vals[] = { 0, 1, 3 };
  double spacing = 2;
  assertClose(interpolate(1, vals, 3, spacing), 0.5, "Wrong interpolation at x=1.");
  assertClose(interpolate(0, vals, 3, spacing), 0., "Wrong interpolation at x=0.");
  assertClose(interpolate(3, vals, 3, spacing), 2., "Wrong interpolation at x=3.");
  // different spacing
  assertClose(interpolate(0.3, vals, 3, 0.2), 2., "Wrong interpolation at x=0.3 with non-integer spacing.");
  // extrapolate linearly to the left
  assertClose(interpolate(-4, vals, 3, spacing), -2., "Wrong interpolation at x=-4.");
  // extrapolate linearly to the right
  assertClose(interpolate(8, vals, 3, spacing), 7., "Wrong interpolation at x=8.");
  
  assertClose(sixVolumeCorrFactor(0), 1., "Wrong six corr factor at x=0.");
  assertClose(sixVolumeCorrFactor(4), 1.785045060669792516,  "Wrong six corr factor at x=4.");

  // NN search tests
  Farray testCrd;
  testCrd.push_back(1); testCrd.push_back(0); testCrd.push_back(0);
  testCrd.push_back(2); testCrd.push_back(0); testCrd.push_back(0);
  testCrd.push_back(3); testCrd.push_back(0); testCrd.push_back(0);
  Farray testQ;
  testQ.push_back(1); testQ.push_back(0); testQ.push_back(0); testQ.push_back(0);
  testQ.push_back(0); testQ.push_back(1); testQ.push_back(0); testQ.push_back(0);
  testQ.push_back(0); testQ.push_back(1); testQ.push_back(0); testQ.push_back(0);

  // test that finds trans nearest neighbor
  assertClose(searchTrans(Vec3(0.9, 0, 0), 1, 0, 0, 0, testCrd, testQ, -1), 0.01, "Wrong translational distance");
  // test that finds orient nearest neighbor
  assertClose(searchSix(Vec3(1, 0, 0), 0.99875026, 0., 0.049979169, 0, testCrd, testQ, -1), 0.01, "Wrong orientational distance");
  // shift test point by one. Still the same nearest neighbor because of orientation
  assertClose(searchSix(Vec3(2, 0, 0), 0.99875026, 0., 0.049979169, 0, testCrd, testQ, -1), 1.01, "Wrong orientational distance");
  // rotate 180 deg around x axis. Now the second quat is the closest.
  assertClose(searchSix(Vec3(2, 0, 0), 0.049979169, 0.99875026, 0, 0, testCrd, testQ, -1), 0.01, "Wrong orientational distance");
  // Now, omit the second molecule.
  assertClose(searchSix(Vec3(2, 0, 0), 0.049979169, 0.99875026, 0, 0, testCrd, testQ, 1), 1.01, "Wrong orientational distance");


  assertEqual(voxel_num(1, 0, 0, 3, 3, 2), 6, "Wrong voxel_num for (1, 0, 0)");
  assertEqual(voxel_num(2, 1, 0, 3, 3, 2), 14, "Wrong voxel_num for (1, 0, 0)");
  std::vector<Farray> grid_crd(27);
  std::vector<Farray> grid_Q(27);
  int grid_Nx = 3;
  int grid_Ny = 3;
  int grid_Nz = 3;
  Vec3 grid_origin(0, 0, 0);
  double grid_spacing = 0.5;
  assertClose(searchGridNearestNeighbors6D(
    Vec3(1, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 2, -1).first,
    GIST_HUGE, "Wrong trans. dist. with empty grid.");

  // after adding a point in the same voxel, it should be found.
  grid_crd[9].push_back(0.6); grid_crd[9].push_back(0); grid_crd[9].push_back(0);
  grid_Q[9].push_back(1); grid_Q[9].push_back(0); grid_Q[9].push_back(0); grid_Q[9].push_back(0);
  std::pair<double, double>nbrs;
  nbrs = searchGridNearestNeighbors6D(Vec3(1, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 1, -1);
  assertClose(nbrs.first, 0.16, "Wrong trans. dist. with empty grid.");
  assertClose(nbrs.second, PI_SQR + 0.16, "Wrong trans. dist. with empty grid.");
  
  nbrs = searchGridNearestNeighbors6D(Vec3(1, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 0, -1);
  assertClose(nbrs.first, GIST_HUGE, "Does not omit outer shells in NN search.");

  nbrs = searchGridNearestNeighbors6D(Vec3(-1, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 10, -1);
  assertClose(nbrs.first, 1.6*1.6, "Does not work for out of bounds center.");

  // add another point one voxel higher, with different rotation;
  // Depending on rotation, this might be the NN. But might not be found if n_layers is < 2;
  grid_crd[18].push_back(1.1); grid_crd[18].push_back(0); grid_crd[18].push_back(0);
  grid_Q[18].push_back(0); grid_Q[18].push_back(1); grid_Q[18].push_back(0); grid_Q[18].push_back(0);

  nbrs = searchGridNearestNeighbors6D(Vec3(0.1, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 2, -1);
  assertClose(nbrs.first, 0.5 * 0.5, "Does not find correct trans NN.");
  assertClose(nbrs.second, 1.0 * 1.0, "Does not find correct six NN with n_layers=2.");

  nbrs = searchGridNearestNeighbors6D(Vec3(0.1, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 1, -1);
  assertClose(nbrs.first, 0.5*0.5, "Does not find correct trans NN.");
  assertClose(nbrs.second, 0.5*0.5 + PI_SQR, "Does not find correct six NN with n_layers=1.");

  // test that a molecule can be omitted in the central voxel.
  nbrs = searchGridNearestNeighbors6D(Vec3(0.6, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 1, 0);
  assertClose(nbrs.first, 0.5*0.5, "Does not find correct trans NN.");

  // test that no molecule is always omitted in the central voxel.
  nbrs = searchGridNearestNeighbors6D(Vec3(0.6, 0, 0), 0, 1, 0, 0, grid_crd, grid_Q, grid_Nx, grid_Ny, grid_Nz, grid_origin, grid_spacing, 1, -1);
  assertClose(nbrs.first, 0.0, "Does not find correct trans NN.");

  // Test quaternion distance with slightly non-normalized quaternions
  assertClose(quaternion_angle(
    0.453887561488233, 0.440500369544494, 0.615810865052675, 0.469811115705563,
    0.233334999318213, 0.253147448327817, 0.580026938731094, 0.738268174747172), 0.7979849037416394, "Wrong distance between normalized quaternions.");

  const double FAC = 1.0002;
  assertClose(quaternion_angle(
    0.453887561488233*FAC, 0.440500369544494*FAC, 0.615810865052675*FAC, 0.469811115705563*FAC,
    0.233334999318213, 0.253147448327817, 0.580026938731094, 0.738268174747172), 0.7979849037416394, "Wrong distance between non-normalized quaternions.");

  assertClose(quaternion_angle(
    0.453887561488233, 0.440500369544494, 0.615810865052675, 0.469811115705563,
    0.233334999318213*FAC, 0.253147448327817*FAC, 0.580026938731094*FAC, 0.738268174747172*FAC), 0.7979849037416394, "Wrong distance between non-normalized quaternions.");
  return 0;
}
