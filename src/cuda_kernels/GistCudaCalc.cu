#include "GistCudaCalc.cuh"

#include <cstdio>

#define ELECTOAMBER_2 332.05221729


 /**
  * Calculate the squared distance in an orthorhombic box. See cpptraj implementation.
  * @param vec1: The first point of the distance calculation.
  * @param vec2: The seconf point of the distance calculation
  * @param box: The boxinfo of the object.
  * @return: The minimal distance in an orthorhombic box.
  */
__device__ 
float dist2_imageOrtho(float *vec1, float *vec2, BoxInfo box) {
  if (box[0] == 0 || box[1] == 0 || box[2] == 0) {
    return -1;
  }
  float x = abs(vec1[0] - vec2[0]);
  float y = abs(vec1[1] - vec2[1]);
  float z = abs(vec1[2] - vec2[2]);

  while (x > box[0]) {
    x -= box[0];
  }
  while (y > box[1]) {
    y -= box[1];
  }
  while (z > box[2]) {
    z -= box[2];
  }

  x = min(x, box[0] - x);
  y = min(y, box[1] - y);
  z = min(z, box[2] - z);

  return x * x + y * y + z * z;
}

/**
 * Calculate M * v.
 * @param vec: The vector v.
 * @param mat3x3: The matrix M.
 * @param ret: The values to be returned. If null, returns into vec.
 */
__device__
void scalarProd(float* vec, BoxInfo mat3x3, float *ret) {
  if (ret != NULL) {
    ret[0] = vec[0] * mat3x3[0] + vec[1] * mat3x3[1] + vec[2] * mat3x3[2];
    ret[1] = vec[0] * mat3x3[3] + vec[1] * mat3x3[4] + vec[2] * mat3x3[5];
    ret[2] = vec[0] * mat3x3[6] + vec[1] * mat3x3[7] + vec[2] * mat3x3[8];
  } else {
    float x = vec[0];
    float y = vec[1];
    float z = vec[2];
    vec[0] = x * mat3x3[0] + y * mat3x3[1] + z * mat3x3[2];
    vec[1] = x * mat3x3[3] + y * mat3x3[4] + z * mat3x3[5];
    vec[2] = x * mat3x3[6] + y * mat3x3[7] + z * mat3x3[8];
  }
}

/**
 * Find the result_O squared distance in a non-orthogonal box.
 * @param vec1: Position vector of point 1.
 * @param vec2: Position vector of point 2.
 * @param recip: The inverse of the ucell.
 * @param ucell: The box matrix.
 * @return: The minimal squared distance between two atoms, also considering the images.
 */
__device__
float dist2_imageNonOrtho(float *vec1, float *vec2, BoxInfo recip, UnitCell ucell) {
  float vecRecip1[3];
  float vecRecip2[3];
  scalarProd(vec1, recip, vecRecip1);
  scalarProd(vec2, recip, vecRecip2);
  float r_2 = dist2_imageNonOrthoRecip(vecRecip1, vecRecip2, ucell);
  
  return r_2;
}

/**
 * Calculate if the distance to an image is smaller than a given distance.
 * @param f: The first vector in the reciprocal space.
 * @param vec2Cartesian: The second vector in cartesian coordinates.
 * @param nx: Which neighbour in x direction?
 * @param ny: Which neighbour in y direction?
 * @param nz: Which neighbour in z direction?
 * @param ucell: The box matrix.
 * @param finalMin: The already calculated minimum.
 * @return: The new minimum, if it is smaller than finalMin, finalMin otherwise.
 */
__device__
float calcIfDistIsSmaller(float *f, float *vec2Cartesian, int nx, int ny, int nz, UnitCell ucell, float finalMin) {
  float fx = f[0] + nx;
  float fy = f[1] + ny;
  float fz = f[2] + nz;
  // Bring f back in Cartesian coordinates
  float x = fx * ucell[0] + fy * ucell[3] + fz * ucell[6];
  float y = fx * ucell[1] + fy * ucell[4] + fz * ucell[7];
  float z = fx * ucell[2] + fy * ucell[5] + fz * ucell[8];
  x -= vec2Cartesian[0];
  y -= vec2Cartesian[1];
  z -= vec2Cartesian[2];
  float min = x * x + y * y + z * z;
  if ( min < finalMin || finalMin < 0) {
    return min;
  }
  return finalMin;
}

/**
 * Calculate the distance in a non-orthorhombic box, if the vectors are already
 * multiplied by the inverse box matrix.
 * @param vec1: The first position vector.
 * @param vec2: The second position vector.
 * @param ucell: The box cell.
 * @return: The minimal distance between images.
 */
__device__
float dist2_imageNonOrthoRecip(float * vec1, float * vec2, UnitCell ucell) {
    
  // Bring the points back into the main unit cell
  float fx = vec1[0] - floor(vec1[0]);
  float fy = vec1[1] - floor(vec1[1]);
  float fz = vec1[2] - floor(vec1[2]);
  float f2x = vec2[0] - floor(vec2[0]);
  float f2y = vec2[1] - floor(vec2[1]);
  float f2z = vec2[2] - floor(vec2[2]);

  float vec[3] = {fx, fy, fz};

  // Bring f2 back in cartesian coordinates
  float xFactor = f2x * ucell[0] + f2y * ucell[3] + f2z * ucell[6];
  float yFactor = f2x * ucell[1] + f2y * ucell[4] + f2z * ucell[7];
  float zFactor = f2x * ucell[2] + f2y * ucell[5] + f2z * ucell[8];
  float vec2Real[3] = {xFactor, yFactor, zFactor};

  // Now the different cases, and always store the minimum
  // Define the finalMinimum as a negative value, since it can never
  // actually be negative this is fine.
  float finalMinimum = -1.0;

  // Run through all cells
  for (int ix = -1; ix <= 1; ++ix) {
    for (int iy = -1; iy <= 1; ++iy) {
      for (int iz = -1; iz <= 1; ++iz) {
        finalMinimum = calcIfDistIsSmaller(vec, vec2Real, ix, iy, iz, ucell, finalMinimum);
      }
    }
  }

  return finalMinimum;
}

/**
 * Calculates the distance between two vectors, without imaging.
 * @param vec1: The position vector of the first atom.
 * @param vec2: The position vector of the second atom.
 * @return: The squared distance between the two positions.
 */
__device__
float dist2_noImage(float *vec1, float *vec2) {
  float x = vec1[0] - vec2[0];
  float y = vec1[1] - vec2[1];
  float z = vec1[2] - vec2[2];

  return x*x + y*y + z*z;
}

/**
 * Caclulate the total energy between two atoms.
 * @param vec1: The position vector of atom 1.
 * @param vec2: The position vector of atom 2.
 * @param q1: The charge of atom 1.
 * @param q2: The charge of atom 2.
 * @param LJA: The lennard jones A parameter.
 * @param LJB: The lennard jones B parameter.
 * @param boxinfo: Which kind of box, 0 not periodic, 1 orthorhombic, 2 otherwise.
 * @param recip_o_box: Holds either the inverse of the cell matrix, if boxinfo is 2,
 *                        or the box dimensions, if boxinfo is 1, or is NULL, if boxinfo is 0.
 * @param ucell: Holds the cell matrix, if boxinfo is 2, NULL otherwise.
 * @return: The total interaction energy between the two atoms.
 */
__device__
float calcTotalEnergy(float q1, float q2, 
                            float LJA, float LJB, float r_2) {
#ifdef DEBUG_GIST_CUDA
  if (r_2 <= 0.000001 && r_2 >= -0.000001) {
    printf("(%8.3f, %8.3f, %8.3f) (%8.3f, %8.3f, %8.3f) %d\n", vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2], boxinfo);
  }
#endif
  return calcVdWEnergy(r_2, LJA, LJB) + calcElectrostaticEnergy(r_2, q1, q2);
}

/**
 * Calculate the distance between two different points.
 * @param vec1: The first vector to calculate the distance.
 * @param vec2: The second vector to calculate the distance.
 * @param recip_o_box: The boxinfo, either the box or the reciprocal.
 * @param ucell: The unitcell of a box.
 * @return: The squared distance between two points.
 */
__device__
float calcDist(float *vec1, float *vec2, BoxInfo recip_o_box,
                    UnitCell ucell) {
  float r_2 = 0;
  switch(recip_o_box.boxinfo) {
    case 0:
      r_2 = dist2_noImage(vec1, vec2);
      break;
    case 1:
      // Uses recip for box info as well;
      r_2 = dist2_imageOrtho(vec1, vec2, recip_o_box);
      break;
    case 2:
      r_2 = dist2_imageNonOrtho(vec1, vec2, recip_o_box, ucell);
      break;
    default:
      r_2 = 0;
  }
  return r_2;
}

/**
 * Calculate the Van der Waals energy between two atoms.
 * @param r_2: The squared distance between the two atoms.
 * @param LJA: The A part of the lennard jones potential.
 * @param LJB: The B part if the Lennard Jones potential.
 * @return: The Van der Waals energy.
 */
__device__
float calcVdWEnergy(float r_2, float LJA, float LJB) {
  float r_6 = r_2 * r_2 * r_2;
  float r_12 = r_6 * r_6;
  float LJ =  LJA / r_12 - LJB / r_6;
  return LJ;
}

/**
 * Calculate the electrostatic energy between two different atoms.
 * @param r_2: The square distance between the two atoms.
 * @param q1: The charge of atom 1.
 * @param q2: The charge of atom 2.
 * @return: The electrostatic energy between the two atoms.
 */
__device__
float calcElectrostaticEnergy(float r_2, float q1, float q2) {
  double charge = q1 * q2 * ELECTOAMBER_2;
  double r = sqrt(r_2);
  float ELE = charge / r;
  return ELE;
}

/**
 * Get the index into the lennard jones parameter array.
 * @param a1: The atom type index of atom 1.
 * @param a2: The atom type index of atom 2.
 * @param NBindex: The arrays holding the indices into the LJ array.
 * @param ntypes: The number of atom types.
 * @return: The index into the parameter arrays.
 */
__device__
int getLJIndex(int a1, int a2, int *NBindex, int ntypes) {
  return NBindex[a1 * ntypes + a2];
}

/**
 * Get the LJ parameters from a parameter array.
 * @param a1: The atom type index of atom 1.
 * @param a2: The atom type index of atom 2.
 * @param NBindex: The indices into the parameter array.
 * @param ntypes: The number of atom types.
 * @param paramsLJ: The parameter array.
 * @return: The LJ parameter belonging to the atom type pair a1, a2.
 */
__device__
ParamsLJ getLJParam(int a1, int a2, int *NBindex, int ntypes, ParamsLJ *paramsLJ) {
  int idx = getLJIndex(a1, a2, NBindex, ntypes);
  if (idx < 0) {
    return ParamsLJ();
  }
  return paramsLJ[idx];
}

/**
 * Checks if atom is on grid.
 * @param vec: Position vector to the atom.
 * @param min: The grid starting position.
 * @param max: The end of the grid.
 * @return: True if vector is on grid.
 */
__device__
bool isOnGrid(float *vec, float *min, float *max) {
  return ( ( (vec[0] >= min[0]) && (vec[0] <= max[0]) ) &&
           ( (vec[1] >= min[1]) && (vec[1] <= max[1]) ) &&
           ( (vec[2] >= min[2]) && (vec[2] <= max[2]) ) );
}

/**
 * Calculate the energy on the GPU.
 * @param coords: An array holding all the coordinates of all atoms.
 * @param NBindex: An array holding indices into the LJ parameter arrays.
 * @param ntypes: The number of atom types.
 * @param paramsLJA: The A LJ parameters.
 * @param paramsLJB: The B LJ parameters.
 * @param charges: The charges of the atoms.
 * @param boxinfo: Which kind of box, 0 not periodic, 1 orthorhombic, 2 otherwise.
 * @param recip_o_box: Holds either the inverse of the cell matrix, if boxinfo is 2,
 *                        or the box dimensions, if boxinfo is 1, or is NULL, if boxinfo is 0.
 * @param ucell: Holds the cell matrix, if boxinfo is 2, NULL otherwise.
 * @param maxAtoms: The number of atoms in the system.
 * @param a_types: The different atom types of the atoms.
 * @param solvent: True if atom is a solvent atom, false otherwise.
 * @param molecule: The number of the molecule this atom belong to.
 * @param result_ww: The result of the water - water interactions.
 * @param result_sw: The result of the solute - water interactions.
 * @param min: The minimum values of the grid.
 * @param max: The maximum values of the grid.
 */
__global__
void cudaCalcEnergy(Coordinates *coords, int *NBindex, int ntypes, ParamsLJ *parameterLJ, AtomProperties *atomProps, 
                          BoxInfo recip_o_box, UnitCell ucell, int maxAtoms, float *result_ww, float *result_sw, 
                          float *min, float *max, int headAtomType, float neighbourCut2, int *result_O, int *result_N) {
  

  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  int a2 = blockIdx.y * blockDim.y + threadIdx.y;
  if ( (a1 >= maxAtoms) || (a2 >= maxAtoms) || (a1 == a2)) {
    return;
  }

  AtomProperties atom1 = atomProps[a1];
  AtomProperties atom2 = atomProps[a2];

  // Do not calculate if the two values are the same or they belong to the same molecule.
  if ( (atom1.molecule != atom2.molecule)) {

    Coordinates t1 = coords[a1];
    Coordinates t2 = coords[a2];
    ParamsLJ lj = getLJParam(atom1.atomType, atom2.atomType, NBindex, ntypes, parameterLJ);

    float vec1[3] = {t1.x, t1.y, t1.z};
    float vec2[3] = {t2.x, t2.y, t2.z};
    float r_2 = calcDist(vec1, vec2, recip_o_box, ucell);
    float energy = calcTotalEnergy(atom1.charge, atom2.charge, lj.A, lj.B, r_2);
    
    if (atom2.solvent) {
      atomicAdd(&(result_ww[a1]), energy * 0.5f);
    } else {
      atomicAdd(&(result_sw[a1]), energy);
    }
  
  }
}

/**
 * Calculate the energy on the GPU. This implementation is somewhat slower,
 * but is able to calculate the order parameters as well as the rest.
 * @param coords: An array holding all the coordinates of all atoms.
 * @param NBindex: An array holding indices into the LJ parameter arrays.
 * @param ntypes: The number of atom types.
 * @param paramsLJA: The A LJ parameters.
 * @param paramsLJB: The B LJ parameters.
 * @param charges: The charges of the atoms.
 * @param boxinfo: Which kind of box, 0 not periodic, 1 orthorhombic, 2 otherwise.
 * @param recip_o_box: Holds either the inverse of the cell matrix, if boxinfo is 2,
 *                        or the box dimensions, if boxinfo is 1, or is NULL, if boxinfo is 0.
 * @param ucell: Holds the cell matrix, if boxinfo is 2, NULL otherwise.
 * @param maxAtoms: The number of atoms in the system.
 * @param a_types: The different atom types of the atoms.
 * @param solvent: True if atom is a solvent atom, false otherwise.
 * @param molecule: The number of the molecule this atom belong to.
 * @param result_ww: The result of the water - water interactions.
 * @param result_sw: The result of the solute - water interactions.
 * @param min: The minimum values of the grid.
 * @param max: The maximum values of the grid.
 */
__global__
void cudaCalcEnergySlow(Coordinates *coords, int *NBindex, int ntypes, ParamsLJ *parameterLJ, AtomProperties *atomProps, 
  BoxInfo recip_o_box, UnitCell ucell, int maxAtoms, float *result_ww, float *result_sw,
  float *min, float *max, int headAtomType, float neighbourCut2, int *result_O, int *result_N) {
  
  int a1 = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (a1 >= maxAtoms) {
    return;
  }
  
  AtomProperties atom1 = atomProps[a1];
  float distances[4] = {HUGE_C, HUGE_C, HUGE_C, HUGE_C};
  result_N[a1] = 0;
  result_O[4 * a1 + 3] = 0;
  result_O[4 * a1 + 2] = 0;
  result_O[4 * a1 + 1] = 0;
  result_O[4 * a1    ] = 0;
  float energy_ww = 0.0f;
  float energy_sw = 0.0f;
  for (int a2 = 0; a2 < maxAtoms; ++a2) {
    AtomProperties atom2 = atomProps[a2];
    // Do not calculate if the two values are the same or they belong to the same molecule.
    if ((a1 != a2) && (atom1.molecule != atom2.molecule)) {
      Coordinates t1 = coords[a1];
      Coordinates t2 = coords[a2];
      ParamsLJ lj = getLJParam(atom1.atomType, atom2.atomType, NBindex, ntypes, parameterLJ);
      float vec1[3] = {t1.x, t1.y, t1.z};
      float vec2[3] = {t2.x, t2.y, t2.z};
      float r_2 = calcDist(vec1, vec2, recip_o_box, ucell);
      float energy = calcTotalEnergy(atom1.charge, atom2.charge, lj.A, lj.B, r_2);
      if ((atom2.atomType == headAtomType) && atom2.solvent && atom1.solvent) {
        if (r_2 < distances[0]) {
          distances[3] = distances[2];
          distances[2] = distances[1];
          distances[1] = distances[0];
          distances[0] = r_2;
          result_O[4 * a1 + 3] = result_O[4 * a1 + 2];
          result_O[4 * a1 + 2] = result_O[4 * a1 + 1];
          result_O[4 * a1 + 1] = result_O[4 * a1    ];
          result_O[4 * a1    ] = a2;
        } else if (r_2 < distances[1]) {
          distances[3] = distances[2];
          distances[2] = distances[1];
          distances[1] = r_2;
          result_O[4 * a1 + 3] = result_O[4 * a1 + 2];
          result_O[4 * a1 + 2] = result_O[4 * a1 + 1];
          result_O[4 * a1 + 1] = a2;
        } else if (r_2 < distances[2]) {
          distances[3] = distances[2];
          distances[2] = r_2;
          result_O[4 * a1 + 3] = result_O[4 * a1 + 2];
          result_O[4 * a1 + 2] = a2;
        } else if (r_2 < distances[3]) {
          distances[3] = r_2;
          result_O[4 * a1 + 3] = a2;
        }
        if (r_2 < neighbourCut2) {
          result_N[a1] += 1;
        }
      }
      if (atom2.solvent) {
        energy_ww += energy * 0.5f;
      } else {
        energy_sw += energy;
      }
    }
  }
  result_ww[a1] = energy_ww;
  result_sw[a1] = energy_sw;
}
