#ifndef GIST_CUDA_CALC_CUH
#define GIST_CUDA_CALC_CUH
#if defined(__HIP_PLATFORM_HCC__)
#include <hip/hip_runtime.h>
#endif

#define HUGE_C 1e30f
#define BLOCKSIZE 16
#define SLOW_BLOCKSIZE 512

/**
 * Coordinates class.
 *
 * Holds three different values for the coordinates, to achieve
 * a faster access on the GPU (only a single read from the memory
 * instead of three different reads, also the memory is better
 * aligned).
 * @author Johannes Kraml
 */
class Coordinates {
public:
  float x;
  float y;
  float z;

  // Empty Constructor
  __host__ __device__
  Coordinates(): x(0), y(0), z(0) {}

  /**
   * Constructor using an array of three values.
   * @param array: The array that holds the x, y and z coordinates
   */
  __host__ __device__
  Coordinates(const double *array) {
    this->x = array[0];
    this->y = array[1];
    this->z = array[2];
  }

  /**
   * Copy Constructor.
   * @param other: The other object.
   */
  __host__ __device__
  Coordinates(const Coordinates &other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
  }

  /**
   * The assignment operator for an object of this class.
   * @param other: The other Coordinate object.
   * @return This object.
   */
  __host__ __device__
  Coordinates &operator=(const Coordinates &other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    return *this;
  }
  
};

/**
 * An implementation holding the A and B values for a
 * Lennard-Jones van der Waals interaction energy calculation.
 * @author Johannes Kraml
 */
class ParamsLJ{
public:
  float A;
  float B;

  // Empty constructor
  __host__ __device__
  ParamsLJ(): A(0), B(0) {}

  /**
   * Constructor from an array of values.
   * @param arr: The array that holds the A and B values.
   */
  __host__ __device__
  ParamsLJ(float *arr) {
    this->A = arr[0];
    this->B = arr[1];
  }

  /**
   * Constructor using two different floats.
   * @param A: The A value in the Lennard-Jones equation.
   * @param B: The B value in the Lennard-Jones equation.
   */
  __host__ __device__
  ParamsLJ(float A, float B) {
    this->A = A;
    this->B = B;
  }

  /**
   * Copy constructor for this class.
   * @param other: The other object of this class
   */
  __host__ __device__
  ParamsLJ(const ParamsLJ &other) {
    this->A = other.A;
    this->B = other.B;
  }

  /**
   * The assignment operator of this class.
   * @param other: The other object of this type.
   * @return this object.
   */
  __host__ __device__
  ParamsLJ &operator=(const ParamsLJ &other) {
    this->A = other.A;
    this->B = other.B;
    return *this;
  }
};

/**
 * The implementation of the Unit Cell on the GPU.
 * @author Johannes Kraml
 */
class UnitCell {
public:
  float array[9];
  
  // The empty constuctor.
  __host__ __device__
  UnitCell() {}

  /**
   * Constructor using an array.
   * @param arr: The array from which the values should be
   *             taken.
   */
  __host__ __device__
  UnitCell(float *arr) {
    this->array[0] = arr[0];
    this->array[1] = arr[1];
    this->array[2] = arr[2];
    this->array[3] = arr[3];
    this->array[4] = arr[4];
    this->array[5] = arr[5];
    this->array[6] = arr[6];
    this->array[7] = arr[7];
    this->array[8] = arr[8];
  }

  /**
   * Copy constructor.
   * @param other: The other object.
   */
  __host__ __device__
  UnitCell(const UnitCell &other) {
    this->array[0] = other.array[0];
    this->array[1] = other.array[1];
    this->array[2] = other.array[2];
    this->array[3] = other.array[3];
    this->array[4] = other.array[4];
    this->array[5] = other.array[5];
    this->array[6] = other.array[6];
    this->array[7] = other.array[7];
    this->array[8] = other.array[8];
  }

  /**
   * Assignment operator.
   * @param other: The other object.
   * @return This object.
   */
  __host__ __device__
  UnitCell &operator=(const UnitCell &other){
    this->array[0] = other.array[0];
    this->array[1] = other.array[1];
    this->array[2] = other.array[2];
    this->array[3] = other.array[3];
    this->array[4] = other.array[4];
    this->array[5] = other.array[5];
    this->array[6] = other.array[6];
    this->array[7] = other.array[7];
    this->array[8] = other.array[8];
    return *this;
  }

  /**
   * Access operator.
   * @param idx: The index which should be accessed.
   *
   * @exceptsafe Not safe.
   */
  __host__ __device__
  float operator[](int idx) {
    if (idx >= 0 && idx < 9) {
      return this->array[idx];
    }
    return 1;
  }

};

/**
 * Class holding the box information, i.e., which kind of box
 * and the box dimensions.
 */
class BoxInfo {
public:
  float array[9]; ///< box dimensions
  int boxinfo;
  
  // Empty constructor
  __host__ __device__
  BoxInfo(): boxinfo(0) {}

  /**
   * Constructor using an array and the boxinfo.
   * @param arr: The array holding the values for the box dimensions.
   * @param boxinfo: Which kind of box this is.
   */
  __host__ __device__
  BoxInfo(float *arr, int boxinfo) {
    this->array[0] = arr[0];
    this->array[1] = arr[1];
    this->array[2] = arr[2];
    this->array[3] = arr[3];
    this->array[4] = arr[4];
    this->array[5] = arr[5];
    this->array[6] = arr[6];
    this->array[7] = arr[7];
    this->array[8] = arr[8];
    this->boxinfo = boxinfo;
  }

  /**
   * Copy constuctor.
   * @param other: The other object.
   */
  __host__ __device__
  BoxInfo(const BoxInfo &other) {
    this->array[0] = other.array[0];
    this->array[1] = other.array[1];
    this->array[2] = other.array[2];
    this->array[3] = other.array[3];
    this->array[4] = other.array[4];
    this->array[5] = other.array[5];
    this->array[6] = other.array[6];
    this->array[7] = other.array[7];
    this->array[8] = other.array[8];
    this->boxinfo = other.boxinfo;
  }

  /**
   * Assignement operator.
   * @param other: The other object.
   * @return This object.
   */
  __host__ __device__
  BoxInfo &operator=(const BoxInfo &other){
    this->array[0] = other.array[0];
    this->array[1] = other.array[1];
    this->array[2] = other.array[2];
    this->array[3] = other.array[3];
    this->array[4] = other.array[4];
    this->array[5] = other.array[5];
    this->array[6] = other.array[6];
    this->array[7] = other.array[7];
    this->array[8] = other.array[8];
    this->boxinfo = other.boxinfo;
    return *this;
  }

  /**
   * Access operator.
   * @param idx: The index at which to access the array.
   * @return The value stored at that point in the array or 0 if none is found.
   * @exceptsafe Is not exception safe.
   */
  __host__ __device__
  float operator[](int idx) {
    if (idx < 9 && idx >= 0) {
      return this->array[idx];
    }
    return 0;
  }

};

/**
 * Class to store different Atom properties, like charge, atomType, solvent and the
 * molecule this atom belongs to.
 * @author Johannes Kraml
 */
class AtomProperties {
public:
  float charge;
  int atomType;
  bool solvent;
  int molecule;

  // Empty constructor
  __host__ __device__
  AtomProperties() {}

  /**
   * Constructor.
   * @param charge: The charge of the atom.
   * @param atomType: The atom type for access of the lennard-jones parameters.
   * @param solvent: Is this atom a solvent atom?
   * @param molecule: The molecule this atom belongs to.
   */
  __host__ __device__
  AtomProperties(float charge, int atomType, bool solvent, int molecule) {
    this->charge = charge;
    this->atomType = atomType;
    this->solvent = solvent;
    this->molecule = molecule;
  }

  /**
   * Copy constructor.
   * @param other: The other object.
   */
  __host__ __device__
  AtomProperties(const AtomProperties &other) {
    this->charge = other.charge;
    this->atomType = other.atomType;
    this->solvent = other.solvent;
    this->molecule = other.molecule;
  }

  /**
   * Assignement operator.
   * @param other: The other object.
   * @return This object.
   */
  __host__  __device__
  AtomProperties &operator=(const AtomProperties &other) {
    this->charge = other.charge;
    this->atomType = other.atomType;
    this->solvent = other.solvent;
    this->molecule = other.molecule;
    return *this;
  }
};


// Device functions
__device__ float dist2_imageOrtho(float *, float *, BoxInfo);
__device__ void scalarProd(float* , BoxInfo , float *);
__device__ float dist2_imageNonOrtho(float *, float *, BoxInfo, UnitCell);
__device__ float calcIfDistIsSmaller(float *, float *, int , int , int , UnitCell, float );
__device__ float dist2_imageNonOrthoRecip(float * , float * , UnitCell);
__device__ float dist2_noImage(float *, float *);
__device__ float calcTotalEnergy(float , float , float , float , float);
__device__ float calcVdWEnergy(float , float , float );
__device__ float calcElectrostaticEnergy(float, float, float);
__device__ int getLJIndex(int , int , int *, int );
__device__ ParamsLJ getLJParam(int , int , int *, int , ParamsLJ *);
__device__ bool isOnGrid(float *, float *, float *);
__device__ float calcDist(float *, float *, BoxInfo, UnitCell);

// Global functions
__global__ void cudaCalcEnergy    (Coordinates *, int *, int, ParamsLJ *, AtomProperties *, BoxInfo, UnitCell, int, float *, float *, float *, float *, int, float, int *, int *);
__global__ void cudaCalcEnergySlow(Coordinates *, int *, int, ParamsLJ *, AtomProperties *, BoxInfo, UnitCell, int, float *, float *,	float *, float *, int, float, int *, int *);
#endif
