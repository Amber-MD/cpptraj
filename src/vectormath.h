#ifndef INC_VECTORMATH_H
#define INC_VECTORMATH_H
/*! \file vectormath.h
    \brief Collection of vector and matrix math routines.
 */
#define CROSS_PRODUCT(TX,TY,TZ,UX,UY,UZ,VX,VY,VZ) \
  TX = (UY*VZ) - (UZ*VY); \
  TY = (UZ*VX) - (UX*VZ); \
  TZ = (UX*VY) - (UY*VX)

void normalize(double a[3]);
double vector_norm(double V[3], double *);
void vector_sub(double V[3], double U[3], double W[3]);
void vector_sum(double V[3], double U[3], double W[3]);
double dot_product(double V[3], double U[3]);
double dot_product_angle(double V[3], double U[3]);
void cross_product(double V[3], double U[3], double W[3]);
double dot_product_sign(double *, double *, double *);

void matrix_transpose(double M[9], double U[9]);
void matrix_transpose_3x3(double M[9]);
double *matrix_transpose(double *M, int mrows, int ncols);
void matrix_times_vector(double U[3], double R[9], double V[3]);
void matrixT_times_vector(double U[3], double R[9], double V[3]);
void matrix_multiply_3x3(double M[9], double Row[9], double Col[9]);
double matrix_to_angle(double U[9]);
int axis_of_rotation(double V[3], double U[9], double theta);
void calcRotationMatrix(double T[9], double V[3], double theta);
void calcRotationMatrix(double T[9], double psiX, double psiY, double psiZ);

int diagEsort(double *mat, double *Emat, double *Evec[], double *Eigenvalue);

void printVector(const char *Name, double V[3]);
void printMatrix_3x3(const char *Title, double U[9]);
void printMatrix(const char *Title, double *U, int mrows, int ncols);
void printRotTransInfo(double U[9], double trans[6]);
#endif
