#ifndef INC_NETCDFROUTINES_H
#define INC_NETCDFROUTINES_H
/*! \file NetcdfRoutines.h
    \brief Collection of useful subroutines for any code using NetCDF format.
 */
// Defines
#ifdef BINTRAJ
#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_LENGTHS "cell_lengths"
#define NCCELL_ANGULAR "cell_angular"
#define NCCELL_ANGLES "cell_angles"
#define NCCOORDS "coordinates"
#define NCVELO "velocities"
#define NCTEMPERATURE "temp0"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5
#define NCREMD_DIMENSION "remd_dimension"
#define NCREMD_GROUPNUM "remd_groupnum"
#define NCREMD_DIMTYPE "remd_dimtype"
#define NCREMD_INDICES "remd_indices"

// Routines
void NetcdfDebug(int);
int checkNCerr(int, const char *, ...);
int GetDimInfo(int, const char *, int *);
char *GetAttrText(int, int, const char *);
#endif
char *GetNetcdfConventions(char *);
void FloatToDouble(double *, float *, int);
void DoubleToFloat(float *, double *, int);
#endif
