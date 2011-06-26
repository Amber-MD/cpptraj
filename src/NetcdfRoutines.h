#ifndef INC_NETCDFROUTINES_H
#define INC_NETCDFROUTINES_H
/*
 * NetcdfRoutines
 * Collection of useful subroutines for any code using NetCDF format.
 */
// Defines
#ifdef BINTRAJ
#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_ANGULAR "cell_angular"
#define NCCOORDS "coordinates"
#define NCVELO "velocities"
#define NCTEMPERATURE "temp0"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5

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
