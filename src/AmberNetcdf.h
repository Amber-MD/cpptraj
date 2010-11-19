#ifndef INC_AMBERNETCDF_H
#define INC_AMBERNETCDF_H
#ifdef HASNETCDF
/* amber_netcdf.h
 * Dan Roe 10-2008
 */
#include "TrajFile.h"

class AmberNetcdf : public TrajFile {
  int ncid;
  int frameID;
  int ncframe;
  int atomID;
  int ncatom; 
  float *Coord;
  int coordID;
  int cellAngleID;
  int cellLengthID;

  int spatialID;
  int labelID;
  int cell_spatialID;
  int cell_angularID;
  int spatialVID;
  int timeVID;
  int cell_spatialVID;
  int cell_angularVID;
  int TempVarID;

  int checkNCerr(int, const char *, ...); 
  int GetDimInfo(const char *, int *);

  int SetupRead();
  int SetupWrite();

  public:
  
  AmberNetcdf();
  ~AmberNetcdf();
  void Debug();

  int open();
  void close();
  int getFrame(int);
  int writeFrame(int);
  void Info();
};
#endif
#endif
