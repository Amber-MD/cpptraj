#ifndef INC_AMBERNETCDF_H
#define INC_AMBERNETCDF_H
#ifdef BINTRAJ
/* amber_netcdf.h
 * Dan Roe 10-2008
 */
#include "TrajFile.h"

class AmberNetcdf : public TrajFile {
  int ncid;
  int frameDID;
  int ncframe;
  int atomDID;
  int ncatom; 
  float *Coord;
  int coordVID;
  int cellAngleVID;
  int cellLengthVID;

  int spatialDID;
  int labelDID;
  int cell_spatialDID;
  int cell_angularDID;
  int spatialVID;
  int timeVID;
  int cell_spatialVID;
  int cell_angularVID;
  int TempVID;

  int SetupRead();
  int SetupWrite();
  int WriteArgs(ArgList *A);

  public:
  
  AmberNetcdf();
  ~AmberNetcdf();

  int open();
  void close();
  int getFrame(int);
  int writeFrame(int);
  void Info();
};
#endif
#endif
