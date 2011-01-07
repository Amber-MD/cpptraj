#ifndef INC_AMBERRESTARTNC_H
#define INC_AMBERRESTARTNC_H
#ifdef BINTRAJ
// AmberRestartNC
#include "TrajFile.h"

class AmberRestartNC : public TrajFile {
  int ncid;
  int atomDID;
  int ncatom;
  int coordVID;
  int hasVelocity;
  int velocityVID;
  double velocityScale;
  int cellAngleVID;
  int cellLengthVID;

  int spatialDID;
  int labelDID;
  int cell_spatialDID;
  int cell_angularDID;
  int spatialVID;
  int timeVID;
  double restartTime;
  int cell_spatialVID;
  int cell_angularVID;
  int TempVID;

  int SetupRead();
  int SetupWrite();
  int setupWriteForSet(int);

  public:

  AmberRestartNC();
  ~AmberRestartNC();

  int open();
  void close();
  int getFrame(int);
  int writeFrame(int);
  void Info();
};

#endif
#endif  
