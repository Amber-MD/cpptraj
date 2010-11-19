#ifndef INC_AMBERRESTART_H
#define INC_AMBERRESTART_H
// AmberRestart.h
#include "TrajFile.h"

class AmberRestart : public TrajFile {
  int frameSize;        // Size of 1 coord frame in bytes, inc box & velo if present
  char *frameBuffer;    // Used to read in restart coord
  int numBoxCoords;     // Number of box coords (3 or 6)
  int hasVelocity;      // This restart file has velocity info
  int restartAtoms;     // Number of atoms in restart file, read in
  double restartTime;   // Time in restart file, read in
  double restartTemp;   // (Optional) replica temperature, read in.

  int SetupRead();
  int SetupWrite();

  public:

  //Frame *V;             // Hold velocities

  AmberRestart();
  ~AmberRestart();
  int open();
  void close();
  int getFrame(int);
  int writeFrame(int);
  void Info();
};
#endif
