#ifndef INC_AMBERTRAJ_H
#define INC_AMBERTRAJ_H
#include "TrajFile.h"

/*
 * Class: AmberTraj
 * A type of TrajFile used for reading in formatted (text) Amber trajectories.
 * NOTE: Use size_t and off_t for anything involving file calcs to avoid
 * implicit type conversions?
 */
class AmberTraj: public TrajFile {
  int titleSize;        // Size of title in bytes 
  int frameSize;        // Size of 1 coord frame in bytes, inc box coords/REMD header if present
  int hasREMD;          // Size of REMD header if present
  //size_t frameSize;   // Use size_t because this will be used with fread?
  char *frameBuffer;    // Hold 1 coord frame (inc. box coords)
  int numBoxCoords;     // Number of box coords, 3 (Amber spec.) or 6 (future development)

  int SetupRead();
  int SetupWrite();
  int WriteArgs(ArgList *);

  public:

  AmberTraj();
  ~AmberTraj();
  int open();
  void close();
  int getFrame(int);
  int writeFrame(int);
  void Info();
};

#endif
