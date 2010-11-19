#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
// PDBfile
#include <cstdio>
#include "TrajFile.h"

class PDBfile: public TrajFile {
    char buffer[256];
  public:

    PDBfile();
    ~PDBfile();

    int open();
    void close();
    int getFrame(int);
    int SetupRead();
    int SetupWrite();
    int writeFrame(int);
    void Info();   
};
#endif
