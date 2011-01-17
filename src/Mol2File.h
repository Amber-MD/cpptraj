#ifndef INC_MOL2FILE_H
#define INC_MOL2FILE_H
// Mol2File
#include "TrajFile.h"

class Mol2File : public TrajFile {
    int mol2atom;
  public :

    Mol2File();
    ~Mol2File();

    int open();
    void close();
    int getFrame(int);
    int SetupRead();
    int SetupWrite();
    int writeFrame(int);
    void Info();
};
#endif
