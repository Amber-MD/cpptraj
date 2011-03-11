#ifndef INC_MOL2FILE_H
#define INC_MOL2FILE_H
// Mol2File
#include "TrajFile.h"

class Mol2File : public TrajFile {
    int mol2atom;
    int mol2bonds;
    AmberParm::NAME *Types;
    int writeMode; // 0=single mol2, 1=single mol2 with multi MOLECULE, 2=multi mol2s
  public :

    Mol2File();
    ~Mol2File();

    int open();
    void close();
    int getFrame(int);
    int SetupRead();
    int WriteArgs(ArgList *);
    int SetupWrite();
    int writeFrame(int);
    void Info();
};
#endif
