#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
// PDBfile
#include "TrajFile.h"

class PDBfile: public TrajFile {
    char buffer[256];
    int pdbAtom;
  public:
    int writeMode; // 0=single pdb, 1=single pdb with MODEL keyword, 2=multiple pdbs

    PDBfile();
    ~PDBfile();

    int open();
    void close();
    int getFrame(int);
    int SetupRead();
    int WriteArgs(ArgList*);
    int SetupWrite();
    int writeFrame(int);
    void Info();   
};
#endif
