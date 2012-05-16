#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
#include "CpptrajFile.h"
#include "PDBtype.h"
class PDBfile : CpptrajFile, PDBtype {
  public:
    PDBfile();
    ~PDBfile();

    int OpenPDB(std::string const&);
    int OpenPDB(const char*);
    void WriteHET(int, double, double, double);
  private:
    int anum_;
};
#endif
