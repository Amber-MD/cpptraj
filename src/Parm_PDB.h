#ifndef INC_PARM_PDB_H
#define INC_PARM_PDB_H
#include "ParmIO.h"
class Parm_PDB : public ParmIO {
  public :
    Parm_PDB() : readAsPQR_(false) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_PDB(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&);
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
  private:
    bool readAsPQR_; ///< If true get charge and radius from occ/b factor cols
};
#endif
