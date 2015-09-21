#ifndef INC_PARM_PDB_H
#define INC_PARM_PDB_H
#include "ParmIO.h"
class Parm_PDB : public ParmIO {
  public :
    Parm_PDB() : readAsPQR_(false), readBox_(false), readConect_(true) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_PDB(); }
    static void ReadHelp();
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&);
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&) { return 1; }
    void SetDebug(int) {}
    int processWriteArgs(ArgList&) { return 0; }
    bool NeedsBondSearch() const { return true; }
  private:
    bool readAsPQR_;  ///< If true get charge and radius from occ/b factor cols
    bool readBox_;    ///< If true try to read CRYST1 record as box info.
    bool readConect_; ///< If true read bonds from CONECT records.
};
#endif
