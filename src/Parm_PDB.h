#ifndef INC_PARM_PDB_H
#define INC_PARM_PDB_H
#include "ParmIO.h"
class Parm_PDB : public ParmIO {
  public :
    Parm_PDB() : ConectMode_(UNSPECIFIED), LinkMode_(UNSPECIFIED), readAsPQR_(false), readBox_(false) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_PDB(); }
    static void ReadHelp();
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&);
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&) { return 1; }
    int processWriteArgs(ArgList&) { return 0; }
  private:
    enum ReadType { UNSPECIFIED = 0, READ, SKIP };
    ReadType ConectMode_; ///< Specify how to handle CONECT records.
    ReadType LinkMode_;   ///< Specify how to handle LINK records.
    bool readAsPQR_;      ///< If true get charge and radius from occ/b factor cols
    bool readBox_;        ///< If true try to read CRYST1 record as box info.
};
#endif
