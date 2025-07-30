#ifndef INC_PARM_CIF_H
#define INC_PARM_CIF_H
#include "ParmIO.h"
// Forward declarations
class CIFfile;
class Parm_CIF : public ParmIO {
  public :
    Parm_CIF();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CIF(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&) { return 1;   }
    int processWriteArgs(ArgList&) { return 0; }
  private:
    enum EntryType { ANAME=0, RNAME, X, Y, Z, RNUM, CHAINID, NENTRY };
    static const char* Entries_[];
    /// Read _atom_site data block
    int read_atom_site(CIFfile&, FileName const&, Topology&) const;
    /// Read _chem_comp_atom
    int read_chem_comp_atom(CIFfile&, FileName const&, Topology&) const;

    bool read_struct_conn_;
};
#endif
