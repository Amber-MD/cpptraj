#ifndef INC_PARM_GROMACS_H
#define INC_PARM_GROMACS_H
#include "ParmIO.h"
class Parm_Gromacs : public ParmIO {
  public :
    Parm_Gromacs() : debug_(0), numOpen_(0) { }
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Gromacs(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1;   }
    void SetDebug(int i)                               { debug_ = i; }
    int processWriteArgs(ArgList&) { return 0; }
  private:
    int ReadGmxFile(std::string const&);

    typedef std::vector<Atom> AtomArray;
    typedef std::vector<AtomArray> MolArray;
    MolArray gmx_molecules_;
    typedef std::vector<std::string> Sarray;
    Sarray gmx_molnames_;

    int debug_;
    int numOpen_;
};
#endif
