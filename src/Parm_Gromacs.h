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

    class gmx_atom;
    typedef std::vector<gmx_atom> AtomArray;
    typedef std::vector<AtomArray> MolArray;
    MolArray gmx_molecules_;
    typedef std::vector<std::string> Sarray;
    Sarray gmx_molnames_;
    Sarray mols_; ///< Molecules present
    std::vector<int> nums_; ///< How much of each mol is present.

    int debug_;
    int numOpen_;
    std::string title_;
};
// ----- CLASSES ---------------------------------------------------------------
class Parm_Gromacs::gmx_atom {
  public:
    gmx_atom() : charge_(0.0), mass_(0.0), rnum_(0) {}
    gmx_atom(NameType const& an, NameType const& at, NameType const& rn,
             double c, double m, int r) :
             aname_(an), atype_(at), rname_(rn), charge_(c), mass_(m), rnum_(r) {}
    
    NameType aname_;
    NameType atype_;
    NameType rname_;
    double charge_;
    double mass_;
    int rnum_;
};
#endif
