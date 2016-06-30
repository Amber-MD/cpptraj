#ifndef INC_PARM_GROMACS_H
#define INC_PARM_GROMACS_H
#include "ParmIO.h"
#include "BufferedLine.h"
class Parm_Gromacs : public ParmIO {
  public :
    Parm_Gromacs() : numOpen_(0), directive_err_(0) { }
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Gromacs(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&) { return 1;   }
    int processWriteArgs(ArgList&) { return 0; }
  private:
    enum KeyType { G_UNKNOWN_KEY = 0, G_MOLECULE_TYPE, G_ATOMS, G_BONDS,
                   G_SYSTEM, G_MOLECULES, G_SETTLES, G_VIRTUAL_SITES3 };
    static const char* SEP;

    bool LineContainsKey(std::string const&, std::string const&) const;
    int LineContainsDirective(std::string const&, std::string const&, std::string&);
    int Defined(std::string const&) const;
    int AdvanceToElse(BufferedLine&) const;
    KeyType FindKey(std::string const&) const;
    int ReadAtomsSection(BufferedLine&);
    int ReadBondsSection(BufferedLine&);
    int ReadSettles(BufferedLine&);
    int ReadVsite3(BufferedLine&);
    int ReadMolsSection(BufferedLine&);
    int ReadGmxFile(FileName const&);

    class gmx_atom;
    class gmx_mol;
    typedef std::vector<gmx_atom> AtomArray;
    typedef std::vector<gmx_mol> MolArray;
    typedef std::vector<int> BondArray;
    MolArray gmx_molecules_; ///< Hold each molecule defined in topololgy
    typedef std::vector<std::string> Sarray;
    Sarray mols_;           ///< Which molecules are present
    Sarray defines_;        ///< Which preprocessor defines present.
    std::vector<int> nums_; ///< How much of each mol is present.

    int numOpen_;
    int directive_err_;
    std::string title_;
    std::string currentWorkDir_;
    FileName infileName_;
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

class Parm_Gromacs::gmx_mol {
  public:
    gmx_mol() {}
    gmx_mol(std::string const& s) : mname_(s) {}
    const char* Mname() { return mname_.c_str(); }

    AtomArray atoms_;
    BondArray bonds_;
    std::string mname_;
};
#endif
