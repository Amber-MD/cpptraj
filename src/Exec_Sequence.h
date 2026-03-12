#ifndef INC_EXEC_SEQUENCE_H
#define INC_EXEC_SEQUENCE_H
#include "Exec.h"
namespace Cpptraj {
namespace Structure { 
class Creator;
}
}
/// Create a molecule from a sequence of units 
class Exec_Sequence : public Exec {
  public:
    Exec_Sequence();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Sequence(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<std::string> Sarray;
    typedef std::vector<DataSet_Coords*> Uarray;
    enum ModeType { UNSPECIFIED = 0, NEW, OLD };

    int get_units(Uarray&, Sarray const&, Cpptraj::Structure::Creator const&) const;

    int old_generate_sequence(DataSet_Coords*, Uarray const&) const;

    int generate_sequence(DataSet_Coords*, std::string const& title, Uarray const&, Cpptraj::Structure::Creator const&) const;

    int debug_;
    ModeType mode_;
};
#endif
