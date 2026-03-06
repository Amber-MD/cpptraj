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
    Exec_Sequence() : Exec(COORDS), debug_(0) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Sequence(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<std::string> Sarray;

    int generate_sequence(DataSet_Coords*, DataSetList const&,
                          Sarray const&, Cpptraj::Structure::Creator const&) const;

    int debug_;
};
#endif
