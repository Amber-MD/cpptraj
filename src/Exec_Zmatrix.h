#ifndef INC_EXEC_ZMATRIX_H
#define INC_EXEC_ZMATRIX_H
#include "Exec.h"
/// Calculate a Zmatrix from a COORDS data set 
class Exec_Zmatrix : public Exec {
  public:
    Exec_Zmatrix() : Exec(COORDS), debug_(0) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Zmatrix(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int getZmatrix(DataSet_Coords*, int, int, std::string const&, DataFile*, CpptrajState&) const;
    int putZmatrix(DataSet_Coords*, Topology*, std::string const&, std::string const&, CpptrajState&) const;

    int debug_;
};
#endif
