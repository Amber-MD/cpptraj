#ifndef INC_ASSOCIATEDDATA_NOE_H
#define INC_ASSOCIATEDDATA_NOE_H
#include "AssociatedData.h"
class ArgList;
/// For Analysis_Statistics DISTANCE NOE
class AssociatedData_NOE : public AssociatedData {
  public:
    AssociatedData_NOE() : AssociatedData(NOE), l_bound_(0.0), u_bound_(0.0), rexp_(-1.0) {}
    AssociatedData_NOE(double l, double u, double r) :
      AssociatedData(NOE), l_bound_(l), u_bound_(u), rexp_(r) {}
    static const char* HelpText;
    int NOE_Args(ArgList&);
    double NOE_bound()  const { return l_bound_;  }
    double NOE_boundH() const { return u_bound_; }
    double NOE_rexp()   const { return rexp_;   }

    AssociatedData* Copy() const { return new AssociatedData_NOE(*this); }
    void Ainfo() const;
  private:
    double l_bound_; ///< Lower bound
    double u_bound_; ///< Upper bound
    double rexp_;    ///< Expected distance
};
#endif
