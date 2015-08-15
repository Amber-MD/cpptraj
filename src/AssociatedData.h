#ifndef INC_ASSOCIATEDDATA_H
#define INC_ASSOCIATEDDATA_H
#include "ArgList.h"
class AssociatedData {
  public:
    enum AssociatedType { NOE = 0 };
    AssociatedData(AssociatedType t) : type_(t) {}
    AssociatedType Type() { return type_; }
    virtual AssociatedData* Copy() const = 0;
  private:
    AssociatedType type_;
};

/// For Analysis_Statistics DISTANCE NOE
class AssociatedData_NOE : public AssociatedData {
  public:
    AssociatedData_NOE() : AssociatedData(NOE), l_bound_(0.0), u_bound_(0.0), rexp_(-1.0) {}
    static const char* HelpText;
    int NOE_Args(ArgList&);
    //void SetNOE(double l, double u, double r) { l_bound_=l; u_bound_=u; rexp_=r;}
    double NOE_bound()  const { return l_bound_;  }
    double NOE_boundH() const { return u_bound_; }
    double NOE_rexp()   const { return rexp_;   }

    AssociatedData* Copy() const { return new AssociatedData_NOE(*this); }
  private:
    double l_bound_; ///< Lower bound
    double u_bound_; ///< Upper bound
    double rexp_;    ///< Expected distance
};
#endif
