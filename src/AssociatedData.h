#ifndef INC_ASSOCIATEDDATA_H
#define INC_ASSOCIATEDDATA_H
class AssociatedData {
  public:
    enum AssociatedType { NOE = 0 };
    AssociatedData(AssociatedType t) : type_(t) {}
    AssociatedType Type() { return type_; }
  private:
    AssociatedType type_;
};

/// For Analysis_Statistics DISTANCE NOE
class AssociatedData_NOE : public AssociatedData {
  public:
    AssociatedData_NOE() : AssociatedData(NOE), bound_(0.0), boundh_(0.0), rexp_(-1.0) {}
    void SetNOE(double b, double bh, double r) { bound_=b; boundh_=bh; rexp_=r;}
    double NOE_bound()  const { return bound_;  }
    double NOE_boundH() const { return boundh_; }
    double NOE_rexp()   const { return rexp_;   }
  private:
    double bound_; ///< Lower bound
    double boundh_; ///< Upper bound
    double rexp_; ///< Expected distance
};
#endif
