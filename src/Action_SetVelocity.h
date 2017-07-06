#ifndef INC_ACTION_SETVELOCITY_H
#define INC_ACTION_SETVELOCITY_H
#include "Action.h"
#include "Random.h"
/// Calculate the temperature of parts of a system.
class Action_SetVelocity : public Action {
  public:
    Action_SetVelocity();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_SetVelocity(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    int AddBonds(BondArray const&, Topology const&, CharMask const&);
    int Rattle2(Frame&) const;

    /// Hold atom indices and bond force constant for each constrained bond
    class Cbond {
      public:
        Cbond() : rk_(0.0), at1_(-1), at2_(-1) {}
        Cbond(int a1, int a2, double rk) : rk_(rk), at1_(a1), at2_(a2) {}
        bool operator<(Cbond const& rhs) const {
          if (at1_ == rhs.at1_)
            return (at2_ < rhs.at2_);
          else
            return (at1_ < rhs.at1_);
        }
        double rk_;
        int at1_;
        int at2_;
    };
    typedef std::vector<Cbond> Carray;
    typedef std::vector<double> Darray;

    AtomMask Mask_;
    Darray SD_;      ///< Hold sqrt(kB*(1/mass)) for each atom
    Carray Bonds_;   ///< Hold constrained bonds
    double tempi_;
    double EPS_;
    int ntc_;        ///< Constraint type: 1=None, 2=Hydrogen, 3=All
    Random_Number RN_;
    CoordinateInfo cInfo_;
    Frame newFrame_;
};
#endif
