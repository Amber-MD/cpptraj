#ifndef INC_CONSTRAINTS_H
#define INC_CONSTRAINTS_H
#include "Topology.h"
#include "ArgList.h"
/// Used to apply constraints
class Constraints {
  public:
    Constraints();
    enum ShakeType { OFF = 0, BONDS_TO_H, ALL_BONDS };
    static const char* constraintArgs;
    int InitConstraints(ArgList&);
    static const char* rattleArgs;
    int InitRattle(ArgList&);
    int SetupConstraints(AtomMask const&, Topology const&);
    int Rattle2(Frame&) const;
    ShakeType Type()       const { return shakeType_;          }
    int DegreesOfFreedom() const { return degrees_of_freedom_; }
    double DT()            const { return dt_;                 }
    double Epsilon()       const { return epsilon_;            }
    const char* shakeString() const;
  private:
    int AddBonds(BondArray const&, Topology const&, CharMask const&);

    /// Hold atom indices and bond eq. length for each constrained bond
    class Cbond {
      public:
        Cbond() : req_(0.0), at1_(-1), at2_(-1) {}
        Cbond(int a1, int a2, double req) : req_(req), at1_(a1), at2_(a2) {}
        bool operator<(Cbond const& rhs) const {
          if (at1_ == rhs.at1_)
            return (at2_ < rhs.at2_);
          else
            return (at1_ < rhs.at1_);
        }
        double req_;
        int at1_;
        int at2_;
    };
    typedef std::vector<Cbond> Carray;

    Carray Bonds_;           ///< Hold constrained bonds
    double dt_;              ///< Time step
    double epsilon_;         ///< epsilon
    double EPS_;             ///< Hold epsilon / dt
    ShakeType shakeType_;    ///< What bonds constraints being applied to.
    int degrees_of_freedom_; ///< Unconstrained deg. of freedom.
};
#endif
