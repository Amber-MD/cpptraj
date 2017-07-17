#ifndef INC_ACTION_LIPIDORDER_H
#define INC_ACTION_LIPIDORDER_H
#include <map>
#include "Action.h"
/// <Enter description of Action_LipidOrder here>
class Action_LipidOrder : public Action {
  public:
    Action_LipidOrder();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_LipidOrder(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    static const unsigned int MAX_H_;

    enum AxisType { DX = 0, DY, DZ };

    /// Hold information for a single lipid carbon site
    class CarbonSite;
    typedef std::vector<CarbonSite> Carray;

    /// Hold information for all carbons with the same name.
    class CarbonList;

    /// Pair a carbon name with all carbon sites with that name.
    typedef std::pair<NameType, CarbonList> Cpair;

    /// Map carbon name to all carbon sites with that name.
    typedef std::map<NameType, CarbonList> Cmap;

    /// Create carbon site with hydrogens, then follow remaining carbon chain
    void FollowChain(int, Topology const&, std::vector<bool>, int);

    CharMask mask_;
    Cmap Carbons_;
    AxisType axis_;
};

/// Hold information for a lipid carbon site.
class Action_LipidOrder::CarbonSite {
  public:
    CarbonSite();
    /// Create site with carbon atom index
    CarbonSite(int);
    bool operator<(CarbonSite const& rhs) const { return (c_idx_ < rhs.c_idx_); }
    unsigned int NumH() const { return nH_; }
    int Cidx()          const { return c_idx_;  }
    int Hidx(int i)     const { return h_idx_[i]; }
    void AddHindex(int);
  private:
    int c_idx_;       ///< Carbon atom index.
    int h_idx_[3];    ///< Hydrogen atom indices.
    unsigned int nH_; ///< Number of hydrogen atom indices.
};

/// Hold information for all carbons with the same name.
class Action_LipidOrder::CarbonList {
  public:
    CarbonList();
    CarbonList(NameType const&);
    const char* resName() const { return *resname_; }
    unsigned int Nsites() const { return sites_.size(); }
    unsigned int Nvals()  const { return nvals_; }
    CarbonSite& AddSite(int);
    Carray::const_iterator begin() const { return sites_.begin(); }
    Carray::const_iterator end()   const { return sites_.end();   }
    Carray::iterator begin()             { return sites_.begin(); }
    Carray::iterator end()               { return sites_.end();   }
    void UpdateAngle(int i, double val) {
      sum_[i] += val;
      sum2_[i] += val * val;
    }
    void UpdateNvals() { nvals_++; }
    double Avg(int, double&) const;
  private:
    NameType resname_;
    Carray sites_;
    double sum_[3];
    double sum2_[3];
    unsigned int nvals_; ///< Number of times list has been updated, used to calc avg.
};
#endif
