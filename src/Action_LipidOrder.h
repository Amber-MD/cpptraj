#ifndef INC_ACTION_LIPIDORDER_H
#define INC_ACTION_LIPIDORDER_H
#include <map>
#include "Action.h"
/// <Enter description of Action_LipidOrder here>
class Action_LipidOrder : public Action {
  public:
    Action_LipidOrder() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_LipidOrder(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

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
};

/// Hold information for a lipid carbon site.
class Action_LipidOrder::CarbonSite {
  public:
    CarbonSite() : c_idx_(-1), h1_idx_(-1), h2_idx_(-1), h3_idx_(-1) {}
    CarbonSite(int c) : c_idx_(c), h1_idx_(-1), h2_idx_(-1), h3_idx_(-1) {}
    int Cidx()  const { return c_idx_;  }
    int H1idx() const { return h1_idx_; }
    int H2idx() const { return h2_idx_; }
    int H3idx() const { return h3_idx_; }
    void AddHindex(int h) {
      if      (h1_idx_ == -1) h1_idx_ = h;
      else if (h2_idx_ == -1) h2_idx_ = h;
      else                    h3_idx_ = h; // TODO check for 4+ h?
    }
  private:
    int c_idx_;  ///< Carbon atom index
    int h1_idx_; ///< First hydrogen index
    int h2_idx_; ///< Second hydrogen index if present.
    int h3_idx_; ///< Third hydrogen index if present.
};

/// Hold information for all carbons with the same name.
class Action_LipidOrder::CarbonList {
  public:
    CarbonList() {}
    CarbonList(NameType const& n) : resname_(n) {}
    const char* resName() const { return *resname_; }
    unsigned int Nsites() const { return sites_.size(); }
    CarbonSite& AddSite(int cidx) {
      sites_.push_back( CarbonSite(cidx) );
      return sites_.back();
    }
    Carray::const_iterator begin() const { return sites_.begin(); }
    Carray::const_iterator end()   const { return sites_.end();   }
  private:
    NameType resname_;
    Carray sites_;
};
#endif
