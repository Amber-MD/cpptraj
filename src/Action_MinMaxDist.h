#ifndef INC_ACTION_MINMAXDIST_H
#define INC_ACTION_MINMAXDIST_H
#include "Action.h"
#include "ImageOption.h"
/// Record the min/max distance between atoms/residues/molecules 
class Action_MinMaxDist : public Action {
  public:
    /// CONSTRUCTOR
    Action_MinMaxDist();
    /// ALLOCATOR
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MinMaxDist(); }
    /// HELP
    void Help() const;
  private:
    enum ModeType {BY_ATOM = 0, BY_RES, BY_MOL, NO_MODE};
    enum DistType { MIN_DIST = 0, MAX_DIST, BOTH_DIST, NO_DIST };

    static const char* modeStr_[];
    static const char* distTypeStr_[];

    /// Used to track residues/molecules
    class Entity {
      public:
        Entity() : num_(-1) {}
        bool IsSet() const { return (num_ > -1); }
        AtomMask emask_;   ///< Selected atoms in entity.
        std::string name_; ///< residue/molecule name
        int num_;          ///< residue/molecule number
    };

    typedef std::vector<Entity> Earray;

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    /// For debug, print entities to STDOUT
    void printEntities(Earray const&, AtomMask const&) const;
    /// Set up entity array for current mode
    int setupEntityArray(Earray&, AtomMask const&, Topology const&) const;

    AtomMask mask1_;
    AtomMask mask2_;
    ModeType mode_;
    DistType distType_;
    ImageOption imageOpt_; ///< Used to determine if imaging should be used.
    Earray entities1_; ///< Entities corresponding to mask1_
    Earray entities2_; ///< Entities corresponding to mask2_
};
#endif
