#ifndef INC_ACTION_MINMAXDIST_H
#define INC_ACTION_MINMAXDIST_H
#include "Action.h"
#include "ImageOption.h"
#include "InteractionData.h"
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
        Entity(std::string const& na, int nu) : name_(na), num_(nu) {}
        AtomMask emask_;   ///< Selected atoms in entity.
        std::string name_; ///< residue/molecule name
        int num_;          ///< residue/molecule number
    };
    typedef std::vector<Entity> Earray;

    /// Used to track active sets and corresponding entities
    class ActiveSet {
      public:
        ActiveSet(DataSet* ds, Earray::const_iterator const& it1, Earray::const_iterator const& it2) :
          ds_(ds), it1_(it1), it2_(it2) {}
        DataSet* ds_;
        Earray::const_iterator it1_;
        Earray::const_iterator it2_;
    };
    typedef std::vector<ActiveSet> DSarray;

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    /// For debug, print entities to STDOUT
    void printEntities(Earray const&, AtomMask const&) const;
    /// Set up entity array for current mode
    int setupEntityArray(Earray&, AtomMask const&, Topology const&) const;
    /// Get minimum distance between atoms in masks
    double get_min_dist(AtomMask const&, AtomMask const&, Frame const&) const;
    /// Get minimum distance between atoms in mask
    double get_min_dist(AtomMask const&, Frame const&) const;
    /// Get maximum distance between atoms in masks
    double get_max_dist(AtomMask const&, AtomMask const&, Frame const&) const;
    /// Get maximum distance between atoms in mask
    double get_max_dist(AtomMask const&, Frame const&) const;

    AtomMask mask1_;
    AtomMask mask2_;
    ModeType mode_;
    DistType distType_;
    ImageOption imageOpt_; ///< Used to determine if imaging should be used.
    Earray entities1_; ///< Entities corresponding to mask1_
    Earray entities2_; ///< Entities corresponding to mask2_
    std::string dsname_; ///< Data set name
    DataFile* outfile_;  ///< Output data file
    DataSet* byAtomSet_; ///< By atom output data set
    DataSetList* masterDSL_; ///< Pointer to master data set list
    Cpptraj::InteractionData interactionSets_; ///< Hold interaction sets
    DSarray activeSets_;                       ///< Hold active interaction sets
    int debug_;
    int resOffset_;          ///< If BY_RES, do not look at pairs separated by less than resOffset_ in sequence
};
#endif
