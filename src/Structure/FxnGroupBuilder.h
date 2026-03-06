#ifndef INC_STRUCTURE_FXNGROUPBUILDER_H
#define INC_STRUCTURE_FXNGROUPBUILDER_H
#include "FunctionalGroup.h"
// Forward declares
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
class Sugar;
/// Class to help prepare recognized functional groups
// TODO Determine if we need GetGroup() 
class FxnGroupBuilder {
  public:
    class FxnResType;
    /// CONSTRUCTOR - take debug level
    FxnGroupBuilder(int);
    /// Determine if sugar has functional groups that need to be split off into separate residues
    int CheckForFunctionalGroups(Sugar&, Topology&, Frame&);
    /// Make any necessary modifications to found functional groups in Topology
    int ModifyFxnGroups(Topology&, bool&) const;
    /// Print LEAP-related warnings for found functional groups
    void LeapFxnGroupWarning(Topology const&) const;

    typedef std::vector<FxnResType>::const_iterator const_iterator;
    const_iterator begin() const { return foundGroups_.begin(); }
    const_iterator end() const { return foundGroups_.end(); }
    unsigned int NgroupsFound() const { return foundGroups_.size(); }
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<FunctionalGroup> FGarray;
    typedef std::vector<FxnResType> FRarray;
    typedef std::vector<FunctionalGroup::Type> FTarray;
    /// \return identity of the group bonded to given atom
    FunctionalGroup::Type IdFunctionalGroup_Silent(Iarray&, int, int, int, Topology const&) const;
    /// \return identity of the group bonded to given atom, print to stdout
    FunctionalGroup::Type IdFunctionalGroup(Iarray&, int, int, int, Topology const&) const;

    /// Add recognized functional groups to array
    int AddGroups();
    /// \return Functional group if one is recognized
    FGarray::const_iterator GetGroup(Iarray&, Iarray const&, int, int, Topology const&) const;

    /// Determine if sugar anomeric carbon has a terminal functional group
    int CheckIfSugarIsTerminal(Sugar&, Topology&, Frame&);

    FGarray functionalGroups_; ///< Functional groups to search for (FunctionalGroupType).
    FRarray foundGroups_;      ///< Groups that have been found
    FTarray resTypes_;         ///< Groups that have been found in current call of CheckForFunctionalGroups()
    int debug_;                ///< Debug level

};
// -----------------------------------------------------------------------------
/// Used to keep track of found functional groups
class FxnGroupBuilder::FxnResType {
  public:
    FxnResType() : ires_(-1), ftype_(FunctionalGroup::UNRECOGNIZED_GROUP) {}
    FxnResType(int ir, FunctionalGroup::Type ft) : ires_(ir), ftype_(ft) {}
    int Ires() const { return ires_; }
    FunctionalGroup::Type Ftype() const { return ftype_; }
  private:
    int ires_; ///< Residue # of the functional group
    FunctionalGroup::Type ftype_; ///< Type of functional group
};
}
}
#endif
