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
// TODO change function names, update help text
// TODO combine CheckForFunctionalGroups and CheckIfSugarIsTerminal
class FxnGroupBuilder {
 public:
   /// CONSTRUCTOR - take debug level
   FxnGroupBuilder(int);
    /// Determine if sugar has sulfates that need SO3 residue(s)
    int CheckForFunctionalGroups(Sugar&, Topology&, Frame&) const;
    /// Determine if sugar is terminal and need an ROH residue
    int CheckIfSugarIsTerminal(Sugar&, Topology&, Frame&) const;
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<FunctionalGroup> FGarray;
    /// \return identity of the group bonded to given atom
    FunctionalGroup::Type IdFunctionalGroup_Silent(Iarray&, int, int, int, Topology const&) const;
    /// \return identity of the group bonded to given atom, print to stdout
    FunctionalGroup::Type IdFunctionalGroup(Iarray&, int, int, int, Topology const&) const;

    /// Add recognized functional groups to array
    int AddGroups();
    /// \return Functional group if one is recognized
    FGarray::const_iterator GetGroup(Iarray&, Iarray const&, int, int, Topology const&) const;

    FGarray functionalGroups_; ///< Recognized functional groups (FunctionalGroupType).
    int debug_;                ///< Debug level

};
}
}
#endif
