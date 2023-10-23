#ifndef INC_STRUCTURE_BUILDATOM_H
#define INC_STRUCTURE_BUILDATOM_H
namespace Cpptraj {
namespace Structure {
/// Hold information for an atom used when building/modelling new coordinates.
class BuildAtom {
  public:
    BuildAtom();
  private:
    bool positionKnown_; ///< True if position is "known", i.e. can be used to build.
    
#endif
