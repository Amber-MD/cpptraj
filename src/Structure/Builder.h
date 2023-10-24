#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
  public:
    /// CONSTRUCTOR
    Builder();
    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int);
  private:
    typedef std::vector<bool> Barray;

    int debug_;
};
}
}
#endif
