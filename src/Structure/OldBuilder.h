#ifndef INC_STRUCTURE_OLDBUILDER_H
#define INC_STRUCTURE_OLDBUILDER_H
#include <vector>
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
/// This is the old (<v7) builder used to attach different topology/frame combos using internal coordinates.
/** Retaining this class for backwards compatibility with modXNA TODO. */
class OldBuilder {
  public:
    /// CONSTRUCTOR
    OldBuilder();
    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int);
    /// Set debug
    void SetDebug(int d) { debug_ = d; }
  private:
    typedef std::vector<bool> Barray;

    int debug_;
};


}
}
#endif
