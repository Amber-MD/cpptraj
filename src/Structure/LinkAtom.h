#ifndef INC_STRUCTURE_LINKATOM_H
#define INC_STRUCTURE_LINKATOM_H
namespace Cpptraj {
namespace Structure {
/// Hold topology index and chain position for sugar link atom.
class LinkAtom {
  public:
    /// CONSTRUCTOR
    LinkAtom() : idx_(-1), position_(-1) {}
    /// CONSTRUCTOR - index, position (starting from 1 at the anomeric carbon)
    LinkAtom(int i, int p) : idx_(i), position_(p) {}
    /// \return Index in topology
    int Idx() const { return idx_; }
    /// \return Index in carbon chain (starting from 1 at the anomeric carbon)
    int Position() const { return position_; }
    /// First sort by position, then absolute index
    bool operator<(LinkAtom const& rhs) const {
      if (position_ == rhs.position_) {
        return (idx_ < rhs.idx_);
      } else {
        return position_ < rhs.position_;
      }
    }
  private:
    int idx_;      ///< Atom index in topology
    int position_; ///< Position in sugar carbon chain, starting from 1 at the anomeric carbon
};

}
}
#endif
