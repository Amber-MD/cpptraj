#ifndef INC_ATOMTOPTYPE_H
#define INC_ATOMTOPTYPE_H
/// Hold atom topology information which can be used when re-ordering
class AtomTopType {
  public:
    enum PdbType { ATOM = 0, HETATM };
    /// Type, index, resnum, icode, chain
    AtomTopType(PdbType t, int i, int r, char e, char c) :
      type_(t), index_(i), resnum_(r), icode_(e), chainid_(c) {}

    /// \return true if less than (Chain, resnum, icode, index)
    bool operator<(const AtomTopType& rhs) const {
      if (type_ == rhs.type_) {
        if (chainid_ == rhs.chainid_) {
          if (resnum_ == rhs.resnum_) {
            if (icode_ == rhs.icode_) {
              return (index_ < rhs.index_);
            } else {
              return (icode_ < rhs.icode_);
            }
          } else {
            return (resnum_ < rhs.resnum_);
          }
        } else {
          return (chainid_ < rhs.chainid_);
        }
      } else {
        return (type_ < rhs.type_);
      }
    }
    /// \return True if equal
    bool operator==(const AtomTopType& rhs) const {
      if (type_    != rhs.type_) return false;
      if (chainid_ != rhs.chainid_) return false;
      if (resnum_ != rhs.resnum_) return false;
      if (icode_ != rhs.icode_) return false;
      if (chainid_ != rhs.chainid_) return false;
      return true;
    }
    /// \return True if not equals
    bool operator!=(const AtomTopType& rhs) const {
      if (type_ != rhs.type_) return true;
      if (chainid_ != rhs.chainid_) return true;
      if (resnum_ != rhs.resnum_) return true;
      if (icode_ != rhs.icode_) return true;
      if (chainid_ != rhs.chainid_) return true;
      return false;
    }
    /// \return Atom index
    int AtomIdx() const { return index_; }
  private:
    PdbType type_; ///< The PDB atom type (ATOM, HETATM)
    int index_;    ///< The atom index
    int resnum_;   ///< The residue number
    char icode_;   ///< The residue insertion code.
    char chainid_; ///< The chain
};
#endif
