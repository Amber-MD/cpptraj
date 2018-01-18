#ifndef INC_PARAMETERHOLDERS_H
#define INC_PARAMETERHOLDERS_H
#include <map>
#include <vector>
#include "NameType.h"
#include "ParameterTypes.h"
class AtomTypeHolder {
  public:
    typedef std::vector<NameType> Narray;
    AtomTypeHolder() {}
    AtomTypeHolder(Narray const& namesIn) : types_(namesIn) {}
    AtomTypeHolder(int size) { types_.clear(); types_.reserve(size); }
    void AddName(NameType const& n) { types_.push_back( n ); }
/*    int SetTypes(Narray const& namesIn) {
      if (namesIn.size() != types_.size()) return 1;
      for (unsigned int idx = 0; idx != namesIn.size(); idx++)
        types_[idx] = namesIn[idx];
      return 0;
    }*/
    /// \return number of types in holder
    unsigned int Size() const { return types_.size(); }
    /// \return Type name at index
    NameType const& operator[](int idx) const { return types_[idx]; }
    /// \return true if either direction is a match.
    bool operator==(AtomTypeHolder const& rhs) const {
      // Sanity check
      if (types_.size() != rhs.types_.size()) return false;
      // Forwards direction
      bool match = true;
      for (unsigned int idx = 0; idx != types_.size(); idx++)
        if (types_[idx] != rhs.types_[idx]) {
          match = false;
          break;
        }
      if (match) return true;
      // Reverse direction
      match = true;
      unsigned int idx2 = types_.size() - 1;
      for (unsigned int idx = 0; idx != types_.size(); idx++, idx2--)
        if (types_[idx] != rhs.types_[idx2]) {
          match = false;
          break;
        }
      return match;
    }
    /// The lowest type name of either end is used for sorting.
    bool operator<(AtomTypeHolder const& rhs) const {
      if (types_.size() != rhs.types_.size()) {
        return (types_.size() < rhs.types_.size());
      }
      int lastidx = types_.size() - 1;
      int idx0;
      if (types_[0] < types_[lastidx])
        idx0 = 0;
      else
        idx0 = lastidx;
      int idx1;
      if (rhs.types_[0] < rhs.types_[lastidx])
        idx1 = 0;
      else
        idx1 = lastidx;
      return (types_[idx0] < rhs.types_[idx1]);
    }
  private:
    Narray types_;
};
// -----------------------------------------------------------------------------
class BondParmHolder {
    typedef std::map<AtomTypeHolder,BondParmType> Bmap;
    typedef std::pair<AtomTypeHolder,BondParmType> Bpair;
  public:
    BondParmHolder() {}
    int AddBondParm(AtomTypeHolder const&, BondParmType const&, bool);
  private:
    Bmap bpmap_;
};
#endif
