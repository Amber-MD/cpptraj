#ifndef INC_PARAMETERHOLDERS_H
#define INC_PARAMETERHOLDERS_H
#include <vector>
#include <utility>
#include "NameType.h"
#include "ParameterTypes.h"
//#inc lude "CpptrajStdio.h" // DEBUG
/// Used to hold two or more atom type names.
class AtomTypeHolder {
  public:
    typedef std::vector<NameType> Narray;
    typedef Narray::const_iterator const_iterator;
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
    const_iterator begin() const { return types_.begin(); }
    const_iterator end() const { return types_.end(); }
    /// \return number of types in holder
    size_t Size() const { return types_.size(); }
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
    /// \return size used in memory
    size_t DataSize() const { return types_.size() * NameType::DataSize(); } 
/*
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
    }*/
  private:
    Narray types_;
};
// -----------------------------------------------------------------------------
/// Used to associate atom type names with an object (parameter etc)
template <class T> class ParmHolder {
    // TODO may want to actually use a map one day for performance reasons.
    typedef std::pair<AtomTypeHolder,T> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    ParmHolder() {}
    void clear() { bpmap_.clear(); }
    unsigned int size() const { return bpmap_.size(); }
/*
    static inline void PrintTypes(AtomTypeHolder const& types) {
      for (AtomTypeHolder::const_iterator it = types.begin(); it != types.end(); ++it)
        mprintf(" %s", *(*it));
      mprintf("\n");
    }
*/
    int AddParm(AtomTypeHolder const& types, T const& bp, bool allowUpdate) {
//       if (types.Size() != 2) {
//    mprinterr("Internal Error: ParmHolder::AddParm(): # types is not 2 (%zu)\n",
//              types.Size());
//    return -1;
//  }
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.begin();
      for (; it != bpmap_.end(); ++it)
        if (it->first == types) break;
      if (it == bpmap_.end()) {
        // New parm
        bpmap_.push_back( Bpair(types, bp) );
        //mprintf("\tAdded new bond params for ");
        //PrintTypes(types);
      } else {
        if (allowUpdate) {
          //mprintf("\tUpdating bond parameters for ");
          //PrintTypes(types);
          it->second = bp;
        } else {
          //mprinterr("Error: Update of bond params not allowed for ");
          //PrintTypes( types);
          return 1;
        }
      }
      return 0;
    }

    typedef typename Bmap::const_iterator const_iterator;
    const_iterator begin() const { return bpmap_.begin(); }
    const_iterator end()   const { return bpmap_.end();   }

    size_t DataSize() const {
      if (bpmap_.empty()) return 0;
      const_iterator elt0 = begin();
      // Assume all AtomTypeHolders are the same size
      return (bpmap_.size() * elt0->first.DataSize()) +
             (bpmap_.size() * sizeof(T)) +
             sizeof(Bmap);
    }
  private:
    Bmap bpmap_;
};
#endif
