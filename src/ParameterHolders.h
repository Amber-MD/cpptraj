#ifndef INC_PARAMETERHOLDERS_H
#define INC_PARAMETERHOLDERS_H
#include <vector>
#include <utility>
#include "NameType.h"
#include "ParameterTypes.h"
#include "AtomType.h"
namespace ParameterHolders {
  enum RetType { ADDED = 0, SAME, UPDATED, ERR };
} /* END namespace ParameterHolders */

/// Used to hold two or more atom type names. TODO rename to TypeNameHolder?
class AtomTypeHolder {
  public:
    typedef std::vector<NameType> Narray;
    typedef Narray::const_iterator const_iterator;
    AtomTypeHolder() {}
    /// CONSTRUCTOR - Take single atom type name
    AtomTypeHolder(NameType const& nameIn) : types_(1, nameIn) {}
    /// CONSTRUCTOR - Take array of atom type names
    AtomTypeHolder(Narray const& namesIn) : types_(namesIn) {}
    /// CONSTRUCTOR - Reserve space for given number of type names
    AtomTypeHolder(int size) { types_.clear(); types_.reserve(size); }
    /// CONSTRUCTOR - Set wildcard name and reserve space for given number of type names.
    AtomTypeHolder(int size, NameType const& wc) : wildcard_(wc) { types_.clear(); types_.reserve(size); }
    /// Add atom type name.
    void AddName(NameType const& n) { types_.push_back( n ); }
    /// \return Iterator to beginning of type name array.
    const_iterator begin() const { return types_.begin(); }
    /// \return Iterator to end of type name array.
    const_iterator end() const { return types_.end(); }
    /// \return number of types in holder
    unsigned int Size() const { return types_.size(); }
    /// \return size used in memory
    size_t DataSize() const { return (types_.size() * NameType::DataSize()) + NameType::DataSize(); } 
    /// \return Type name at index
    NameType const& operator[](int idx) const { return types_[idx]; }
    /// \return true if either direction is a match, taking into account wildcard.
    bool operator==(AtomTypeHolder const& rhs) const {
      // Sanity check
      if (types_.size() != rhs.types_.size()) return false;
      // Forwards direction
      bool match = true;
      for (unsigned int idx = 0; idx != types_.size(); idx++)
        if (types_[idx] != rhs.types_[idx] && types_[idx] != wildcard_) {
          match = false;
          break;
        }
      if (match) return true;
      // Reverse direction
      match = true;
      unsigned int idx2 = types_.size() - 1;
      for (unsigned int idx = 0; idx != types_.size(); idx++, idx2--)
        if (types_[idx] != rhs.types_[idx2] && types_[idx] != wildcard_) {
          match = false;
          break;
        }
      return match;
    }
    /// Will sort by type names in ascending order.
    bool operator<(AtomTypeHolder const& rhs) const {
      if (types_.size() != rhs.types_.size()) {
        return (types_.size() < rhs.types_.size());
      }
      for (unsigned int idx = 0; idx != types_.size(); idx++)
        if (types_[idx] < rhs.types_[idx])
          return true;
        else if (types_[idx] != rhs.types_[idx])
          return false;
      return false;
    }
    /// \return string containing atom type names.
    std::string TypeString() const {
      std::string tstr;
      for (Narray::const_iterator it = types_.begin(); it != types_.end(); ++it)
        tstr.append( " " + std::string( *(*it) ) );
      return tstr;
    }
  private:
    Narray types_;
    NameType wildcard_;
};

// -----------------------------------------------------------------------------
/// Used to associate atom type names with an object (parameter etc)
template <class T> class ParmHolder {
    // TODO may want to actually use a map one day for performance reasons.
    typedef std::pair<AtomTypeHolder,T> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    ParmHolder() {}
    void clear()              { bpmap_.clear(); }
    unsigned int size() const { return bpmap_.size(); }
    bool empty()        const { return bpmap_.empty(); }
    /// Add (or update if allowed) given parameter to holder.
    ParameterHolders::RetType AddParm(AtomTypeHolder const& types, T const& bp, bool allowUpdate) {
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.begin();
      for (; it != bpmap_.end(); ++it)
        if (it->first == types) break;
      if (it == bpmap_.end()) {
        // New parm
        bpmap_.push_back( Bpair(types, bp) );
      } else {
        if (bp < it->second || it->second < bp) {
          if (allowUpdate) {
            it->second = bp;
            return ParameterHolders::UPDATED;
          } else {
            return ParameterHolders::ERR;
          }
        } else
          return ParameterHolders::SAME;
      }
      return ParameterHolders::ADDED;
    }
    /// Constant iterator
    typedef typename Bmap::const_iterator const_iterator;
    /// \return constant iterator to beginning
    const_iterator begin() const { return bpmap_.begin(); }
    /// \return constant iterator to end.
    const_iterator end()   const { return bpmap_.end();   }
    /// Iterator
    typedef typename Bmap::iterator iterator;
    /// \return iterator to beginning
    iterator begin() { return bpmap_.begin(); }
    /// \return iterator to end
    iterator end()   { return bpmap_.end();   }
    /// \return Parameter matching given types, or empty parameter if not found.
    T FindParam(AtomTypeHolder const& types, bool& found) const { // TODO only use GetParam()?
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first == types) return it->second;
      found = false;
      return T();
    }
    /// \return iterator to parameter matching the given types.
    iterator GetParam(AtomTypeHolder const& types) {
      for (iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first == types) return it;
      return bpmap_.end();
    }
    /// \return const iterator to parameter matching the given types.
    const_iterator GetParam(AtomTypeHolder const& types) const {
      for (const_iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first == types) return it;
      return bpmap_.end();
    }
    /// \return size in memory in bytes
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

// -----------------------------------------------------------------------------
/// Specialized class for associating atom types with dihedral parameters.
/** NOTE: Instead of using a specialize template here I'm creating a new
  *       class because while I want AddParm() to accept DihedralParmType,
  *       I want FindParam to return an array of DihedralParmType, one for
  *       each unique multiplicity.
  */
class DihedralParmHolder {
    typedef std::pair<AtomTypeHolder,DihedralParmArray> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    DihedralParmHolder() {}
    void clear()              { bpmap_.clear();        }
    unsigned int size() const { return bpmap_.size();  }
    bool empty()        const { return bpmap_.empty(); }
    /** Add (or update) a single dihedral parameter for given atom types. */
    ParameterHolders::RetType
    AddParm(AtomTypeHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first == types)
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        bpmap_.push_back( Bpair(types, DihedralParmArray(1, dp)) );
      } else {
        // If we are here types match - check multiplicity.
        DihedralParmArray::iterator it1 = it0->second.begin();
        for (; it1 != it0->second.end(); ++it1)
        {
          if (it1->Pn() == dp.Pn())
            break;
        }
        if (it1 == it0->second.end()) {
          // Brand new multiplicity for this dihedral.
          it0->second.push_back( dp );
        } else {
          if (dp < *it1 || *it1 < dp) {
            if (allowUpdate) {
              *it1 = dp;
              return ParameterHolders::UPDATED;
            } else {
              return ParameterHolders::ERR;
            }
          } else
            return ParameterHolders::SAME;
        }
      }
      return ParameterHolders::ADDED;
    }

    /** This version takes an array of dihedral parameters. */
    ParameterHolders::RetType
    AddParm(AtomTypeHolder const& types, DihedralParmArray const& dpa, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first == types)
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        bpmap_.push_back( Bpair(types, dpa) );
      } else {
        if (!allowUpdate) return ParameterHolders::ERR;
        // Check if sizes are the same.
        bool update = false;
        if (it0->second.size() != dpa.size())
          update = true;
        else {
          // Sizes are the same. See if parameters are the same.
          for (unsigned int i = 0; i != it0->second.size(); i++) {
            if (it0->second[i] < dpa[i] || dpa[i] < it0->second[i]) {
              update = true;
              break;
            }
          }
        }
        if (update) {
          it0->second = dpa;
          return ParameterHolders::UPDATED;
        } else
          return ParameterHolders::SAME;
      }
      return ParameterHolders::ADDED;
    }

    typedef typename Bmap::const_iterator const_iterator;
    const_iterator begin() const { return bpmap_.begin(); }
    const_iterator end()   const { return bpmap_.end();   }
    /// \return Array of dihedral parameters matching given atom types.
    DihedralParmArray FindParam(AtomTypeHolder const& types, bool& found) const {
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first == types) return it->second;
      found = false;
      return DihedralParmArray();
    }
    /// \return size in memory in bytes
    size_t DataSize() const {
      if (bpmap_.empty()) return 0;
      const_iterator elt0 = begin();
      // Assume all AtomTypeHolders are the same size
      return (bpmap_.size() * elt0->first.DataSize()) +
             (bpmap_.size() * sizeof(DihedralParmArray)) +
             sizeof(Bmap);
    }
  private:
    Bmap bpmap_;
};
#endif
