#ifndef INC_TYPENAMEHOLDER_H
#define INC_TYPENAMEHOLDER_H
#include <cstddef> // size_t
#include "NameType.h"
/// Used to hold one or more atom type names.
class TypeNameHolder {
  public:
    typedef std::vector<NameType> Narray;
    typedef Narray::const_iterator const_iterator;
    TypeNameHolder() {}
    /// CONSTRUCTOR - Take single atom type name
    TypeNameHolder(NameType const& nameIn) : types_(1, nameIn) {}
    /// CONSTRUCTOR - Take array of atom type names
    TypeNameHolder(Narray const& namesIn) : types_(namesIn) {}
    /// CONSTRUCTOR - Reserve space for given number of type names
    TypeNameHolder(int size) { types_.clear(); types_.reserve(size); }
    /// CONSTRUCTOR - Set wildcard name and reserve space for given number of type names.
    TypeNameHolder(int size, NameType const& wc) : wildcard_(wc) { types_.clear(); types_.reserve(size); }
    /// Add atom type name.
    void AddName(NameType const& n) { types_.push_back( n ); }
    /// \return Iterator to beginning of type name array.
    const_iterator begin() const { return types_.begin(); }
    /// \return Iterator to end of type name array.
    const_iterator end() const { return types_.end(); }
    /// \return number of types in holder
    size_t Size() const { return types_.size(); }
    /// \return Type name at index
    NameType const& operator[](int idx) const { return types_[idx]; }
    /// \return true if either direction is a match, taking into account wildcard.
    bool operator==(TypeNameHolder const& rhs) const {
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
    bool operator<(TypeNameHolder const& rhs) const {
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
    /// \return size in bytes
    size_t DataSize() const { return (types_.size()*NameType::DataSize()) + NameType::DataSize(); }
  private:
    Narray types_;
    NameType wildcard_;
};
#endif
