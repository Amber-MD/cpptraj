#ifndef INC_NAMETYPE_H
#define INC_NAMETYPE_H
#include <string>
/// Class for holding strings of small size, like atom names or type names.
class NameType {
    /** The size of the char array. Max name that can fit is ArraySize_ - 1 */
    static const unsigned int ArraySize_ = 8;
  public:
    NameType();
    NameType(const NameType&);
    NameType(const char*);
    NameType(std::string const&);
    NameType& operator=(const NameType&);
    /// Assign char buffer to this name
    void Assign(const char*);
    /// Copy this name to give char buffer.
    void ToBuffer(char*) const;
    /// \return true if given NameType (with optional wildcard chars) matches.
    bool Match(NameType const&) const;
    /// \return true if given NameType matches exactly.
    bool operator==(const NameType&) const;
    /// \return true if given NameType matches exactly.
    bool operator==(const char*) const;
    /// \return true if given NameType does not match exactly.
    bool operator!=(const NameType&) const;
    /// \return true if given NameType does not match exactly.
    bool operator!=(const char*) const;
    /// \return Const pointer to internal array.
    const char* operator*() const { return c_array_; }
    /// \return Character at specified position, or null if out of range
    char operator[](int) const;
    /// \return Name as a string.
    std::string Truncated() const;
    /// \return Name with minimal given width, padded with spaces if necessary.
    std::string Formatted(int) const;
    /// \return non-space length of name
    int len() const;
    /// \return true if name comes before given name alphabetically
    bool operator<(NameType const& rhs) const {
      for (unsigned int i = 0; i != ArraySize_; i++)
      {
        if      (c_array_[i] == '\0' && rhs.c_array_[i] == '\0') return false;
        else if (c_array_[i] == '\0' && rhs.c_array_[i] != '\0') return true;
        else if (c_array_[i] != '\0' && rhs.c_array_[i] == '\0') return false;
        else if (c_array_[i] < rhs.c_array_[i]) return true;
        else if (c_array_[i] > rhs.c_array_[i]) return false;
      }
      return false;
    }
    /// \return true if name comes after given name alphabetically
    bool operator>(NameType const& rhs) const {
      for (unsigned int i = 0; i != ArraySize_; i++)
      {
        if      (c_array_[i] == '\0' && rhs.c_array_[i] == '\0') return false;
        else if (c_array_[i] != '\0' && rhs.c_array_[i] == '\0') return true;
        else if (c_array_[i] == '\0' && rhs.c_array_[i] != '\0') return false;
        else if (c_array_[i] > rhs.c_array_[i]) return true;
        else if (c_array_[i] < rhs.c_array_[i]) return false;
      }
      return false;
    }
    /// \return size taken by this NameType in bytes.
    static size_t DataSize() { return ArraySize_ * sizeof(char); }
  private:
    char c_array_[ArraySize_];
};
#endif
