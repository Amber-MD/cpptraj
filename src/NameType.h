#ifndef INC_NAMETYPE_H
#define INC_NAMETYPE_H
#include <string>
class NameType {
    static const unsigned int NameSize_ = 6;
  public:
    NameType();
    NameType(const NameType&);
    NameType(const char*);
    NameType(std::string const&);
    NameType& operator=(const NameType&);

    void ToBuffer(char*) const;
    bool Match(NameType const&) const;
    bool operator==(const NameType&) const;
    bool operator==(const char*) const;
    bool operator!=(const NameType&) const;
    bool operator!=(const char*) const;
    const char* operator*() const { return c_array_; }
    char operator[](int) const;
    std::string Truncated() const;
    void ReplaceAsterisk();
    bool operator<(NameType const& rhs) const {
      for (unsigned int i = 0; i != NameSize_; i++)
      {
        if      (c_array_[i] == '\0' && rhs.c_array_[i] == '\0') return false;
        else if (c_array_[i] == '\0' && rhs.c_array_[i] != '\0') return true;
        else if (c_array_[i] != '\0' && rhs.c_array_[i] == '\0') return false;
        else if (c_array_[i] < rhs.c_array_[i]) return true;
        else if (c_array_[i] > rhs.c_array_[i]) return false;
      }
      return false;
    }
    bool operator>(NameType const& rhs) const {
      for (unsigned int i = 0; i != NameSize_; i++)
      {
        if      (c_array_[i] == '\0' && rhs.c_array_[i] == '\0') return false;
        else if (c_array_[i] != '\0' && rhs.c_array_[i] == '\0') return true;
        else if (c_array_[i] == '\0' && rhs.c_array_[i] != '\0') return false;
        else if (c_array_[i] > rhs.c_array_[i]) return true;
        else if (c_array_[i] < rhs.c_array_[i]) return false;
      }
      return false;
    }
  private:
    char c_array_[NameSize_];

    void FormatName();
};
#endif
