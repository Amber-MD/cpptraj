#ifndef INC_NAMETYPE_H
#define INC_NAMETYPE_H
#include <cstddef> // size_t
#include <string>
class NameType {
  public:
    NameType();
    NameType(const NameType&);
    NameType(const char*);
    NameType(std::string const&);
    //NameType( char*);
    NameType& operator=(const NameType&);

//    void AssignNoFormat(const char*);
    void ToBuffer(char*);
    bool Match(NameType &);

    //NameType& operator=(const char*);
    bool operator==(const NameType&);
    bool operator==(const char*);
    bool operator!=(const NameType&);
    bool operator!=(const char*);
    //char *Assign(char*);

    const char * operator*() const;

    char operator[](int) const;
    void ReplaceAsterisk();

  private:
    const size_t NameSize_;
    char c_array_[6];

//    void FormatName(bool);
    void FormatName();
};  
#endif
