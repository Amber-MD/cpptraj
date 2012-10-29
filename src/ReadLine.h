#ifndef INC_READLINE_H
#define INC_READLINE_H
#include <string>
/// Wrapper around GNU readline library
class ReadLine {
  public:
    ReadLine() {}
    int GetInput();
    const char* c_str()            { return input_.c_str(); }
    std::string const& operator*() { return input_;         }
    bool empty()                   { return input_.empty(); }
  private:
    std::string input_;
};
#endif
