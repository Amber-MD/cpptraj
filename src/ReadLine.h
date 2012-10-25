#ifndef INC_READLINE_H
#define INC_READLINE_H
#include <string>
/// Wrapper around GNU readline library
class ReadLine {
  public:
    ReadLine() {}
   void GetInput();
   const char* c_str()            { return input_.c_str(); }
   std::string const& operator*() { return input_;         }
  private:
    std::string input_;
};
#endif
