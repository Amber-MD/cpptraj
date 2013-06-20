#ifndef INC_READLINE_H
#define INC_READLINE_H
#include <string>
/// Wrapper around GNU readline library
class ReadLine {
  public:
    ReadLine() {}
    int GetInput();
    bool YesNoPrompt(const char*);
    const char* c_str()            const { return input_.c_str(); }
    std::string const& operator*() const { return input_;         }
    bool empty()                   const { return input_.empty(); }
  private:
    std::string input_;
};
#endif
