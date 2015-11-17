#ifndef INC_READLINE_H
#define INC_READLINE_H
#include "CmdInput.h"
/// Wrapper around GNU readline library
class ReadLine {
  public:
    ReadLine();
    int GetInput();
    void AddHistory(const char*);
    bool YesNoPrompt(const char*);
    const char* c_str()            const { return input_.str();   }
    std::string const& operator*() const { return input_.Str();   }
    bool empty()                   const { return input_.Empty(); }
  private:
    CmdInput input_;
};
#endif
