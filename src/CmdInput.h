#ifndef INC_CMDINPUT_H
#define INC_CMDINPUT_H
#include <string>
/// This class is used to format command input lines.
class CmdInput {
  public:
    CmdInput() {}
    /// Add and format given input; \return 1 if more input needed.
    int AddInput(const char*);
    std::string const& Str() const { return input_;         }
    const char* str()        const { return input_.c_str(); }
    void Clear()                   { input_.clear();        }
    bool Empty()             const { return input_.empty(); }
  private:
    std::string input_;
};
#endif
