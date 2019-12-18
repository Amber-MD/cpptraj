#ifndef INC_FORLOOP_H
#define INC_FORLOOP_H
#include <string>
class CpptrajState;
class ArgList;
/// Abstract base class for all for loops
class ForLoop {
  public:
    //ForLoop() : varType_(UNKNOWN), Niterations_(-1) {}
    ForLoop() : Niterations_(-1) {}
    virtual ~ForLoop() {}

    virtual int SetupFor(CpptrajState&, std::string const&, ArgList&)  = 0;
  protected:
    void SetDescription(std::string const& descIn) { description_ = descIn; }
    //void SetType(ForType f)                        { varType_ = f;          }
    void SetVarName(std::string const& v)          { varname_ = v;          }
    void SetNiterations(int n)                     { Niterations_ = n;      }

    std::string const& VarName() const { return varname_; }
    //ForType VarType()            const { return varType_; }
  private:
    std::string description_; ///< For loop long description
    std::string varname_;     ///< Variable over which for loop is iterating
    //ForType varType_;         ///< For loop type
    int Niterations_;         ///< # of iterations the loops will execute if known.
};
#endif
