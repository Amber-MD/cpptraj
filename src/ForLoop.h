#ifndef INC_FORLOOP_H
#define INC_FORLOOP_H
#include <string>
class CpptrajState;
class ArgList;
class DataSetList;
class DataSet;
/// Abstract base class for all for loops
class ForLoop {
  public:
    //ForLoop() : varType_(UNKNOWN), Niterations_(-1) {}
    ForLoop() {}
    virtual ~ForLoop() {}
    /// Set up the loop, ensure syntax is correct
    virtual int SetupFor(CpptrajState&, std::string const&, ArgList&)  = 0;
    /// Start the loop
    /** \return Number of iterations in loop, or NITERATIONS_UNKNOWN. */
    virtual int BeginFor(DataSetList const&) = 0;
    /// \return True if loop is done, otherwise increment the loop
    virtual bool EndFor(DataSetList const&) = 0;
    /// Return value when number of iterations is not known
    static const int NITERATIONS_UNKNOWN = -1;
    /// Return value when an error has occurred
    static const int LOOP_ERROR = -2;

    std::string const& Description() const { return description_; }
    //std::string const& VarName() const { return varname_; }
    std::string const& VarName() const;
  protected:
    void SetDescription(std::string const& descIn) { description_ = descIn; }
    //void SetType(ForType f)                        { varType_ = f;          }
    /// Set up the loop variable with given name.
    int SetupLoopVar(DataSetList&, std::string const&);

    //ForType VarType()            const { return varType_; }
  private:
    std::string description_; ///< For loop long description
    DataSet* loopvar_;        ///< Variable over which loop is iterating.
    //std::string varname_;     ///< Variable over which for loop is iterating
    //ForType varType_;         ///< For loop type
};
#endif
