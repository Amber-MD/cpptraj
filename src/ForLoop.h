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
    ForLoop() : loopvar_(0) {}
    virtual ~ForLoop() {}

    /// Set up the loop, ensure syntax is correct
    virtual int SetupFor(CpptrajState&, ArgList&)  = 0;
    /// Start the loop
    /** \return Number of iterations in loop, or NITERATIONS_UNKNOWN. */
    virtual int BeginFor(DataSetList const&) = 0;
    /// \return True if loop is done, otherwise increment the loop
    virtual bool EndFor(DataSetList&) = 0;

    /// Return value when number of iterations is not known
    static const int NITERATIONS_UNKNOWN = -1;
    /// Return value when an error has occurred
    static const int LOOP_ERROR = -2;

    /// \return For loop description TODO make virtual fn?
    std::string const& Description() const { return description_; }
    /// \return For loop variable name
    std::string const& VarName() const;
    /// \return true if For loop is set up and has a variable
    bool IsSetup() const { return (loopvar_ != 0); }
  protected:
    /// Set loop description. Should be called inside SetupFor.
    void SetDescription(std::string const& descIn) { description_ = descIn; }
    /// Set up the loop variable with given name. Should be called inside SetupFor
    int SetupLoopVar(DataSetList&, std::string const&);

    //ForType VarType()            const { return varType_; }
  private:
    std::string description_; ///< For loop long description
    DataSet* loopvar_;        ///< Variable over which loop is iterating.
};
#endif
