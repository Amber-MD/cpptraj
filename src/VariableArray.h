#ifndef INC_VARIABLEARRAY_H
#define INC_VARIABLEARRAY_H
#include "ArgList.h"
/// Hold script variables and their values.
class VariableArray {
  public:
    /// CONSTRUCTOR
    VariableArray() {}
    /// Add/update variable with given value.
    void UpdateVariable(std::string const&, std::string const&);
    /// Add/append variable with given value.
    void AppendVariable(std::string const&, std::string const&);
    /// Replace all variables in given ArgList with their values.
    ArgList ReplaceVariables(ArgList const&);
    /// Print all variable/value pairs to stdout
    void PrintVariables() const;
  private:
    /// Hold variable and corresponding value.
    typedef std::pair<std::string, std::string> Vpair;
    /// Hold list of variables and corresponding values.
    typedef std::vector<Vpair> Varray;

    Varray CurrentVars_; ///< Hold all current variables
};
#endif
