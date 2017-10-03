#ifndef INC_VARIABLEARRAY_H
#define INC_VARIABLEARRAY_H
#include "ArgList.h"
/// Hold script variables and their values.
class VariableArray {
  public:
    /// Hold variable and corresponding value.
    typedef std::pair<std::string, std::string> Vpair;
    /// Hold list of variables and corresponding values.
    typedef std::vector<Vpair> Varray;
    /// CONSTRUCTOR
    VariableArray() {}
    /// Add variable with initial value.
    void AddVariable(std::string const&, std::string const&);
    /// Replace all variables in given ArgList with their values.
    ArgList ReplaceVariables(ArgList const&);
  private:
    Varray CurrentVars_; ///< Hold all current variables
};
#endif
