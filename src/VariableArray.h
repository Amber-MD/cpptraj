#ifndef INC_VARIABLEARRAY_H
#define INC_VARIABLEARRAY_H
#include "DataSetList.h"
/// Hold script variables and their values.
class VariableArray {
    /// Hold variable and corresponding value.
    typedef std::pair<std::string, std::string> Vpair;
    /// Hold list of variables and corresponding values.
    typedef std::vector<Vpair> Varray;
  public:
    /// CONSTRUCTOR
    VariableArray() {}
    /// Add/update variable with given value.
    void UpdateVariable(std::string const&, std::string const&);
    /// Add/append variable with given value.
    void AppendVariable(std::string const&, std::string const&);
    /// Replace all variables in given ArgList with their values.
    ArgList ReplaceVariables(ArgList const&, DataSetList const&, int);
    /// Print all variable/value pairs to stdout
    void PrintVariables() const;

    typedef Varray::const_iterator const_iterator;
    const_iterator begin() const { return CurrentVars_.begin(); }
    const_iterator end()   const { return CurrentVars_.end(); }
  private:
    Varray CurrentVars_; ///< Hold all current variables
};
#endif
