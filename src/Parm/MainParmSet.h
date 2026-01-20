#ifndef INC_PARM_MAINPARMSET_H
#define INC_PARM_MAINPARMSET_H
#include "../Timer.h"
class ArgList;
class DataSetList;
class DataSet_Parameters;
namespace Cpptraj {
namespace Parm {
/// Hold any combination of parameter sets 
class MainParmSet {
  public:
    /// CONSTRUCTOR
    MainParmSet();
    /// DESTRUCTOR
    ~MainParmSet();
    /// Associated parameter keywords for InitMainParmSet
    static const char* parm_keywords_;
    /// Initialize the MainParmSet
    int InitMainParmSet(ArgList&, DataSetList const&, int);
    /// Write timing info to stdout
    void TimingInfo(double, int) const;

    /// \return True if a parameter set is defined
    bool HasMainParmSet() const { return (mainParmSet_ != 0); }
    /// \return Main parm set
    DataSet_Parameters const* MainParmSetPtr() const { return mainParmSet_; }
  private:
    /// Get parameter sets
    int getParameterSets(ArgList&, DataSetList const&);

    DataSet_Parameters* mainParmSet_; ///< Hold optional parameter set.
    int debug_;                       ///< Debug level
    bool free_parmset_mem_;           ///< True if main parm set is combined and should be freed
    Timer t_total_;
};
}
}
#endif
