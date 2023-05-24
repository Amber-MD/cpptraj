#ifndef INC_DIFFUSIONRESULTS_H
#define INC_DIFFUSIONRESULTS_H
#include <string>
// Fwd declare
class DataSet;
class DataFile;
class DataSetList;
class DataFileList;
class ArgList;
namespace Cpptraj {
/// Calculate and hold calculated diffusion constants from MSD data
class DiffusionResults {
  public:
    DiffusionResults();

    int InitDiffusionResults(DataSetList&, DataFileList&, ArgList& actionArgs, std::string const&);

    void Info() const;

    void CalcDiffusionConst(unsigned int&, DataSet*, int, std::string const&) const;
  private:
    DataSet* diffConst_;   ///< Hold diffusion constants.
    DataSet* diffLabel_;   ///< Hold diffusion constant labels.
    DataSet* diffSlope_;   ///< Hold MSD vs time line slopes.
    DataSet* diffInter_;   ///< Hold MSD vs time line intercepts.
    DataSet* diffCorrl_;   ///< Hold MSD vs time line correlation.
    DataFile* diffout_;    ///< Diffusion results output file.
};
}
#endif
