#ifndef INC_DIFFUSIONRESULTS_H
#define INC_DIFFUSIONRESULTS_H
#include <string>
// Fwd declare
class DataSet;
class DataFile;
class DataSetList;
class DataFileList;
namespace Cpptraj {
/// Calculate and hold calculated diffusion constants from MSD data
class DiffusionResults {
  public:
    DiffusionResults();

    int AddDiffOut(DataFileList&, std::string const&);

    int CreateDiffusionSets(DataSetList&, std::string const&);

    void Info() const;

    void CalcDiffusionConst(unsigned int&, DataSet*, int, std::string const&) const;

    DataFile* DiffOut() const { return diffout_; }
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
