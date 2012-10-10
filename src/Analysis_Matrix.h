#ifndef INC_ANALYSIS_MATRIX_H
#define INC_ANALYSIS_MATRIX_H
#include "Analysis.h"
#include "DataSet_Matrix.h"
#include "ModesInfo.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_Matrix : public Analysis {
  public:
    Analysis_Matrix();
    ~Analysis_Matrix();

    int Setup(DataSetList*);
    int Analyze();
  private:
    DataSet_Matrix* minfo_;
    std::string outfilename_; // carg2
    std::string outthermo_;
    bool thermopt_; // iarg1
    int nevec_; // iarg2
    bool reduce_; // iarg3
    ModesInfo* modinfo_; // carg3

    // Workspace for eigenvector/eigenvalue calcs, stored in ModesInfo.
    bool freeWork_; ///< If true error occured, free workspace vars
    double* eigval_;
    double* vout_; ///< Hold output eigenvectors

    /// Indicate workspace has been stored in ModesInfo 
    void WorkspaceStored();
};
#endif
