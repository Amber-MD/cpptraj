#ifndef INC_ANALYSIS_MATRIX_H
#define INC_ANALYSIS_MATRIX_H
#include "Analysis.h"
#include "MatrixType.h"
class Analysis_Matrix : public Analysis {
  public:
    Analysis_Matrix();
    ~Analysis_Matrix();

    int Setup(DataSetList*);
    int Analyze();
  private:
    enum ThermoFlag { OFF=0, THERMO, ORDER };

    MatrixType* minfo_;
    std::string outfilename_; // carg2
    std::string orderparamfile_; // carg4
    ThermoFlag thermopt_; // iarg1
    int nevec_; // iarg2
    bool reduce_; // iarg3
    ModesInfo* modinfo_; // carg3

    // Workspace for eigenvector/eigenvalue calcs, stored in ModesInfo.
    bool freeWork_; ///< If true error occured, free workspace vars
    double* vect_;
    double* eigval_;
    double* vout_; ///< Hold output eigenvectors

    /// Indicate workspace has been stored in ModesInfo 
    void WorkspaceStored();
};
#endif
