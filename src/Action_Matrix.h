#ifndef INC_ACTION_MATRIX_H
#define INC_ACTION_MATRIX_H
#include "Action.h"
#include "DataSet_Matrix.h"
#include "DataSet_Vector.h"
/// Calculate various types of matrices
class Action_Matrix : public Action {
  public:
    Action_Matrix();
    void print();
  private:
    enum OutputType { BYATOM=0, BYRESIDUE, BYMASK };

    int init();
    int setup();
    int action();

    DataSet_Matrix* Mat_;
    DataFile* outfile_;
    AtomMask mask1_;
    AtomMask mask2_;
    std::string filename_;
    OutputType outtype_;
    DataSet_Matrix::MatrixType type_;
    // TODO: Put start, stop, offset into its own class
    int start_;
    int stop_;
    int offset_;
    // IRED only
    int order_;
    std::vector<DataSet_Vector*> IredVectors_;
    // MWcovar only
    std::vector<double> mass1_;
    std::vector<double> mass2_;
    
    bool useMask2_;

    void FillMassArray(std::vector<double>&, AtomMask&);
    void CalcIredMatrix();
    void CalcDistanceMatrix();
    void StoreVec(DataSet_Matrix::iterator&,DataSet_Matrix::iterator&,const double*);
    void CalcCovarianceMatrix();
    void CalcIdeaMatrix();
    void CalcCorrelationMatrix();
    void CalcDistanceCovarianceMatrix();

    void FinishCovariance();
    void FinishCorrelation();
    void FinishDistanceCovariance();
};
#endif
