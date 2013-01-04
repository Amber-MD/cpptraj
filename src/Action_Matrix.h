#ifndef INC_ACTION_MATRIX_H
#define INC_ACTION_MATRIX_H
#include "Action.h"
#include "DataSet_Matrix.h"
#include "DataSet_Vector.h"
#include "ActionFrameCounter.h"
/// Calculate various types of matrices
class Action_Matrix : public Action, ActionFrameCounter {
  public:
    Action_Matrix();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Matrix(); }
    static void Help();

    void Print();
  private:
    enum OutputType { BYATOM=0, BYRESIDUE, BYMASK };

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    DataSet_Matrix* Mat_;
    DataFile* outfile_;
    AtomMask mask1_;
    AtomMask mask2_;
    std::string filename_;
    OutputType outtype_;
    DataSet_Matrix::MatrixType type_;
    int snap_;
    // IRED only
    int order_;
    std::vector<DataSet_Vector*> IredVectors_;
    // MWcovar only
    std::vector<double> mass1_;
    std::vector<double> mass2_;
    
    bool useMask2_;
    bool useMass_;
    Topology* CurrentParm_; // For ByResidue output

    void FillMassArray(Topology*, std::vector<double>&, AtomMask&);
    void CalcIredMatrix();
    void CalcDistanceMatrix(Frame*);
    void StoreVec(DataSet_Matrix::iterator&,DataSet_Matrix::iterator&,const double*);
    void CalcCovarianceMatrix(Frame*);
    void CalcIdeaMatrix(Frame*);
    void CalcCorrelationMatrix(Frame*);
    void CalcDistanceCovarianceMatrix(Frame*);

    void FinishCovariance();
    void FinishCorrelation();
    void FinishDistanceCovariance();
};
#endif
