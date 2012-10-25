#ifndef INC_ANALYSIS_STATISTICS_H
#define INC_ANALYSIS_STATISTICS_H
#include "Analysis.h"
class Analysis_Statistics : public Analysis {
  public:
    // The following 2 are also used in Analysis_Crankshaft
    static const char torsion_ss[][8];
    static const double torsion_offset[6];

    Analysis_Statistics();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Statistics(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,int);
    Analysis::RetType Analyze();
    void Print(DataFileList*) { return; }
  private:
    std::vector<DataSet*> datasets_;
    std::string filename_;
    CpptrajFile outfile_;
    double shift_;
    int debug_;

    static const char pucker_ss[][9];
    void PuckerAnalysis( DataSet*, int );
    void TorsionAnalysis( DataSet*, int );
    static const char distance_ss[][8];
    void DistanceAnalysis( DataSet*, int, double, double );
};
#endif
