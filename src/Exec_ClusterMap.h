#ifndef INC_EXEC_CLUSTERMAP_H
#define INC_EXEC_CLUSTERMAP_H
#include "Exec.h"
#include "DataSet_MatrixFlt.h"
// EXPERIMENTAL ALPHA CODE
class Exec_ClusterMap : public Exec {
  public:
    Exec_ClusterMap();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ClusterMap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;

    void RegionQuery(Iarray&, double, int, DataSet_2D const&);
    void AddCluster(Iarray const&, DataSet_MatrixFlt&);

    double epsilon_;
    double epsilon2_;
    double Avg_;
    int minPoints_;
    int nClusters_;
};
#endif
