#ifndef INC_EXEC_CLUSTERMAP_H
#define INC_EXEC_CLUSTERMAP_H
#include "Exec.h"
#include "DataSet_MatrixFlt.h"
#ifdef TIMER
#include "Timer.h"
#endif
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
    void AddCluster(Iarray const&, DataSet_2D const&, DataSet_MatrixFlt&);

    class Cluster;
    typedef std::vector<Cluster> Carray;
    Carray clusters_;
#   ifdef _OPENMP
    std::vector<Iarray> thread_neighbors_;
    int mythread_;
#   endif

    double epsilon_;
    double epsilon2_;
    double Avg_;
    int minPoints_;
    int nClusters_;
#   ifdef TIMER
    Timer t_overall_;
    Timer t_query1_;
    Timer t_query2_;
#   endif
};
// -----------------------------------------------
class Exec_ClusterMap::Cluster {
  public:
    Cluster() : avg_(0.0), cnum_(-1), min_col_(-1), max_col_(-1), min_row_(-1), max_row_(-1) {}
    Cluster(Iarray const& p, double A, int C, int minc, int maxc, int minr, int maxr) :
      points_(p), avg_(A), cnum_(C), min_col_(minc), max_col_(maxc), min_row_(minr), max_row_(maxr)
      {}
    Iarray const& Points() const { return points_; }
    double Avg()           const { return avg_; }
    int Cnum()             const { return cnum_; }
    int MinCol()           const { return min_col_; }
    int MaxCol()           const { return max_col_; }
    int MinRow()           const { return min_row_; }
    int MaxRow()           const { return max_row_; }
  private:
    Iarray points_;
    double avg_;
    int cnum_;
    int min_col_;
    int max_col_;
    int min_row_;
    int max_row_;
};
#endif
