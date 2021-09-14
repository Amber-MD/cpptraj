#ifndef INC_CLUSTER_OUTPUT_H
#define INC_CLUSTER_OUTPUT_H
class CpptrajFile;
namespace Cpptraj {
namespace Cluster {
class Algorithm;
class Cframes;
class List;
class MetricArray;
class BestReps;
/// Cluster output routines.
class Output {
  public:
    static void PrintClustersToFile(CpptrajFile&, List const&, Algorithm const&, MetricArray&,
                                    int, Cframes const&);
    static int PrintSilhouetteFrames(CpptrajFile&, List const&);
    static int PrintSilhouettes(CpptrajFile&, List const&);
    static int Summary(CpptrajFile&, List const&, Algorithm const&, MetricArray&,
                        bool, bool, Cframes const&);
    static void Summary_Part(CpptrajFile&, unsigned int, Cframes const&, List const&,
                             BestReps const&, MetricArray&, Cframes const&);
  private:
    static unsigned int DetermineNameWidth(List const&);
};

}
}
#endif
