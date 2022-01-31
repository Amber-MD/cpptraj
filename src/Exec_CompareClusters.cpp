#include "Exec_CompareClusters.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

/** CONSTRUCTOR. */
Exec_CompareClusters::Exec_CompareClusters() :
  Exec(GENERAL)
{
  SetHidden( true );
}

/** Get cluster number vs time data set. */
DataSet* Exec_CompareClusters::getClusterSet(std::string const& dsname, DataSetList const& DSL) {
  DataSet* ds = DSL.GetDataSet( dsname );
  if (ds == 0) {
    mprinterr("Error: '%s' does not correspond to a data set.\n", dsname.c_str());
    return 0;
  }
  if (ds->Group() != DataSet::SCALAR_1D) {
    mprinterr("Error: '%s' is not a 1D scalar data set.\n", dsname.c_str());
    return 0;
  }
  return ds;
}

/** For tracking same/different frames in clusters. */
class ClusterStat {
  public:
    ClusterStat() : n_in_common_(0) {}
    /// Both sets had frame assigned to same cluster
    void FrameInCommon() { n_in_common_++; }
    /// Increment counter for specified set only
    void IncrementSet(unsigned int i) {
      if (i >= n_in_set_.size())
        n_in_set_.resize(i + 1, 0);
      n_in_set_[i]++;
    }
    /// Header to stdout
    static void Header() {
      mprintf("%-12s %12s %12s %12s\n", "#Cluster", "Common", "Set0", "Set1");
    }
    /// Print to stdout
    void Print(unsigned int cnum) {
      if (n_in_set_.size() < 2) n_in_set_.resize(2, 0);
      mprintf("%12i %12u %12u %12u\n", cnum, n_in_common_, n_in_set_[0], n_in_set_[1]);
    }
  private:
    typedef std::vector<unsigned int> Uarray;

    unsigned int n_in_common_; ///< Number of frames in common
    Uarray n_in_set_;          ///< Number of frames just in specified set
};

/** Compare cluster # vs time sets. */
int Exec_CompareClusters::CompareClusters(DataSet* c0, DataSet* c1) const {
  std::vector<ClusterStat> clusterStats;

  DataSet_1D const& set0 = static_cast<DataSet_1D const&>( *c0 );
  DataSet_1D const& set1 = static_cast<DataSet_1D const&>( *c1 );

  unsigned int n0 = set0.Size();
  unsigned int n1 = set1.Size();

  if (n0 != n1) {
    mprintf("Warning: Size of '%s' (%u) != size of '%s' (%u)\n",
            c0->legend(), n0, c1->legend(), n1);
  }

  unsigned int noise0 = 0;
  unsigned int noise1 = 0;
  unsigned int idx0 = 0;
  unsigned int idx1 = 0;
  while (idx0 < n0 && idx1 < n1)
  {
    // TODO check for fractional component
    int cnum0 = (int)set0.Dval(idx0);
    int cnum1 = (int)set1.Dval(idx1);

    if (cnum0 == -1) noise0++;
    if (cnum1 == -1) noise1++;

    if (cnum0 != -1 && cnum0 == cnum1) {
      // This frame was assigned the same cluster in both sets
      if ((unsigned int)cnum0 >= clusterStats.size())
        clusterStats.resize(cnum0+1);
      clusterStats[cnum0].FrameInCommon();
    } else {
      // This frame was assigned different clusters
      if (cnum0 != -1) {
        if ((unsigned int)cnum0 >= clusterStats.size())
          clusterStats.resize(cnum0+1);
        clusterStats[cnum0].IncrementSet( 0 );
      }
      if (cnum1 != -1) {
        if ((unsigned int)cnum1 >= clusterStats.size())
          clusterStats.resize(cnum1+1);
        clusterStats[cnum1].IncrementSet( 1 );
      }
    }

    idx0++;
    idx1++;
  }

  mprintf("\tSet 0: %u noise frames\n", noise0);
  mprintf("\tSet 1: %u noise frames\n", noise1);

  ClusterStat::Header();
  for (std::vector<ClusterStat>::iterator it = clusterStats.begin();
                                          it != clusterStats.end(); ++it)
    it->Print(it - clusterStats.begin());
  return 0;
}

// Exec_CompareClusters::Help()
void Exec_CompareClusters::Help() const
{
  mprintf("\tset <set0> set <set1>\n");
}

// Exec_CompareClusters::Execute()
Exec::RetType Exec_CompareClusters::Execute(CpptrajState& State, ArgList& argIn)
{
  DataSet* c0 = getClusterSet( argIn.GetStringKey("set"), State.DSL() );
  if (c0 == 0) return CpptrajState::ERR;
  DataSet* c1 = getClusterSet( argIn.GetStringKey("set"), State.DSL() );
  if (c1 == 0) return CpptrajState::ERR;
  mprintf("\tCluster set 0: '%s'\n", c0->legend());
  mprintf("\tCluster set 1: '%s'\n", c1->legend());

  if (CompareClusters(c0, c1)) return CpptrajState::ERR;

  return CpptrajState::OK;
}
