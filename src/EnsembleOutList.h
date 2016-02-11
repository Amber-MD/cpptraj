#ifndef INC_ENSEMBLEOUTLIST_H
#define INC_ENSEMBLEOUTLIST_H
#include "EnsembleOut.h"
/// Hold output ensembles.
class EnsembleOutList {
  public:
    EnsembleOutList() : debug_(0) {}
    ~EnsembleOutList() { Clear(); }
    void SetDebug(int d) { debug_ = d; }
    void Clear();
    /// Add output ensemble of given size to the list and associate with given Topology
    int AddEnsembleOut(std::string const&, ArgList const&, Topology*, int);
    /// Set up ensembles for given topology.
    int SetupEnsembleOut(Topology*, CoordinateInfo const&, int); //TODO const Topology?
    /// Write array of Frames to ensemble
    int WriteEnsembleOut(int, FramePtrArray const&);
    /// Close all ensembles
    void CloseEnsembleOut();
    /// List output ensembles.
    void List(std::vector<int> const&) const;
    /// \return true if no output ensembles present.
    bool Empty() const { return ensout_.empty(); }
#   ifdef MPI
    int ParallelSetupEnsembleOut(Topology*, CoordinateInfo const&, int, Parallel::Comm const&);
#   endif
  private:
    int debug_;
    typedef std::vector<Topology*> TopArray;
    typedef std::vector<EnsembleOut*> EnsArray;
    EnsArray ensout_;        ///< Array of output ensembles
    EnsArray active_;        ///< Array of active output ensembles
    TopArray ensTops_;       ///< Array of associated topology files
    std::vector<bool> open_; ///< True if corresponding ensemble is open
};
#endif
