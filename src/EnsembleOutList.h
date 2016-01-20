#ifndef INC_ENSEMBLEOUTLIST_H
#define INC_ENSEMBLEOUTLIST_H
#include "EnsembleOut.h"
/// Hold output ensembles.
class EnsembleOutList {
  public:
    EnsembleOutList() {}
    ~EnsembleOutList() { Clear(); }
    void Clear();
    int AddEnsembleOut(std::string const&, ArgList const&, Topology*,
                       int, TrajectoryFile::TrajFormatType);
    int SetupEnsembleOut(Topology*, CoordinateInfo const&, int); //TODO const, Topology array?
    int WriteEnsembleOut(int, FramePtrArray const&);
    void CloseEnsembleOut();
    void List(std::vector<int> const&) const;
  private:
    typedef std::vector<Topology*> TopArray;
    typedef std::vector<EnsembleOut*> EnsArray;
    EnsArray ensout_; ///< Array of output ensembles
    EnsArray active_; ///< Array of active output ensembles
    TopArray ensTops_; ///< Array of associated topology files
    std::vector<bool> open_; ///< True if corresponding ensemble is open
};
#endif
