#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include "Trajout_Single.h"
#include "EnsembleOut.h"
/// Hold output trajectories for a run.
/** When a 'trajout' or equivalent command is given, the trajectory will
  * be added to TrajoutList as a Trajout_Single class by default. This is
  * done so that some error checking of arguments etc can occur, even if
  * eventually the list will be used in ensemble output.
  */
class TrajoutList {
  public:
    typedef std::vector<EnsembleOut*> EnsembleArray;

    TrajoutList() : debug_(0) {}
    void SetDebug(int);
    void Clear();
    /// Add output trajectory to the list and associate with given topology.
    int AddTrajout(std::string const&, ArgList const&, Topology*);
    /// \return Array with current output trajectories converted to ensemble output trajectories.
    EnsembleArray MakeEnsembleTrajout() const; // FIXME should this clear the trajout array?
  private:
    int debug_;
    typedef std::vector<Trajout_Single> ListType;
    typedef std::vector<ArgList> ArgsArray;
    typedef std::vector<Topology*> TopArray;
    typedef std::vector<std::string> Sarray;
    ListType trajout_; ///< Hold output trajectories.
    ListType active_;  ///< Hold only active output trajectories.
    ArgsArray trajoutArgs_; ///< Array of trajout args for potentially setting up ensemble.
    TopArray  trajoutTops_; ///< Array of associated topologies.
    Sarray trajoutNames_;   ///< Array of trajout file names.
};
#endif
