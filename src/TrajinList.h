#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "Trajin.h"
#include "TopologyList.h"
// Class: TrajinList
/// Hold input trajectories
class TrajinList {
  public:
    enum TrajModeType { UNDEFINED, NORMAL, ENSEMBLE };
    TrajinList();
    ~TrajinList();
    void SetDebug(int dIn) { debug_ = dIn; }
    TrajModeType Mode() { return mode_; }
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(ArgList&, TopologyList&);
    /// Set up frames to be processed 
    int SetupFrames();

    typedef std::vector<Trajin*> ListType;
    typedef ListType::const_iterator const_iterator;
    const_iterator begin() { return trajin_.begin(); }
    const_iterator end()   { return trajin_.end();   }
  private:
    ListType trajin_;
    int debug_;
    TrajModeType mode_; ///< Trajectory processing mode 
};
#endif

