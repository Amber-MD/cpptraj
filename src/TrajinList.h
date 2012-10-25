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
    static void Help();
    void SetDebug(int dIn) { debug_ = dIn; }
    TrajModeType Mode() { return mode_; }
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(ArgList&, TopologyList&);

    typedef std::vector<Trajin*> ListType;
    typedef ListType::const_iterator const_iterator;
    const_iterator begin() { return trajin_.begin(); }
    const_iterator end()   { return trajin_.end();   }
    void List();
    int MaxFrames() { return maxframes_; }

  private:
    ListType trajin_;
    int debug_;
    int maxframes_;
    TrajModeType mode_; ///< Trajectory processing mode 
};
#endif

