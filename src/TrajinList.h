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
    void Clear();
    void SetDebug(int dIn) { debug_ = dIn; }
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(ArgList&, TopologyList&);

    typedef std::vector<Trajin*> ListType;
    typedef ListType::const_iterator const_iterator;
    const_iterator begin() const { return trajin_.begin(); }
    const_iterator end()   const { return trajin_.end();   }
    TrajModeType Mode()    const { return mode_;           }
    const Trajin* front()  const { return trajin_.front(); }
    int MaxFrames()        const { return maxframes_;      }
    void List() const;
  private:
    ListType trajin_;
    int debug_;
    int maxframes_;
    TrajModeType mode_; ///< Trajectory processing mode 
};
#endif

