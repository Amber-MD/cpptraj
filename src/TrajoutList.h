#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include "Trajout.h"
#include "TopologyList.h"
// Class: TrajoutList
/// Hold trajectories for output
class TrajoutList : public FileList {
  public:
    TrajoutList();
    ~TrajoutList();
    static bool CheckCommand(ArgList&);
    int AddEnsembleTrajout(ArgList const&, TopologyList&, int);
    /// Add a traj file to the list with given access and associate with a parm
    int AddTrajout(ArgList&, TopologyList&);
    /// Call write for all trajectories
    int Write(int, Topology*, Frame*);
    /// Call end for all trajectories
    void Close();
    void Info();
  private:
    typedef std::vector<Trajout*> ListType;
    ListType trajout_;

    int AddTrajout(std::string const&, ArgList&, TopologyList&);
};
#endif

