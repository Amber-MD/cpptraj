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
    void Clear();
    int AddEnsembleTrajout(ArgList const&, TopologyList const&, int);
    /// Add a traj file to the list with given access and associate with a parm
    int AddTrajout(ArgList const&, TopologyList const&);
    /// Call write for all trajectories
    int Write(int, Topology*, Frame*);
    /// Call end for all trajectories
    void Close();
    void List() const;
    // The definitions below are for ensemble processing.
    typedef std::vector<ArgList> ArgsArray;
    typedef std::vector<ArgList>::const_iterator ArgIt;
    ArgIt argbegin() const { return trajoutArgs_.begin(); }
    ArgIt argend()   const { return trajoutArgs_.end();   }
  private:
    typedef std::vector<Trajout*> ListType;
    ListType trajout_;
    /// Array of trajout args for setting up ensemble trajout.
    ArgsArray trajoutArgs_;

    int AddTrajout(std::string const&, ArgList&, TopologyList const&);
};
#endif
