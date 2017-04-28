#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "Trajin.h"
#include "EnsembleIn.h"
/// Hold input trajectories
class TrajinList {
    typedef std::vector<Trajin*> tListType;
    typedef std::vector<EnsembleIn*> eListType;
  public:
    TrajinList();
    ~TrajinList();
    void Clear();
    void SetDebug(int dIn) { debug_ = dIn; }
    /// Add a trajectory file to the list.
    int AddTrajin(std::string const&, Topology*, ArgList const&);
    /// Add an ensemble to the list.
    int AddEnsembleIn(std::string const&, Topology*, ArgList const&);

    typedef tListType::const_iterator trajin_it;
    trajin_it trajin_begin() const { return trajin_.begin(); }
    trajin_it trajin_end()   const { return trajin_.end();   }

    typedef eListType::const_iterator ensemble_it;
    ensemble_it ensemble_begin() const { return ensemble_.begin(); }
    ensemble_it ensemble_end()   const { return ensemble_.end();   }
    void FirstEnsembleReplicaInfo() const {
      if (!ensemble_.empty()) ensemble_.front()->PrintReplicaInfo();
    }
#   ifdef MPI
    EnsembleIn* EnsPtr(int idx) { return ensemble_[idx]; }
#   endif
    bool empty()         const { return trajin_.empty() && ensemble_.empty(); }
    int MaxFrames()      const { return maxframes_; }
    int TopFrames(int i) const { return topFrames_[i]; }
    int EnsembleSize()   const { return ensembleSize_; }
    unsigned int Size()  const { return trajin_.size() + ensemble_.size(); }
    std::vector<int> const& PindexFrames() const { return topFrames_; }
    void List() const;
  private:
    void UpdateMaxFrames(InputTrajCommon const&);

    tListType trajin_;
    eListType ensemble_;
    int debug_;
    int maxframes_;
    typedef std::vector<int> Iarray;
    Iarray topFrames_; ///< Record how many frames currently associated with each topology.
    /// CRDIDXARG: Used when processing ensemble and sorting by CRDIDX
    std::string finalCrdIndicesArg_;
    /// Current ensemble size. All input ensembles in given run must be same size.
    int ensembleSize_;
};
#endif
