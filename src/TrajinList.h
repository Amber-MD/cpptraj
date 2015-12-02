#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "Trajin.h"
#include "EnsembleIn.h"
/// Hold input trajectories
class TrajinList {
    typedef std::vector<Trajin*> tListType;
    typedef std::vector<EnsembleIn*> eListType;
  public:
    enum TrajModeType { UNDEFINED = 0, NORMAL, ENSEMBLE };
    TrajinList();
    ~TrajinList();
    void Clear();
    void SetDebug(int dIn) { debug_ = dIn; }
    /// Add a trajectory file to the list.
    int AddTrajin(std::string const&, Topology*, ArgList const&);
    /// Add an ensemble to the list.
    int AddEnsemble(std::string const&, Topology*, ArgList const&);

    typedef tListType::const_iterator trajin_it;
    trajin_it trajin_begin() const { return trajin_.begin(); }
    trajin_it trajin_end()   const { return trajin_.end();   }

    typedef eListType::const_iterator ensemble_it;
    ensemble_it ensemble_begin() const { return ensemble_.begin(); }
    ensemble_it ensemble_end()   const { return ensemble_.end();   }
    void FirstEnsembleReplicaInfo() const {
      if (!ensemble_.empty()) ensemble_.front()->PrintReplicaInfo();
    }

    bool empty()         const { return trajin_.empty() && ensemble_.empty(); }
    TrajModeType Mode()  const { return mode_; }
    int MaxFrames()      const { return maxframes_; }
    int TopFrames(int i) const { return topFrames_[i]; }
    unsigned int Size()  const { return trajin_.size() + ensemble_.size(); }
    std::vector<int> const& PindexFrames() const { return topFrames_; }
    void List() const;
  private:
    void UpdateMaxFrames(InputTrajCommon const&);

    tListType trajin_;
    eListType ensemble_;
    int debug_;
    int maxframes_;
    TrajModeType mode_;
    typedef std::vector<int> Iarray;
    Iarray topFrames_; ///< Record how many frames currently associated with each topology.
    /// CRDIDXARG: Used when processing ensemble and sorting by CRDIDX
    std::string finalCrdIndicesArg_;
};
#endif
