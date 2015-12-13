#ifndef INC_ACTION_ATOMICCORR_H
#define INC_ACTION_ATOMICCORR_H
#include "Action.h"
class Action_AtomicCorr : public Action {
  public:
    Action_AtomicCorr();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AtomicCorr(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    class AtomVector {
      public:
        AtomVector() : idx_(0) {}
        AtomVector(std::string const& sIn, int idxIn) : lbl_(sIn), idx_(idxIn) {}
        void push_back(float fval)                 { vec_.push_back( fval ); }
        int operator-(const AtomVector& rhs) const { return idx_ - rhs.idx_; }
        bool empty()                         const { return vec_.empty();    }
        size_t size()                        const { return vec_.size();     }
        std::string const& Label()           const { return lbl_;            }
        Vec3 VXYZ(int idx) const { return Vec3((double)vec_[idx  ], 
                                               (double)vec_[idx+1], 
                                               (double)vec_[idx+2]); }
      private:
        std::vector<float> vec_;
        std::string lbl_;
        int idx_;
    };

    enum AcorrModeType { ATOM = 0, RES };
    static const char* ModeString[];
    AcorrModeType acorr_mode_;
    double cut_;
    int min_;
    int debug_;
    DataSet* dset_;
    DataFile* outfile_;

    typedef std::vector< AtomVector > ACvector;
    ACvector atom_vectors_;
    AtomMask mask_;
    std::vector<AtomMask> resmasks_;
    Frame refframe_;
};
#endif
