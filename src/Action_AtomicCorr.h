#ifndef INC_ACTION_ATOMICCORR_H
#define INC_ACTION_ATOMICCORR_H
#include "Action.h"
/// Calculate the correlation (average dot product) between motion vectors.
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
    /// Hold a series of position vectors ([X1 - X0], [X2 - X1], ..., [XN - XN-1])
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
        std::vector<float> vec_; ///< Array of position vectors for N-1 steps (XYZ)
        std::string lbl_;        ///< Label for this vector array
        int idx_;                ///< Index for this vector array.
    };

    enum AcorrModeType { ATOM = 0, RES };
    static const char* ModeString[];
    AcorrModeType acorr_mode_;
    double cut_;        ///< Do not record correlations less than cut_
    int min_;           ///< Only calculate correlations for residues separated by at least min_
    int debug_;
    DataSet* dset_;     ///< Output matrix data set
    DataFile* outfile_; ///< Output file

    typedef std::vector< AtomVector > ACvector;
    ACvector atom_vectors_;          ///< AtomVector for each atom/residue
    AtomMask mask_;                  ///< Selected atoms for ATOM
    std::vector<AtomMask> resmasks_; ///< Masks for selected residues for RESDIUE
    Frame previousFrame_;            ///< Hold the previous frame
};
#endif
