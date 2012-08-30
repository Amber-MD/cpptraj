#ifndef INC_ACTION_ATOMICCORR_H
#define INC_ACTION_ATOMICCORR_H
#include "Action.h"
class Action_AtomicCorr : public Action {
  public:
    Action_AtomicCorr();
    
    void print();
  private:
    class AtomVector {
      public:
        AtomVector() : idx_(0) {}
        AtomVector(std::string const& sIn, int idxIn) : lbl_(sIn), idx_(idxIn) {}
        int operator-(const AtomVector& rhs) { return idx_ - rhs.idx_; }
        bool empty()   { return vec_.empty(); }
        size_t size()  { return vec_.size();  }
        void push_back(float fval) { vec_.push_back( fval ); }
        void XYZ(double* Vout, int idx3) {
          Vout[0] = (double)vec_[idx3  ];
          Vout[1] = (double)vec_[idx3+1];
          Vout[2] = (double)vec_[idx3+2];
        }
        std::string const& Label() { return lbl_; }
      private:
        std::vector<float> vec_;
        std::string lbl_;
        int idx_;
    };

    int init();
    int setup();
    int action();

    enum AcorrModeType { ATOM = 0, RES };
    static const char ModeString[][8];
    AcorrModeType acorr_mode_;
    double cut_;
    int min_;

    typedef std::vector< AtomVector > ACvector;
    ACvector atom_vectors_;
    AtomMask mask_;
    std::vector<AtomMask> resmasks_;
    std::string outname_;
    Frame refframe_;
};
#endif
