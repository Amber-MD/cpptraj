#ifndef INC_ACTION_ATOMICCORR_H
#define INC_ACTION_ATOMICCORR_H
#include "Action.h"
class Action_AtomicCorr : public Action {
  public:
    Action_AtomicCorr();
    
    void print();
  private:
    int init();
    int setup();
    int action();

    typedef std::vector< std::vector<float> > ACvector;
    ACvector atom_vectors_;
    AtomMask mask_;
    std::string outname_;
    Frame refframe_;
};
#endif
