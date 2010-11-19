#ifndef INC_ACTION_HBOND_H
#define INC_ACTION_HBOND_H

#include "Action.h"
#include <vector>
#include <map>

class Hbond : public Action {
    struct HbondType {
      int A;
      int H;
      int D;
      int Frames;
      double dist;
      double angle;
    };

    int Nframes;
    char *avgout;

    std::map<int,HbondType> HbondMap;
    std::vector<int> Donor;
    std::vector<int> Acceptor;
    AtomMask Mask;
    std::vector<int>::iterator accept;
    std::vector<int>::iterator donor;
    double acut;
    double dcut;

    DataSet *NumHbonds;

    struct hbond_cmp {
      bool operator()(HbondType first, HbondType second) const {
        if (first.Frames > second.Frames)
          return true;
        else
          return false;
      }
    };

    //bool HbondSort( HbondType, HbondType);

  public:
    Hbond();
    ~Hbond();

    int init();
    int setup();
    int action();
    void print();
};

#endif
