#ifndef INC_ACTION_HBOND_H
#define INC_ACTION_HBOND_H
/// Class: Hbond
/// Action to calculate the Hbonds present in each frame.
#include "Action.h"
#include <vector>
#include <map>
class Hbond : public Action {
    struct HbondType {
      int A;        // Acceptor atom#
      int H;        // Hydrogen atom#
      int D;        // Donor atom#
      int Frames;   // # frames this hbond has been present
      double dist;  // Used to calc avg distance of this hbond
      double angle; // Used to calc avg angle of this hbond
    };

    int Nframes;
    char *avgout;

    std::map<int,HbondType> HbondMap;
    std::vector<int> Donor;
    std::vector<int> Acceptor;
    AtomMask Mask;
    AtomMask DonorMask;
    AtomMask AcceptorMask;
    std::vector<int>::iterator accept;
    std::vector<int>::iterator donor;
    double acut;
    double dcut;

    DataSet *NumHbonds;
    DataSetList *HBavg;

    struct hbond_cmp {
      bool operator()(HbondType first, HbondType second) const {
        if (first.Frames > second.Frames)
          return true;
        else
          return false;
      }
    };

    //bool HbondSort( HbondType, HbondType);
    void SearchAcceptor(AtomMask*,bool);
    void SearchDonor(AtomMask*,bool);
    
  public:
    Hbond();
    ~Hbond();

    int init();
    int setup();
    int action();
    void print();
};
#endif
