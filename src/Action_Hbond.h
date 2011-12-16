#ifndef INC_ACTION_HBOND_H
#define INC_ACTION_HBOND_H
#include <vector>
#include <map>
#include "Action.h"
// Class: Hbond
/// Action to calculate the Hbonds present in each frame.
class Hbond : public Action {
    struct HbondType {
      int A;        ///< Acceptor atom#
      int H;        ///< Hydrogen atom#
      int D;        ///< Donor atom#
      int Frames;   ///< # frames this hbond has been present
      double dist;  ///< Used to calc avg distance of this hbond
      double angle; ///< Used to calc avg angle of this hbond
    };

    int Nframes;
    char *avgout;

    std::map<int,HbondType> HbondMap;
    std::vector<int> Donor;
    std::vector<int> Acceptor;
    AtomMask Mask;
    AtomMask DonorMask;
    AtomMask AcceptorMask;
    bool hasDonorMask;
    bool hasAcceptorMask;
    std::vector<int>::iterator accept;
    std::vector<int>::iterator donor;
    double acut;
    double dcut;
    double dcut2;

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
