#ifndef INC_ACTION_HBOND_H
#define INC_ACTION_HBOND_H
#include <vector>
#include <map>
#include "Action.h"
// Class: Action_Hbond
/// Action to calculate the Hbonds present in each frame.
class Action_Hbond : public Action {
  public:
    Action_Hbond();
    ~Action_Hbond();

    void print();
  private:
    int init();
    int setup();
    int action();

    struct HbondType {
      int A;        ///< Acceptor atom#
      int H;        ///< Hydrogen atom#
      int D;        ///< Donor atom#
      int Frames;   ///< # frames this hbond has been present
      double dist;  ///< Used to calc avg distance of this hbond
      double angle; ///< Used to calc avg angle of this hbond
    };

    int Nframes_;
    char *avgout_;
    std::map<int,HbondType> HbondMap_; ///< Track all hbonds found
    typedef std::vector<int> HBlistType;
    HBlistType Donor_;                 ///< Array of hbond donor atoms (D0, H0, D1, H1, ...)
    HBlistType Acceptor_;              ///< Array of hbond acceptor atoms (A0, A1, ...)
    AtomMask Mask_;
    AtomMask DonorMask_;
    AtomMask AcceptorMask_;
    bool hasDonorMask_;
    bool hasAcceptorMask_;
    double acut_;
    double dcut2_;

    DataSet *NumHbonds_;
    DataSetList *HBavg_;
    /// Return true if the first hbond has more frames than the second.
    struct hbond_cmp {
      inline bool operator()(HbondType first, HbondType second) const {
        if (first.Frames > second.Frames)
          return true;
        else
          return false;
      }
    };

    void SearchAcceptor(AtomMask&,bool);
    void SearchDonor(AtomMask&,bool);
    
};
#endif
