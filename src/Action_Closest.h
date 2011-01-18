#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H

// Should automatically include AmberParm.h from Action.h
#include "Action.h"
#include <vector>

class Closest: public Action {
    DataFile *outFile;
    DataSetList *outList;
    DataSet *framedata;
    DataSet *moldata;
    DataSet *distdata;
    DataSet *atomdata;
    int Nclosest;

    int closestWaters;
    bool noimage, firstAtom;
    int imageType; 
    AtomMask Mask1;
    AtomMask *tempMask;
    AmberParm *newParm;
    AmberParm *oldParm;
    Frame *newFrame;
    struct MolDist {
      int mol;        // Solvent molecule number
      double D;       // Solvent molecule distance
      AtomMask *mask; // Solvent molecule atom mask
    };
    struct moldist_cmp {
      bool operator()(MolDist first, MolDist second) const {
        if (first.D < second.D)
          return true;
        else
          return false;
      }
    };
    std::vector<AtomMask*> MaskList;

    void ClearMaskList();    
  public:
    Closest();
    ~Closest();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
