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
    char *prefix;
    int closestWaters;
    bool noimage, firstAtom;
    int imageType; 
    AtomMask Mask1;
    AtomMask *tempMask;
    AmberParm *newParm;
    AmberParm *oldParm;
    Frame *newFrame;
    // The moldist structure is used in order to preserve the original
    // solvent molecule and atom numbers after sorting.
    struct MolDist {
      int mol;        // Original solvent molecule number
      double D;       // Closest distance of solvent molecule to atoms in Mask1
      AtomMask *mask; // Original solvent molecule atom mask
    };
    // Return true if the first molecule is closer than the second
    struct moldist_cmp {
      inline bool operator()(MolDist first, MolDist second) const {
        if (first.D < second.D)
          return true;
        else
          return false;
      }
    };
    std::vector<MolDist> SolventMols;

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
