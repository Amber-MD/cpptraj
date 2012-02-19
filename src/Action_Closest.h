#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H
#include <vector>
#include "Action.h"
// Class: Closest
/// Modify the state so that only the closest solvent molecules are kept.
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
    bool firstAtom;
    AtomMask soluteMask;
    AtomMask stripMask;
    AmberParm *newParm;
    AmberParm *oldParm;
    Frame newFrame;
    /// The moldist structure is used in order to preserve the original
    /// solvent molecule and atom numbers after sorting.
    struct MolDist {
      int mol;        ///< Original solvent molecule number
      double D;       ///< Closest distance of solvent molecule to atoms in soluteMask
      AtomMask mask;  ///< Original solvent molecule atom mask
    };
    /// Return true if the first molecule is closer than the second
    struct moldist_cmp {
      inline bool operator()(MolDist first, MolDist second) const {
        if (first.D < second.D)
          return true;
        else
          return false;
      }
    };
    std::vector<MolDist> SolventMols;

    //void ClearMaskList();    
  public:
    Closest();
    ~Closest();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
