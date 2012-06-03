#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H
#include <vector>
#include "Action.h"
// Class: Closest
/// Modify the state so that only the closest solvent molecules are kept.
class Closest: public Action {
  public:
    Closest();
    ~Closest();

    void print();
  private:
    int init();
    int setup();
    int action();

    DataSetList outList_;
    DataFile *outFile_;
    DataSet *framedata_;
    DataSet *moldata_;
    DataSet *distdata_;
    DataSet *atomdata_;

    int Nclosest_;        ///< Index into Closest datasets
    char *prefix_;
    int closestWaters_;
    bool firstAtom_;
    //AtomMask soluteMask_;
    AtomMask stripMask_;
    AtomMask distanceMask_;
    Topology *newParm_;
    //Topology *oldParm_;
    int NsolventMolecules_;
    Frame newFrame_;
    std::vector<int> keptWaterAtomNum_;

    /** The moldist structure is used in order to preserve the original
      * solvent molecule and atom numbers after sorting. */
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
    std::vector<MolDist> SolventMols_;
};
#endif  
