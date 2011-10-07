#ifndef INC_ACTION_JCOUPLING_H
#define INC_ACTION_JCOUPLING_H
/// Class: Jcoupling
/// Action to calculate j-coupling values based on dihedrals using the
/// the Karplus equation.
/// Adapted from code written by J. Maier, Stony Brook 2011. 
#include "Action.h"
#include "Name.h"
#include <vector>
#include <map>
#include <string>
class Jcoupling: public Action {
    // Hold Karplus parameters
    struct karplusConstant {
      NAME atomName[4];
      int offset[4];
      double C[4];
      int type;
    };
    typedef std::vector<karplusConstant> *karplusConstantList;
    std::map<std::string,karplusConstantList> KarplusConstants;
    // Hold info for single j-coupling calculation
    struct jcouplingInfo {
      int residue;
      int atom[4];
      double *C;
      int type;
    };
    std::vector<jcouplingInfo> JcouplingInfo;

    AtomMask Mask1;
    int Nconstants;
    // DEBUG
    CpptrajFile outputfile;

    int loadKarplus(char*);
  public:
    Jcoupling();
    ~Jcoupling();

    int init();
    int setup();
    int action();
};
#endif  
