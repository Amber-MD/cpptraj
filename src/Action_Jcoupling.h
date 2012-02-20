#ifndef INC_ACTION_JCOUPLING_H
#define INC_ACTION_JCOUPLING_H
#include <vector>
#include <map>
#include <string>
#include "Action.h"
#include "Name.h"
// Class: Jcoupling
/// Calculate j-coupling values based on dihedrals using the Karplus equation.
/*! \author Original Code: J. Maier, Stony Brook 2011.
    \author Adapted by: D. Roe, Rutgers 2011.

    Two types of J-couplings are calculated by this code. If the type is
    0, the form used is that described by Chou et al. JACS (2003) 125 
    p.8959-8966 and the four constants represent A, B, C, and D. If the
    type is 1 (denoted C in $AMBERHOME/Karplus.txt) the form use is that 
    described by Perez et al. JACS (2001) 123 p.7081-7093 and the first 
    three constants represent C0, C1, and C2.
 */
class Jcoupling: public Action {
    /// Hold Karplus parameters for a dihedral
    struct karplusConstant {
      NAME atomName[4]; ///< Name of each atom involved in dihedral
      int offset[4];    ///< Offsets
      double C[4];      ///< Constants
      int type;         ///< Calculation type (0=Chou, 1=Perez)
    };
    /// Hold all Karplus constants for a given residue
    typedef std::vector<karplusConstant> *karplusConstantList;
    /// Hold Karplus constants for all residues
    std::map<std::string,karplusConstantList> KarplusConstants;
    /// Hold info for single j-coupling calculation
    struct jcouplingInfo {
      int residue; ///< Residue number
      int atom[4]; ///< Atom #s of the diehdral
      double *C;   ///< Pointer to C in associated karplusConstant structure
      int type;    ///< Calculation type (0=Chou, 1=Perez)
    };
    /// Hold info for all j-coupling calcs
    std::vector<jcouplingInfo> JcouplingInfo;

    AtomMask Mask1;
    int Nconstants;
    // DEBUG
    CpptrajFile outputfile;
    /// Load Karplus parameters from a file
    int loadKarplus(std::string);
  public:
    Jcoupling();
    ~Jcoupling();

    int init();
    int setup();
    int action();
};
#endif  
