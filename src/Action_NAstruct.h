#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include <vector>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"
/// Basic Nucleic acid structure analysis. 
/** Calculate nucleic acid base/base pair structural parameters.
  * Algorithms for calculation of base/base pair structural parameters
  * adapted from:
  *   Babcock MS, Pednault EPD, Olson WK, "Nucleic Acid Structure Analysis: 
  *   Mathematics for Local Cartesian and Helical Structure Parameters That
  *   Are Truly Comparable Between Structures", J. Mol. Biol. (1994) 237,
  *   125-156.
  * NA base reference frame coordinates taken from:
  *   Olson WK, Bansal M, Burley SK, Dickerson RE, Gerstein M, Harvey SC,
  *   Heinemann U, Lu XJ, Neidle S, Shekked Z, Sklenar H, Suzuki M, Tung CS,
  *   Westhof E, Wolberger C, Berman H, "A Standard Reference Frame for the 
  *   Description of Nucleic Acid Base-pair Geometry", J. Mol. Biol. (2001)
  *   313, 229-237.
  */
class NAstruct: public Action {
    // Variables
    std::vector<AxisType*> RefCoords;
    std::vector<AxisType*> BaseAxes;
    std::vector<AxisType*> BasePairAxes;
    std::vector<AxisType*> ExpFrames;
    std::vector<AtomMask*> ExpMasks;
    std::vector<int> BasePair;           //< [Base1-0,Base2-0,IsAnti-0], [Base1-1...
    int Nbp;                             //< Total number of base pairs
    int Nbases;                          //< Total number of NA bases
    double HBcut2;                       //< Cutoff^2 for determining hydrogen bonds
    double Ocut2;                        //< Cutoff^2 for determining base-pairing vi origins
    // Datasets
    DataSetList SHEAR;
    DataSetList STRETCH;
    DataSetList STAGGER;
    DataSetList BUCKLE;
    DataSetList PROPELLER;
    DataSetList OPENING;
    DataSetList SHIFT;
    DataSetList SLIDE;
    DataSetList RISE;
    DataSetList TILT;
    DataSetList ROLL;
    DataSetList TWIST;
    int Nframe;                          //< Keep track of # frames for print() function           
    // Init Args
    Range resRange;
    char *outFilename;
    char *naoutFilename;
    bool noheader;
    // Functions
    void ClearLists();
    bool GCpair(AxisType *, AxisType *);
    bool ATpair(AxisType *, AxisType *);
    bool basesArePaired(AxisType *, AxisType *);
    int determineBasePairing();
    //int setupBasePairAxes();
    int setupBaseAxes(Frame *);
    int determineBaseParameters();
    int determineBasepairParameters();
  public:
    NAstruct();
    ~NAstruct();

    int init();
    int setup();
    int action();
    void print();
};
#endif
