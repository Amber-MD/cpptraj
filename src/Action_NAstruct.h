#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// Action_NAstruct
#include <vector>
#include <map>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"
#include "DataSet_1D.h"
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
class Action_NAstruct: public Action {
  public:
    Action_NAstruct();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NAstruct(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    // Functions
    static int TravelBackbone(Topology const&, int, std::vector<int>&); 
    int SetupBaseAxes(Frame const&);

    //int GCpair(NA_Base const&, NA_Base const&);
    //int ATpair(NA_Base const&, NA_Base const&);
    int CalcNumHB(NA_Base const&, NA_Base const&);
    int DetermineBasePairing();

    int calculateParameters(NA_Axis const&, NA_Axis const&, NA_Axis*, double*);
    int helicalParameters(NA_Axis const&, NA_Axis const&, double *);
    int DeterminePairParameters(int);
    void CalcPucker(NA_Base&, int); // TODO: Move to NA_Base
    int DetermineStepParameters(int);

    typedef std::vector<NA_Base> Barray;
    Barray Bases_;        ///< Hold nucleobases
    /// Hold a base pair.
    struct BPtype {
      NA_Axis bpaxis_; ///< Base pair axis.
      DataSet_1D* shear_;
      DataSet_1D* stretch_;
      DataSet_1D* stagger_;
      DataSet_1D* buckle_;
      DataSet_1D* prop_;
      DataSet_1D* opening_;
      DataSet_1D* hbonds_;
      DataSet_1D* major_;
      DataSet_1D* minor_;
      unsigned int bpidx_;
      unsigned int base1idx_; ///< Index of first base in Bases_
      unsigned int base2idx_; ///< Index of second base in Bases_
      int nhb_;               ///< Current # of hydrogen bonds in base bpair.
      bool isAnti_;
    };
    typedef std::pair<int,int> Rpair; ///< Residue pair
    typedef std::map<Rpair,BPtype> BPmap;
    BPmap BasePairs_;     ///< Hold base pairs
    /// Hold a base pair step.
    struct StepType {
      DataSet_1D* shift_;
      DataSet_1D* slide_;
      DataSet_1D* rise_;
      DataSet_1D* tilt_;
      DataSet_1D* roll_;
      DataSet_1D* twist_;
      DataSet_1D* xdisp_;
      DataSet_1D* ydisp_;
      DataSet_1D* hrise_;
      DataSet_1D* incl_;
      DataSet_1D* tip_;
      DataSet_1D* htwist_;
      unsigned int b1idx_; ///< Index of base pair 1 base 1 in Bases_
      unsigned int b2idx_;
      unsigned int b3idx_;
      unsigned int b4idx_;
    };
    typedef std::map<Rpair,StepType> StepMap;
    StepMap Steps_;       ///< Hold base pair steps.
    // Variables
    NA_Base::PmethodType puckerMethod_;
    double HBdistCut2_;                 ///< distance Cutoff^2 for determining hydrogen bonds
    double originCut2_;                 ///< Cutoff^2 for determining base-pairing vi origins
    int maxResSize_;                    ///< Max residue size, used to set up frames for RMS fit.
    int debug_;
    int nframes_;
    Range resRange_;                    ///< Range to search for NA residues.
    bool printheader_;                  ///< If true, print header to naout files.
    bool useReference_;                 ///< If true, use reference to determine base pairing.
    CpptrajFile* bpout_;                ///< Base pair out (BP.<suffix>).
    CpptrajFile* stepout_;              ///< Base pair step out (BPstep.<suffix>).
    CpptrajFile* helixout_;             ///< Helical parameters out (Helix.<suffix>).
    std::string dataname_;              ///< NA DataSet name (default NA).
    /// Map a residue name to an NA base type.
    typedef std::map<std::string, NA_Base::NAType> ResMapType;
    ResMapType CustomMap_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
#   ifdef NASTRUCTDEBUG
    // DEBUG - used to trigger AxisPDBwriter for first call of calculateParameters
    bool calcparam_;
#   endif
};
#endif
