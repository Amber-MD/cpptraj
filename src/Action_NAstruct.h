#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// Action_NAstruct
#include <vector>
#include <map>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"
#include "DataSet_1D.h"
class Trajout_Single;
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
  * Default conventions for determining base pairing etc are those used in 
  * 3DNA:
  *   Lu XJ, Olson WK. "3DNA: a software package for the analysis, rebuilding
  *   and visualization of three-dimensional nucleic acid structures".
  *   Nucleic Acids Res. 2003 Sep 1;31(17):5108-21.
  *   doi: 10.1093/nar/gkg680. PMID: 12930962; PMCID: PMC212791.
  */
class Action_NAstruct: public Action {
  public:
    Action_NAstruct();
    /// DESTRUCTOR - Needed in case axes trajectories are being written
    ~Action_NAstruct();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_NAstruct(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    void NA_Sync( DataSet_1D*, std::vector<int> const&, std::vector<int> const&, int, int ) const;
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print();
    // ----- Enumerations ------------------------
    enum HbondType { WC = 0, HOOG, OTHER };
    enum GrooveType { PP_OO = 0, HASSAN_CALLADINE };
    enum BP_ConventionType { BP_3DNA = 0, BP_BABCOCK };
    /// How to find base pairs: first frame, reference structure, all frames.
    enum FindType { FIRST = 0, REFERENCE, ALL, SPECIFIED };
    // ----- Data Structures ---------------------
    /// Hold consecutive bases
    struct Stype {
      DataSet_1D* dx_;
      DataSet_1D* dy_;
      DataSet_1D* dz_;
      DataSet_1D* rx_;
      DataSet_1D* ry_;
      DataSet_1D* rz_;
      unsigned int strandidx_; ///< Index in Strands_
      //unsigned int base1idx_; ///< Index of first base in Bases_
      //unsigned int base2idx_; ///< Index of second base in Bases_
    };
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
      DataSet_1D* isBP_;
      DataSet_1D* major_;
      DataSet_1D* minor_;
#     ifdef NASTRUCTDEBUG
      //DataSet*    axes_oxyz_; ///< Base pair axes origin vectors
      DataSet*    axes_nxyz_; ///< Base pair axes Z (i.e. mean normal) vectors
#     endif
      unsigned int bpidx_;
      unsigned int base1idx_; ///< Index of first base in Bases_
      unsigned int base2idx_; ///< Index of second base in Bases_
      int nhb_;               ///< Current # of hydrogen bonds in base pair.
      int n_wc_hb_;           ///< Number of WC hydrogen bonds in base pair.
      bool isAnti_;           ///< True if base Z axes are not aligned
      bool isZ_;              ///< True if Z axes are aligned 3' to 5' instead of 5' to 3'
    };
    /// Hold a base pair step.
    struct StepType {
      NA_Axis stepaxis_; ///< Base pair step axis
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
      DataSet_1D* Zp_;
      DataSet_1D* majGroove_; ///< H-C major groove width
      DataSet_1D* minGroove_; ///< H-C minor groove width
      unsigned int b1idx_; ///< Index of base pair 1 base 1 in Bases_
      unsigned int b2idx_; ///< Index of base pair 1 base 2
      unsigned int b3idx_; ///< Index of base pair 2 base 1
      unsigned int b4idx_; ///< Index of base pair 2 base 2
      // For Hassan-Calladine groove calc
      int P_m2_; ///< Phosphate 2 towards 5', strand 1
      int p_p2_; ///< Phosphate 2 towards (5' anti, 3' para), strand 2
      int P_p1_; ///< Phosphate 1 towards 3'. strand 1
      int P_p2_; ///< Phosphate 2 towards 3', strand 1
      int p_m1_; ///< Phosphate 1 towards (3' anti, 5' para), strand 2
      int p_m2_; ///< Phosphate 2 towards (3' anti, 5' para), strand 2
    };
    // ----- Type Definitions --------------------
    typedef std::vector<NA_Base> Barray;      ///< Array of NA bases
    typedef std::pair<int,int> Rpair;         ///< Pair of residue numbers / BP indices
    typedef std::map<Rpair,Stype> Smap;       ///< Map of residue numbers to strand pairs
    typedef std::map<Rpair,BPtype> BPmap;     ///< Map of residue numbers to BP
    typedef std::map<Rpair,StepType> StepMap; ///< Map of BP indices to Steps
    typedef std::vector<Rpair> StrandArray;   ///< Hold indices into Bases_ for strand beg/end
    typedef std::map<std::string, NA_Base::NAType> RefMapType; ///< Map custom res names to target types
    typedef std::vector< std::pair<unsigned int, unsigned int> > PairArray; ///< Specified base pair #s
    // ----- Functions ---------------------------
    /// Initialize an axes pseudo-trajectory
    int init_axes_pseudoTraj(const char*, const char*, const char*,
                             const char*, const char*, 
                             DataSetList const&, ArgList&,
                             Trajout_Single**, Topology**) const;
    /// Set up an axes pseudo-trajectory
    int setup_axes_pseudoTraj(Topology&, Trajout_Single&, Frame&,
                              std::vector<Residue> const&) const;
    /// Recursively travel to 3' terminal base
    static int follow_base_to_3prime(Barray&, unsigned int, std::vector<bool>&, int);
    /// Recursively travel sugar-phosphate backbone to find the next residue in a strand.
    static int TravelBackbone(Topology const&, int, std::vector<int>&);
    /// Set up axes for each base.
    int SetupBaseAxes(Frame const&);
    /// Identify hydrogen bonding between G and C
    static HbondType GCpair(NA_Base const&, int, NA_Base const&, int);
    /// Identify hydrogen bonding between A and T/U
    static HbondType ATpair(NA_Base const&, int, NA_Base const&, int);
    /// Identify the type of hydrogen bonding between two bases.
    static HbondType ID_HBtype(NA_Base const&, int, NA_Base const&, int);
    /// Calculate the total number of hydrogen bonds and WC hbonds between two bases.
    int CalcNumHB(NA_Base const&, NA_Base const&, int&);
    /// \return New/existing base pair corresponding to given bases.
    BPmap::iterator AddBasePair(int, NA_Base const&, int, NA_Base const&);
    /// Determine which bases are paired geometrically, set base pair data.
    int DetermineBasePairing();
    /// Find index in bases for given internal residue #
    int find_index_in_bases(int) const;
    /// Set up base pairs based on user specification
    int SpecifiedBasePairing();
    /// Calculate translational/rotational parameters between two axes.
    int calculateParameters(NA_Axis const&, NA_Axis const&, NA_Axis*, double*);
    /// Calculate helical parameters between two axes.
    int helicalParameters(NA_Axis const&, NA_Axis const&, double *);
    /// \return index of base that is N steps away in specified direction from another base.
    int GetBaseIdxStep(int, int) const;
    /// Determine individual base parameters in single strands.
    int DetermineStrandParameters(int);
    /// Check that base Z axis points 5' to 3'
    int axis_points_5p_to_3p(NA_Base const&) const;
    /// Determine individual base and base pair parameters.
    int DeterminePairParameters(int);
    /// Determine base pair steps and step parameters, including HC groove calc.
    int DetermineStepParameters(int);
    /// Update each time series to nframes_;
    void UpdateSeries();
    /// Calculate number of WC hydrogen bonds for each base pair.
    inline void CalculateHbonds();
    /// Set up data sets for StepType (except HC groove)
    MetaData NewStepType(StepType&, int, int, int, int, int) const;
    // ----- Variables ---------------------------
    NA_Reference refBases_;             ///< Hold reference bases
    RefMapType nameToRef_;              ///< Map residue names to custom reference
    Barray Bases_;                      ///< Hold nucleobases
    Smap StrandPairs_;                  ///< Hold consecutive bases in strands
    BPmap BasePairs_;                   ///< Hold base pairs
    StepMap Steps_;                     ///< Hold base pair steps.
    StrandArray Strands_;               ///< Hold strand info
    PairArray specifiedPairs_;          ///< User-specified base pairing
    NA_Base::PmethodType puckerMethod_; ///< Pucker calculation method.
    double HBdistCut2_;                 ///< distance Cutoff^2 for determining hydrogen bonds
    double originCut2_;                 ///< Cutoff^2 for determining base-pairing vi origins
    double staggerCut_;                 ///< Cutoff for determining base vertical separation
    double z_angle_cut_;                ///< Cutoff for angle between base Z-axes
    int maxResSize_;                    ///< Max residue size, used to set up frames for RMS fit.
    int debug_;
    int nframes_;                       ///< Total number of frames calculated.
    FindType findBPmode_;               ///< How base pairs are to be found.
    GrooveType grooveCalcType_;         ///< Type of groove calc to perform
    BP_ConventionType bpConvention_;    ///< Conventions to use when determining base pairing.
    Range resRange_;                    ///< Range to search for NA residues.
    bool printheader_;                  ///< If true, print header to naout files.
    bool seriesUpdated_;                ///< If false, check that time series data is nframes long
    bool skipIfNoHB_;                   ///< When true, do not calc parameters when BP not present
    bool spaceBetweenFrames_;           ///< If false do not print spaces between frames in naout
    bool sscalc_;                       ///< If true determine params for consecutive bases in strands
    bool wc_hb_only_;                   ///< If true, only report # of WC hydrogen bonds.
    CpptrajFile* bpout_;                ///< Base pair out (BP.<suffix>).
    CpptrajFile* ssout_;                ///< Single strand out (SS.<suffix>).
    CpptrajFile* stepout_;              ///< Base pair step out (BPstep.<suffix>).
    CpptrajFile* helixout_;             ///< Helical parameters out (Helix.<suffix>).
    std::string dataname_;              ///< NA DataSet name (default NA).
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
#   ifdef NASTRUCTDEBUG
    // DEBUG - used to trigger AxisPDBwriter for first call of calculateParameters
    bool calcparam_;
#   endif
    Trajout_Single* axesOut_;     ///< Output trajectory for base axes
    Topology* axesParm_;          ///< Pseudo-topology for base axes
    Frame axesFrame_;             ///< Frame for base axes pseudo traj
    Trajout_Single* bpAxesOut_;   ///< Output trajectory for base pair axes
    Topology* bpAxesParm_;        ///< Pseudo-topology for base pair axes
    Frame bpAxesFrame_;           ///< Frame for base pair axes pseudo traj
    Trajout_Single* stepAxesOut_; ///< Output trajectory for base pair step axes
    Topology* stepAxesParm_;      ///< Pseudo-topology for base pair step axes
    Frame stepAxesFrame_;         ///< Frame for base pair step axes pseudo traj
    int setupNframes_;            ///< Set in Setup(); number of expected frames to write (pseudo-traj)
    Topology* setupTop_;          ///< Set in Setup(); current topology
    NameType axisNameO_;
    NameType axisNameX_;
    NameType axisNameY_;
    NameType axisNameZ_;
};
#endif
