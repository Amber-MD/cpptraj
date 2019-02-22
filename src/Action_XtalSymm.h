#ifndef INC_ACTION_XTALSYMM_H
#define INC_ACTION_XTALSYMM_H
#include "Action.h"
#include "Matrix_3x3.h"
#include "ReferenceAction.h"
#include "Vec3.h"

/** XtalSymm: an action to superimpose symmetry-related parts of a simulation system using
  *           crystallographic symmetry operations
  * \author David S. Cerutti
  */
class Action_XtalSymm : public Action {
  public:
    Action_XtalSymm() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_XtalSymm(); }
    void Help() const;
    ~Action_XtalSymm();
  private:
    // Forward declarations
    class XtalDock;
    class TransOp;

    static const int IASU_GRID_BINS_;
    static const double DASU_GRID_BINS_;

    // The space group
    int nops_;              ///< The number of symmetry operations in the canonical space group
    int nCopyA_;            // The number of unit cell replicas along each of three dimensions.
    int nCopyB_;            //   This generalizes the space group into supercells based on that
    int nCopyC_;            //   space group, and allows us to talk about symmetry operations
                            //   in terms of the supercell.
    std::string spaceGrp_;
    int sgID_;              ///< The space group identification number
 
    // Reference frame
    AtomMask tgtMask_;
    ReferenceAction REF_;  
    Frame RefFrame_;       ///< The reference (or first trajectory) frame
    bool useFirst_;        ///< Flag to use first trajectory frame, true if no reference frame

    // Re-imaging for solvent atoms
    int nMolecule_;        ///< Total number of molecules in the system, taken from the topology
    bool allToFirstASU_;   /**< Flag to have all atoms not explicitly in a designated asymmetric
                                unit re-imaged to the primary ASU volume. */
    bool molCentToASU_;    /**< Flag to use molecule centroids, not individual atoms, in the
                                above re-imaging. */
    int* molLimits_;       /**< Start and end points for each molecule (all molecules are assumed
                                to be contiguous within the topology, but it is not assumed
                                that molecule i+1 starts where molecule i ends). */
    bool* molInSolvent_;   /**< Flags to indicate whether each molecule is part of the non-ASU,
                                free-floating "solvent" component. */
  
    // Masks for the asymmetric units and solvent particles
    int nmasks_;
    AtomMask* Masks_; // TODO should probably use vector here
    AtomMask  SolventParticles_;
    AtomMask  SolventMolecules_;
    std::vector<int> subunitOpID_;
  
    // Rotation matrices and translation vectors
    Matrix_3x3* R;
    Matrix_3x3* Rinv;
    Vec3* T;
    Vec3* RefT;
    bool* rotIdentity;

    // Grid for ASU assignment of loose molecules
    TransOp* AsuGrid;
  
    // Methods
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    Action::RetType LoadSpaceGroupSymOps();
    bool OperationAvailable(XtalDock* leads, int* HowToGetThere, int ncurr);
    bool OriginsAlign(XtalDock* leads, int* HowToGetThere, int ncurr);
    void BestSuperposition(int, int, XtalDock*, int&);
    Vec3 BestOrigin(Frame&, Frame*, std::vector<int>&);
    TransOp DetectAsuResidence(double x, double y, double z);
    void BuildAsuGrid();
    void Print() {}
    void ClearMemory();
    double dmin(double, double);
    double dmin(double, double, double);
    double dmin(double, double, double, double);
    double dmax(double, double);
    double dmax(double, double, double);
    double dmax(double, double, double, double);
    bool PointInPrimaryASU(double x, double y, double z);
};

//---------------------------------------------------------------------------------------------
/** XtalDock: stores the results of crystal symmetry operations on solitary subunits,
  *           attempting to align them back to the original asymmetric unit.  The goal
  *           is to cluster the origins of multiple operations such that one cluster
  *           has all of its required origins in about the same spot and gives
  *           consistently low atom positional RMSDs.
  */
class Action_XtalSymm::XtalDock {
  public:
    int subunit_;  ///< The subunit that this is operating on
    int opID_;     ///< The operation in use, indexing the lists R and T in Action_XtalSymm
    double rmsd_;  ///< The un-fitted rmsd between the original and superimposed subunits
    Vec3 displc_;  /**< The optimal displacement beteen the two subunits' centers of mass,
                        scaled to simulation cell fractional coordinates. */
    Vec3 origin_;  ///< The origin that got the best rmsd

    /// Assignment
    XtalDock& operator=(const XtalDock& rhs) {
      if (this == &rhs) return *this;
      subunit_ = rhs.subunit_;
      opID_    = rhs.opID_;
      rmsd_    = rhs.rmsd_;
      displc_  = rhs.displc_;
      origin_  = rhs.origin_;
      return *this;
    }
};

//---------------------------------------------------------------------------------------------
/** TransOp: stores an initial translation (as three doubles) plus the index of an operation.
  *          This information will help guide points (or entire coordinate sets) from some
  *          starting position into one of the regions of space for a given asymmetric unit.
  */
class Action_XtalSymm::TransOp {
  public:
    int opID_;      ///< Operation index
    double tr_x_;   ///< Initial X translation
    double tr_y_;   ///< Initial Y translation
    double tr_z_;   ///< Initial Z translation
};
#endif
