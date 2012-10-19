#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H
#include <vector>
#include "Action.h"
#include "Range.h"
#include "TrajectoryFile.h"
// Class: Action_Rmsd
/// Action to calculate the RMSD between frame and a reference frame.
class Action_Rmsd: public Action {
  public:
    Action_Rmsd();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Rmsd(); }
    static void Help();

    ~Action_Rmsd();

    void print();
  private:
    int init();
    int setup();
    int action();

    // PerResRMSD -------------
    bool perres_;                      ///< If true calculate per-residue rmsd
    int NumResidues_;                  ///< Total # of residues to calculate per res rmsd for
    std::vector<DataSet*> PerResRMSD_; ///< Hold residue RMSDs
    std::vector<AtomMask> tgtResMask_; ///< Hold target masks for each res in ResRange
    std::vector<AtomMask> refResMask_; ///< Hold reference masks for each res in ResRange
    std::vector<bool> resIsActive_;    ///< True if residue was set up correctly
    Range ResRange_;                   ///< Residues to calculate perRes rmsd for
    Range RefRange_;                   ///< Residues in reference corresponding to those in ResRange
    std::string perresout_;            ///< Per res RMSD data output file name
    std::string perresmask_;           ///< Additional mask to apply to residue masks
    bool perrescenter_;                ///< Move residues to common COM before rms calc
    bool perresinvert_;                ///< If true rows will contain set info instead of cols
    std::string perresavg_;            ///< Hold per residue average filename
    Frame *ResFrame_;                  ///< Hold residue target coords
    Frame *ResRefFrame_;               ///< Hold residue reference coords.
    // ------------------------ 
    AtomMask FrameMask_;               ///< Target mask.
    bool nofit_;                       ///< If true do not calculate best-fit RMSD
    bool rotate_;                      ///< If true and fitting, do not rotate coords.
    bool useMass_;
    double Trans_[6];                  ///< For fit, hold 2 translations: tgt->origin, origin->ref
    Frame SelectedFrame_;              ///< Hold only selected target frame coords.
    DataSet *rmsd_;
    // Reference variables and functions
    enum RefModeType { UNKNOWN_REF=0, FIRST, REF, REFTRAJ };
    RefModeType refmode_;
    Frame RefFrame_;
    Topology* RefParm_; // Needed for PerResSetup
    Frame SelectedRef_;
    AtomMask RefMask_;
    TrajectoryFile RefTraj_;
    int SetRefMask( Topology* );
    void SetRefStructure( Frame& );

    /// Resize per-residue RMSD masks
    void resizeResMasks();
    /// Set up per-residue RMSD calc
    int perResSetup(Topology*);
};
#endif
