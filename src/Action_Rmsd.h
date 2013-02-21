#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H
#include <vector>
#include "Action.h"
#include "Range.h"
#include "ReferenceAction.h"
// Class: Action_Rmsd
/// Action to calculate the RMSD between frame and a reference frame.
class Action_Rmsd: public Action, ReferenceAction {
  public:
    Action_Rmsd();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Rmsd(); }
    static void Help();
    ~Action_Rmsd();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    // PerResRMSD -------------
    bool perres_;                      ///< If true calculate per-residue rmsd
    int NumResidues_;                  ///< Total # of residues to calculate per res rmsd for
    std::vector<DataSet*> PerResRMSD_; ///< Hold residue RMSDs
    std::vector<AtomMask> tgtResMask_; ///< Hold target masks for each res in ResRange
    std::vector<AtomMask> refResMask_; ///< Hold reference masks for each res in ResRange
    std::vector<bool> resIsActive_;    ///< True if residue was set up correctly
    Range ResRange_;                   ///< Residues to calculate perRes rmsd for
    Range RefRange_;                   ///< Residues in reference corresponding to those in ResRange
    DataFile* perresout_;              ///< Per res RMSD data output file
    std::string perresmask_;           ///< Additional mask to apply to residue masks
    bool perrescenter_;                ///< Move residues to common COM before rms calc
    bool perresinvert_;                ///< If true rows will contain set info instead of cols
    DataFile* perresavg_;              ///< Hold per residue average filename
    Frame *ResFrame_;                  ///< Hold residue target coords
    Frame *ResRefFrame_;               ///< Hold residue reference coords.
    Topology* RefParm_;                ///< Needed for mask setup in PerResSetup
    // ------------------------ 
    AtomMask FrameMask_;               ///< Target mask.
    bool fit_;                         ///< If false do not calculate best-fit RMSD
    bool rotate_;                      ///< If true and fitting, do not rotate coords.
    bool useMass_;
    Vec3 Trans_;                       ///< For fit, hold translation tgt->origin
    Frame SelectedFrame_;              ///< Hold only selected target frame coords.
    DataSet *rmsd_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    /// Resize per-residue RMSD masks
    void resizeResMasks();
    /// Set up per-residue RMSD calc
    int perResSetup(Topology*,Topology*);
};
#endif
