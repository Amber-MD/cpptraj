#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H
#include <vector>
#include "Action.h"
#include "Range.h"
#include "TrajectoryFile.h"
// Class: Rmsd
/// Action to calculate the RMSD between frame and a reference frame.
class Rmsd: public Action {
  public:
    DataSet *rmsd;
  private:
    // PerResRMSD -------------
    bool perres;                      ///< If true calculate per-residue rmsd
    int NumResidues;                  ///< Total # of residues to calculate per res rmsd for
    DataSetList *PerResRMSD;          ///< Hold residue RMSDs
    std::vector<AtomMask> tgtResMask; ///< Hold target masks for each res in ResRange
    std::vector<AtomMask> refResMask; ///< Hold reference masks for each res in ResRange
    std::vector<bool> resIsActive;    ///< True if residue was set up correctly
    Range ResRange;                   ///< Residues to calculate perRes rmsd for
    Range RefRange;                   ///< Residues in reference corresponding to those in ResRange
    char *perresout;                  ///< Per res RMSD data output file name
    char *perresmask;                 ///< Additional mask to apply to residue masks
    bool perrescenter;                ///< Move residues to common COM before rms calc
    bool perresinvert;                ///< If true rows will contain set info instead of cols
    char *perresavg;                  ///< Hold per residue average filename
    Frame *ResFrame;                  ///< Hold residue target coords
    Frame *ResRefFrame;               ///< Hold residue reference coords.
    // ------------------------ 
    AtomMask RefMask;                 ///< Reference mask
    AtomMask FrameMask;               ///< Frame masks.
    bool nofit;                       ///< If true do not calculate best-fit RMSD
    bool first;                       ///< If true use first frame read in as reference
    double Trans[6];                  ///< For fit, hold 2 translations: tgt->origin, origin->ref
    Frame RefFrame;                   ///< Hold reference frame coords
    Frame SelectedRef;                ///< Hold only ref coods selected by maskRef
    Frame SelectedFrame;              ///< Hold only frame coords selected by mask0
    TrajectoryFile *RefTraj;          ///< Reference trajectory, each frame used in turn
    AmberParm *RefParm;               ///< Reference frame Parm
    /// Set reference mask based on reference parm
    int SetRefMask();
    /// Set ref structure from ref mask and pre-center if fitting
    void SetRefStructure();
    /// Resize per-residue RMSD masks
    void resizeResMasks();
    /// Set up per-residue RMSD calc
    int perResSetup();

  public:
    Rmsd();
    ~Rmsd();

    int SeparateInit(char*,bool,int);

    int init();
    int setup();
    int action();
    void print();
};
#endif
