#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H
/// Class: Rmsd
/// Action to calculate the RMSD between frame and a reference frame.
#include "Action.h"
#include <vector>
#include "Range.h"
#include "TrajectoryFile.h"
class Rmsd: public Action {
    DataSet *rmsd;
    // PerResRMSD -------------
    int nres;                          // Total # of residues to calculate per res rmsd for
    DataSetList *PerResRMSD;           // Hold residue RMSDs
    std::vector<AtomMask*> tgtResMask; // Hold target masks for each res in ResRange
    std::vector<AtomMask*> refResMask; // Hold reference masks for each res in ResRange
    std::vector<bool> resIsActive;     // True if residue was set up correctly
    Range ResRange;                    // Residues to calculate perRes rmsd for
    Range RefRange;                    // Residues in reference corresponding to those in ResRange
    char *perresout;                   // Per res RMSD data output file name
    char *perresmask;                  // Additional mask to apply to residue masks
    bool perrescenter;                 // Move residues to common COM before rms calc
    bool perresinvert;                 // If true rows will contain set info instead of cols
    Frame *ResFrame;                   // Hold residue target coords
    Frame *ResRefFrame;                // Hold residue reference coords.
    // ------------------------ 
    AtomMask RefMask, FrameMask;        // Frame and reference masks.
    bool nofit, first, perres;          // Action options
    Frame *RefFrame;                    // Hold reference frame coords
    Frame *SelectedRef;                 // Hold only ref coods selected by maskRef
    Frame *SelectedFrame;               // Hold only frame coords selected by mask0
    TrajectoryFile *RefTraj;            // Reference trajectory, each frame used in turn
    AmberParm *RefParm;                 // Reference frame Parm

    int SetRefMask();
    void resizeResMasks();
    int perResSetup();

  public:
    Rmsd();
    ~Rmsd();

    int init();
    int setup();
    int action();
    void print();
};
#endif
