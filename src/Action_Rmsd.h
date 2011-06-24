#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H

#include "Action.h"
#include <vector>
#include "Range.h"
#include "TrajectoryFile.h"

class Rmsd: public Action {
    DataSet *rmsd;
    // PerResRMSD -------------
    DataSetList *PerResRMSD;
    int nres;
    Range ResRange;                     // Residues to calculate perRes rmsd for
    std::vector<AtomMask*> PerResMask;  // Hold AtomMasks for each residue in ResRange
    Range RefRange;                     // Residues in reference corresponding to those in ResRange
    std::vector<AtomMask*> PerRefMask;  // Hold AtomMasks for each residue in RefRange
    char *perresout;                    // PerRes RMSD data output file
    char *perresmask;                   // Additional mask to apply to 
    bool perrescenter;                  // Move residues to common COM before rms calc
    bool perresinvert;                  // If true rows will contain set info instead of cols
    Frame *ResFrame;                    // Hold residue coords
    Frame *ResRefFrame;                 // Hold residue reference coords.
    DataFile *outFile;
    // ------------------------ 
    AtomMask RefMask, FrameMask;        // Frame and reference masks.
    bool nofit, first, perres, useMass; // Action options
    Frame *RefFrame;                    // Hold reference frame coords
    Frame *SelectedRef;                 // Hold only ref coods selected by maskRef
    Frame *SelectedFrame;               // Hold only frame coords selected by mask0
    TrajectoryFile *RefTraj;            // Reference trajectory, each frame used in turn
    AmberParm *RefParm;                 // Reference frame Parm

    int SetRefMask();

  public:
    Rmsd();
    ~Rmsd();

    int init();
    int setup();
    int action();
    void print();
};
#endif
