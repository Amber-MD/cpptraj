#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H

#include "Action.h"
#include <vector>


class Rmsd: public Action {
    DataSet *rmsd;
    // PerResRMSD -------------
    DataSetList *PerResRMSD;
    int nres;
    std::list<int> *ResRange;           // Residues to calculate perRes rmsd for
    std::vector<AtomMask*> PerResMask;  // Hold AtomMasks for each residue in ResRange
    std::list<int> *RefRange;           // Residues in reference corresponding to those in ResRange
    std::vector<AtomMask*> PerRefMask;  // Hold AtomMasks for each residue in RefRange
    char *perresout;                    // PerRes RMSD data output file
    char *perresmask;                   // Additional mask to apply to 
    bool perrescenter;                  // Move residues to common COM before rms calc
    bool perresinvert;                  // If true rows will contain set info instead of cols
    // ------------------------ 
    AtomMask RefMask, FrameMask;        // Frame and reference masks.
    bool nofit, first, perres, useMass; // Action options
    Frame *RefFrame;                    // Hold reference frame coords
    Frame *SelectedRef;                 // Hold only ref coods selected by maskRef
    Frame *SelectedFrame;               // Hold only frame coords selected by mask0
    AmberParm *RefParm;                 // Reference frame Parm

    int SetRefMask();

  public:
    Rmsd();
    ~Rmsd();

    int init();
    int setup();
    int action();

};
#endif
