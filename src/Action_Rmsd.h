#ifndef INC_ACTION_RMSD_H
#define INC_ACTION_RMSD_H
#include <vector>
#include "Action.h"
#include "Range.h"
#include "ActionReference.h"
// Class: Rmsd
/// Action to calculate the RMSD between frame and a reference frame.
class Rmsd: public Action, ActionReference {
  public:
    Rmsd();
    ~Rmsd();

    int SeparateInit(char*,bool,int);

    int init();
    int setup();
    int action();
    void print();

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
    AtomMask FrameMask;               ///< Frame masks.
    bool nofit;                       ///< If true do not calculate best-fit RMSD
    double Trans[6];                  ///< For fit, hold 2 translations: tgt->origin, origin->ref
    Frame SelectedFrame;              ///< Hold only frame coords selected by mask0
    DataSet *rmsd;
    /// Resize per-residue RMSD masks
    void resizeResMasks();
    /// Set up per-residue RMSD calc
    int perResSetup(Topology*);
};
#endif
