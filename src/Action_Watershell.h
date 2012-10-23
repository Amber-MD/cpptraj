#ifndef INC_ACTION_WATERSHELL_H
#define INC_ACTION_WATERSHELL_H
#include "Action.h"
#include "ImagedAction.h"
/// Calculate number of solvent residues in 1st/2nd solvation shell.
class Action_Watershell : public Action, ImagedAction {
  public:
    Action_Watershell();

  private:
    AtomMask soluteMask_;
    AtomMask solventMask_;
    std::string solventmaskexpr_;
    std::vector<int> activeResidues_;
    double lowerCutoff_;
    double upperCutoff_;
    DataSet* lower_;
    DataSet* upper_;

    int init();
    int setup();
    int action();
    void print() {}
};
#endif
