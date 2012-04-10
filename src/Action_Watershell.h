#ifndef INC_ACTION_WATERSHELL_H
#define INC_ACTION_WATERSHELL_H
#include "Action.h"
/// Calculate number of solvent residues in 1st/2nd solvation shell.
class Watershell : public Action {
  public:
    Watershell();

  private:
    AtomMask soluteMask_;
    AtomMask solventMask_;
    char *solventmaskexpr_;
    std::vector<int> activeResidues_;
    std::vector<int> lower_;
    std::vector<int> upper_;
    //int visits_;
    double lowerCutoff_;
    double upperCutoff_;
    char *filename_;

    int init();
    int setup();
    int action();
    void print();
};
#endif
