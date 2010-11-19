#ifndef INC_ACTION_DSSP_H
#define INC_ACTION_DSSP_H

#include "Action.h"
#include <vector>

// From ptraj actions.c:transformSecStruct 0.42*0.20*332
#define DSSP_fac 27.888
class DSSP : public Action {
    enum SStype { SECSTRUCT_NULL, SECSTRUCT_PARA, SECSTRUCT_ANTI, 
                  SECSTRUCT_3_10, SECSTRUCT_ALPHA, SECSTRUCT_PI, 
                  SECSTRUCT_TURN 
                };
    static const char SSchar[];
    static const char SSname[][30];
    struct Residue {
      SStype sstype;            // Assigned secondary structure
      bool isSelected;          // true if residue is being considered for SS calc
      int C;                    // atom # of BB carbon
      int O;                    // atom # of BB oxygen
      int N;                    // atom # of BB nitrogen
      int H;                    // atom # of BB hydrogen
      std::vector<int> CO_HN_Hbond;  // Size=# residues, 1 if this residue CO-HN hbonded to res X
      double SSprob[7];         // Hold probabilities for each SS
    };
    std::vector<Residue> SecStruct;

    DataSet *dssp;
    AtomMask Mask;
    int Nres;
    double Nframe; // DEBUG
    PtrajFile debugout; // DEBUG
    char *sumOut;               // File to output SS avgs (dssp.dat.sum)
    char *SSline;               // Hold SS propensity for frame, each char represents a residue

    int isBonded(int, int);
    void SSassign(int, int, SStype);
  public:
    DSSP();
    ~DSSP();

    int init();
    int setup();
    int action();
    void print();
};

#endif
