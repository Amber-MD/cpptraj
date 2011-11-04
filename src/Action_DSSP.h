#ifndef INC_ACTION_DSSP_H
#define INC_ACTION_DSSP_H
#include "Action.h"
// Class: DSSP
/// Calculate protein secondary structure using DSSP algorithm.
class DSSP : public Action {
    enum SStype { SECSTRUCT_NULL, SECSTRUCT_PARA, SECSTRUCT_ANTI, 
                  SECSTRUCT_3_10, SECSTRUCT_ALPHA, SECSTRUCT_PI, 
                  SECSTRUCT_TURN 
                };
    static const char SSchar[];
    static const char SSname[][6];
    struct Residue {
      double SSprob[7];         ///< Hold probabilities for each SS
      SStype sstype;            ///< Assigned secondary structure
      bool isSelected;          ///< True if C, H, N, and O atom of res selected 
      int C;                    ///< atom idx of BB carbon
      int O;                    ///< atom idx of BB oxygen
      int N;                    ///< atom idx of BB nitrogen
      int H;                    ///< atom idx of BB hydrogen
      std::vector<int> CO_HN_Hbond;  ///< Size=# residues, 1 if this residue CO-HN hbonded to res X
      DataSet *resDataSet;
    };
    std::vector<Residue> SecStruct;

    char *outfilename; ///< Data file name
    DataSet *dssp;    ///< If printString, hold the string dataset
    AtomMask Mask;    ///< Mask used to determine selected residues
    int Nres;         ///< Current total # of residues
    double Nframe;    ///< For calculating SS avg. NOTE: Should be passed in somehow?
    //CpptrajFile debugout; // DEBUG
    std::string sumOut;     ///< File to output SS avgs (dssp.dat.sum)
    char *SSline;     ///< Hold SS propensity for frame, each char represents a residue
    bool printString; ///< If true print 1 char per residue indicating ss type

    // For printString=false, Int dataset, hold SStype for each residue at each frame
    DataSetList *SSdata; ///< Hold data for SS assignment each frame
    DataSetList *dsspData; ///< Used to set up datasets for averaging SS

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
