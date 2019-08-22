#ifndef INC_ACTION_DSSP2_H
#define INC_ACTION_DSSP2_H
#include <vector>
#include <string>
#include "Action.h"
#include "NameType.h"
#include "CharMask.h"
class DataSet;
/// <Enter description of Action_DSSP2 here>
/** Based on protein secondary structure definitions given in:
  *   Kabsch, W.; Sander, C.; \"Dictionary of Protein Secondary Structure:
  *   Pattern Recognition of Hydrogen-Bonded and Geometrical Features.
  *   Biopolymers (1983), V.22, pp.2577-2637.
  * The various secondary structure types that can be assigned (in order
  * of decreasing priority) are:
  *   H - Alpha helix (4)
  *   B - Single bridge (ladder of length 1)
  *   E - Beta strand (ladder length > 1)
  *   G - 3-10 helix (3)
  *   I - Pi helix (5)
  *   T - Turn (3, 4, 5)
  *   S - Bend (Angle i-2, i, i+2 > 70 deg.)
  */
class Action_DSSP2 : public Action {
  public:
    Action_DSSP2();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_DSSP2(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    class ElemHbond;
    class SSres;

    /// Secondary structure types
    enum SStype { NONE=0, EXTENDED, BRIDGE, H3_10, ALPHA, HPI, TURN, BEND };
    static const int NSSTYPE_ = 8;  ///< # of secondary structure types.
    static const char DSSP_char_[]; ///< DSSP 1 char names corresponding to SStype
    static const char* SSname_[];   ///< Full secondary structure names corresponding to SStype
    static const char* SSchar_[];   ///< PTRAJ 1 character names corresponding to SStype
    static const std::string SSzlabels_; ///< Output graph Z labels corresponding to SStype

    /// Elementary hydrogen bond pattern
    enum PatternType {
      NOHBOND = 0, ///< No hydrogen bond
      TURNBEG,     ///< Start of n-Turn, >
      TURNEND,     ///< End of n-Turn, <
      TURNX,       ///< Coincidence of TURNBEG and TURNEND, X
      TURN3,       ///< Inside i to i + 3 turn (T)
      TURN4,       ///< Inside i to i + 4 turn (T)
      TURN5,       ///< Inside i to i + 5 turn (T)
      BRIDGEPARA,  ///< (i-1,j) and (j,i+1) or (j-1,i) and (i,j+1), lower case
      BRIDGEANTI   ///< (i,j) and (j,i) or (i-1,j+1) and (j-1,i+1), upper case
    };

    static const double DSSP_fac_;  ///< Original DSSP factor for calc. H-bond "energy"
    static const double DSSP_cut_;  ///< Original DSSP H-bond energy cutoff in kcal/mol

    typedef std::vector<SSres> SSarrayType;

    int debug_;            ///< Action debug level
    SSarrayType Residues_; ///< Hold SS data for all residues.
    NameType BB_N_;        ///< Protein N atom name ('N')
    NameType BB_H_;        ///< Protein N-H atom name ('H')
    NameType BB_C_;        ///< Protein C atom name ('C')
    NameType BB_O_;        ///< Protein C-O atom name ('O')
    NameType BB_CA_;       ///< Protein alpha C name ('CA')
    CharMask Mask_;        ///< Mask used to determine selected residues.
    DataFile* outfile_;          ///< Output Data file
    DataFile* dsspFile_;         ///< Sum output file
    CpptrajFile* assignout_;     ///< Assignment output file.
    std::string dsetname_;       ///< DSSP data set name
    DataSet* totalDS_[NSSTYPE_]; ///< Hold total SS each frame for each SS type
    ActionInit Init_;            ///< Hold pointers to master DSL/DFL
    bool printString_;           ///< If true print 1 char per residue indicating ss type
};


class Action_DSSP2::SSres {
    typedef std::vector<int> HbArrayType;
  public:
    SSres();
    int Idx()         const { return idx_; }
    bool IsSelected() const { return isSelected_; }
    int C()           const { return C_; }
    int O()           const { return O_; }
    int N()           const { return N_; }
    int H()           const { return H_; }
    int CA()          const { return CA_; }
    DataSet* Dset()   const { return resDataSet_; }

    bool IsMissingAtoms() const { return (C_==-1 || O_==-1 || N_==-1 || H_==-1 || CA_==-1); }
    bool HasCO()          const { return (C_!=-1 && O_!=-1); }
    bool HasNH()          const { return (N_!=-1 && H_!=-1); }
    bool HasCA()          const { return (CA_!=-1); }

    typedef HbArrayType::const_iterator const_iterator;
    const_iterator begin() const { return CO_HN_Hbonds_.begin(); }
    const_iterator end()   const { return CO_HN_Hbonds_.end(); }

    void SetIdx(int i) { idx_ = i; }
    void SetSelected(bool b) { isSelected_ = b; }
    void SetC(int i)         { C_ = i; }
    void SetO(int i)         { O_ = i; }
    void SetN(int i)         { N_ = i; }
    void SetH(int i)         { H_ = i; }
    void SetCA(int i)        { CA_ = i; }
    void SetDset(DataSet* d) { resDataSet_ = d; }
    /// Deselect this residue and reset coordinate indices.
    void Deselect();
    /// Reset hbonds, pattern and SS type assignments
    void Unassign();
    /// Add hbond from this CO to specified NH
    void AddHbond(int i) { CO_HN_Hbonds_.push_back( i ); }
  private:
    HbArrayType CO_HN_Hbonds_; ///< This res C-O hbonded to these res H-N.
    DataSet* resDataSet_;      ///< DataSet for SS assignment each frame for this res.
    double chirality_;         ///< dihedral CA[i-1, i, i+1, i+2]
    int SScount_[NSSTYPE_];    ///< Hold count for each SS type
    PatternType pattern_;      ///< Assigned hbond pattern for this frame
    SStype sstype_;            ///< SS assignment for this frame
    int idx_;                  ///< Residue index in topology
    int C_;                    ///< Coord idx of BB carbon
    int O_;                    ///< Coord idx of BB oxygen
    int N_;                    ///< Coord idx of BB nitrogen
    int H_;                    ///< Coord idx of BB hydrogen
    int CA_;                   ///< Coord idx of BB alpha carbon
    bool isSelected_;          ///< True if calculating SS for this residue.
};
#endif
