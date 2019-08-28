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
    void Print();

    /// Class that will hold SS info for each residue
    class SSres;
    /// Secondary structure types
    enum SStype { NONE=0, EXTENDED, BRIDGE, H3_10, ALPHA, HPI, TURN, BEND };
    static const int NSSTYPE_ = 8;   ///< # of secondary structure types.
    static const char DSSP_char_[];  ///< DSSP 1 char names corresponding to SStype
    static const char* DSSP_name_[]; ///< Full secondary structure names corresponding to SStype
    static const char* SSchar_[];    ///< PTRAJ 1 character names corresponding to SStype
    /// Turn types
    enum TurnType { T3 = 0, T4, T5 };
    static const int NTURNTYPE_ = 3;
    /// Beta types
    enum BetaType { B1 = 0, B2, S };
    static const int NBETATYPE_ = 3;
    /// Bridge direction types
    enum BridgeType { NO_BRIDGE = 0, PARALLEL, ANTIPARALLEL };
    static const int NBRIDGETYPE_ = 3;

    static const double DSSP_fac_;  ///< Original DSSP factor for calc. H-bond "energy"
    static const double DSSP_cut_;  ///< Original DSSP H-bond energy cutoff in kcal/mol

    typedef std::vector<SSres> SSarrayType;
    typedef std::vector<std::string> Sarray;

    int OverHbonds(int, ActionFrame&);

    void AssignBridge(int, int, BridgeType);


    int debug_;            ///< Action debug level
    unsigned int Nframes_; ///< Number of frames processed, for total SS normalization
    SSarrayType Residues_; ///< Hold SS data for all residues.
    Sarray SSname_;        ///< Hold full secondary structure names
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
    bool betaDetail_;            ///< If true use para/anti in place of extended/bridge
};

// =============================================================================
class Action_DSSP2::SSres {
  public:
    SSres();
    SSres(SSres const&);
    SSres& operator=(SSres const&);
    int SScount(SStype t) const { return SScount_[t]; }
    int Bcount(BridgeType t) const { return Bcount_[t]; }
    SStype SS()       const { return sstype_; }
    int Num()         const { return num_; }
    char ResChar()    const { return resChar_; }
    bool IsSelected() const { return isSelected_; }
    int C()           const { return C_; }
    int O()           const { return O_; }
    int N()           const { return N_; }
    int H()           const { return H_; }
    int CA()          const { return CA_; }
    int PrevIdx()     const { return prevIdx_; }
    int NextIdx()     const { return nextIdx_; }
    int Bridge1Idx()  const { return bridge1idx_; }
    BridgeType Bridge1Type() const { return b1type_;}
    int Bridge2Idx()  const { return bridge2idx_; }
    BridgeType Bridge2Type() const { return b2type_;}
    DataSet* Dset()   const { return resDataSet_; }

    bool IsMissingAtoms() const { return (C_==-1 || O_==-1 || N_==-1 || H_==-1 || CA_==-1); }
    bool HasCO()          const { return (C_!=-1 && O_!=-1); }
    bool HasNH()          const { return (N_!=-1 && H_!=-1); }
    bool HasCA()          const { return (CA_!=-1); }

    bool HasTurnStart(TurnType) const;
    bool HasBridge() const;
    bool IsBridgedWith(int) const;
//    char StrandChar() const;
    void PrintSSchar() const;

    void SetTurnBegin(TurnType);
    void SetTurn(TurnType);
    void SetTurnEnd(TurnType);
    /// Set a bridge between this res and other res index into Residues_
    void SetBridge(int, BridgeType);

    void AccumulateData(int, bool, bool);

    void SetSS(SStype);
    void SetNum(int i)       { num_ = i; }
    void SetResChar(char c)  { resChar_ = c; }
    void SetSelected(bool b) { isSelected_ = b; }
    void SetC(int i)         { C_ = i; }
    void SetO(int i)         { O_ = i; }
    void SetN(int i)         { N_ = i; }
    void SetH(int i)         { H_ = i; }
    void SetCA(int i)        { CA_ = i; }
    void SetPrevIdx(int i)   { prevIdx_ = i; }
    void SetNextIdx(int i)   { nextIdx_ = i; }
    void SetDset(DataSet* d) { resDataSet_ = d; }
    /// Deselect this residue and reset coordinate indices.
    void Deselect();
    /// Reset hbonds, pattern and SS type assignments
    void Unassign();
    /// \return Relative priority of currently assigned SS type
    int SSpriority() const;
  private:
    /// \return Relative priority of given SS type
    static inline int ssPriority(SStype);

    DataSet* resDataSet_;       ///< DataSet for SS assignment each frame for this res.
    double chirality_;          ///< Dihedral CA[i-1, i, i+1, i+2]
    double bend_;               ///< Angle CA[i-2, i, i+2]
    int SScount_[NSSTYPE_];     ///< Hold count for each SS type
    int Bcount_[NBRIDGETYPE_];  ///< Hold count for Beta types
    SStype sstype_;             ///< SS assignment for this frame
    int num_;                   ///< Residue index in Topology
    int C_;                     ///< Coord idx of BB carbon
    int O_;                     ///< Coord idx of BB oxygen
    int N_;                     ///< Coord idx of BB nitrogen
    int H_;                     ///< Coord idx of BB hydrogen
    int CA_;                    ///< Coord idx of BB alpha carbon
    int prevIdx_;               ///< Index in Residues_ of previous residue
    int nextIdx_;               ///< Index in Residues_ of next residue
    int bridge1idx_;            ///< Index in Residues_ of res this is bridged to
    BridgeType b1type_;         ///< Type of bridge1
    int bridge2idx_;            ///< Index in Residues_ of res this is bridged to
    BridgeType b2type_;         ///< Type of bridge2
    char resChar_;              ///< Single char residue ID
    char turnChar_[NTURNTYPE_]; ///< Character if part of N turn
    bool isSelected_;           ///< True if calculating SS for this residue.
};
#endif
