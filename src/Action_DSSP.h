#ifndef INC_ACTION_DSSP_H
#define INC_ACTION_DSSP_H
#include <set>
#include <vector>
#include <string>
#include "Action.h"
#include "NameType.h"
#include "CharMask.h"
#include "Timer.h"
class DataSet;
/// Do protein secondary structure assignment. 
/** Based on protein secondary structure definitions given in:
  *   Kabsch, W.; Sander, C.; \"Dictionary of Protein Secondary Structure:
  *   Pattern Recognition of Hydrogen-Bonded and Geometrical Features.
  *   Biopolymers (1983), V.22, pp.2577-2637.
  * This is the second major version of this algorithm in CPPTRAJ; it
  * has been completely rewritten by DRR. It correctly implements 
  * detection of Extended and Bridge regions, including beta bulges.
  * The various secondary structure types that can be assigned (in order
  * of decreasing priority) are:
  *   H - Alpha helix (4)
  *   B - Single bridge (ladder of length 1)
  *   E - Extended beta (ladder length > 1)
  *   G - 3-10 helix (3)
  *   I - Pi helix (5)
  *   T - Turn (3, 4, 5)
  *   S - Bend (Angle CA i-2, i, i+2 > 70 deg.)
  */
class Action_DSSP : public Action {
  public:
    Action_DSSP();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_DSSP(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
#   endif
    void Print();

    /// Class that will hold SS info for each residue
    class SSres;
    /// Secondary structure types. For betaDetail, EXTENDED=PARALLEL and BRIDGE=ANTIPARALLEL
    //            0       1         2       3      4      5    6     7
    enum SStype { NONE=0, EXTENDED, BRIDGE, H3_10, ALPHA, HPI, TURN, BEND };
    static const int NSSTYPE_ = 8;   ///< # of secondary structure types.
    static const char DSSP_char_[];  ///< DSSP 1 char names corresponding to SStype
    static const char* DSSP_name_[]; ///< Full secondary structure names corresponding to SStype
    static const char* SSchar_[];    ///< PTRAJ 1 character names corresponding to SStype
    /// Turn types
    enum TurnType { T3 = 0, T4, T5 };
    static const int NTURNTYPE_ = 3;
    /// Bridge direction types
    enum BridgeType { NO_BRIDGE = 0, PARALLEL, ANTIPARALLEL };
    static const int NBRIDGETYPE_ = 3;

    /// Class for holding bridge info
    class Bridge {
      public:
        Bridge(int i, BridgeType b) : idx_(i), btype_(b) {}
        int Idx()          const { return idx_;   }
        BridgeType Btype() const { return btype_; }
      private:
        int idx_;          ///< Residue to which the bridge is formed.
        BridgeType btype_; ///< The type of bridge being formed.
    };
    typedef std::vector<Bridge> BridgeArray;

    static const double DSSP_fac_;  ///< Original DSSP factor for calc. H-bond "energy"
    static const double DSSP_cut_;  ///< Original DSSP H-bond energy cutoff in kcal/mol

    typedef std::vector<SSres> SSarrayType;
    typedef std::vector<std::string> Sarray;
    typedef std::pair<int,int> HbondPairType;
    typedef std::set<HbondPairType> HbondMapType;

    void CheckBulge(int, int, int);

    int OverHbonds(int, ActionFrame&);

    void AssignBridge(int, int, BridgeType);

    /// Map resi (CO) to resj (NH) hydrogen bonds
#   ifdef _OPENMP
    std::vector<HbondMapType> CO_NH_bondsArray_;
#   else
    HbondMapType CO_NH_bonds_;
#   endif
    int debug_;            ///< Action debug level
    unsigned int Nframes_; ///< Number of frames processed, for total SS normalization
    SSarrayType Residues_; ///< Hold SS data for all residues.
    Sarray SSname_;        ///< Hold full secondary structure names
    NameType BB_N_;        ///< Protein N atom name ('N')
    NameType BB_H_;        ///< Protein N-H atom name ('H')
    NameType BB_C_;        ///< Protein C atom name ('C')
    NameType BB_O_;        ///< Protein C-O atom name ('O')
    NameType BB_CA_;       ///< Protein alpha C name ('CA')
    NameType SG_;          ///< Protein cysteine sulfur name for detecting disulfides
    CharMask Mask_;        ///< Mask used to determine selected residues.
    DataFile* outfile_;          ///< Output Data file
    DataFile* dsspFile_;         ///< Sum output file
    CpptrajFile* assignout_;     ///< Assignment output file.
    std::string dsetname_;       ///< DSSP data set name
    DataSet* totalDS_[NSSTYPE_]; ///< Hold total SS each frame for each SS type
    ActionInit Init_;            ///< Hold pointers to master DSL/DFL
    bool printString_;           ///< If true print 1 char per residue indicating ss type
    bool betaDetail_;            ///< If true use para/anti in place of extended/bridge
    Timer t_total_;
    Timer t_calchb_;
    Timer t_assign_;
};

// =============================================================================
class Action_DSSP::SSres {
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
    BridgeArray const& Bridges() const { return bridges_; }
    DataSet* Dset()   const { return resDataSet_; }

    bool IsMissingAtoms() const { return (C_==-1 || O_==-1 || N_==-1 || H_==-1 || CA_==-1); }
    bool HasCO()          const { return (C_!=-1 && O_!=-1); }
    bool HasNH()          const { return (N_!=-1 && H_!=-1); }
    bool HasCA()          const { return (CA_!=-1); }

    bool HasTurnStart(TurnType) const;
    /// \return Dominant bridge type
    BridgeType DominantBridgeType() const;
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
#   ifdef MPI
    void SyncToMaster(Parallel::Comm const&);
#   endif
  private:
    /// \return Relative priority of given SS type
    static inline int ssPriority(SStype);

    DataSet* resDataSet_;       ///< DataSet for SS assignment each frame for this res.
    double chirality_;          ///< Dihedral CA[i-1, i, i+1, i+2]
    double bend_;               ///< Angle CA[i-2, i, i+2]
    int SScount_[NSSTYPE_];     ///< Hold total count for each SS type over all frames
    int Bcount_[NBRIDGETYPE_];  ///< Hold total count for Beta types over all frames
    SStype sstype_;             ///< SS assignment for this frame
    int num_;                   ///< Residue index in Topology
    int C_;                     ///< Coord idx of BB carbon
    int O_;                     ///< Coord idx of BB oxygen
    int N_;                     ///< Coord idx of BB nitrogen
    int H_;                     ///< Coord idx of BB hydrogen
    int CA_;                    ///< Coord idx of BB alpha carbon
    int prevIdx_;               ///< Index in Residues_ of previous residue
    int nextIdx_;               ///< Index in Residues_ of next residue
    BridgeArray bridges_;       ///< Indices and types of bridges to this residue
    char resChar_;              ///< Single char residue ID
    char turnChar_[NTURNTYPE_]; ///< Character if part of N turn
    bool isSelected_;           ///< True if calculating SS for this residue.
};
#endif
