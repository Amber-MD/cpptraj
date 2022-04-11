#ifndef INC_ACTION_JCOUPLING_H
#define INC_ACTION_JCOUPLING_H
#include <vector>
#include <map>
#include <string>
#include "Action.h"
#include "CharMask.h"
// Class: Action_Jcoupling
/// Calculate j-coupling values based on dihedrals using the Karplus equation.
/*! \author Original Code: J. Maier, Stony Brook 2011.
    \author Adapted by: D. Roe, Rutgers 2011.

    Two types of J-couplings are calculated by this code. If the type is
    0, the form used is that described by Chou et al. JACS (2003) 125 
    p.8959-8966 and the four constants represent A, B, C, and D. If the
    type is 1 (denoted C in $AMBERHOME/Karplus.txt) the form use is that 
    described by Perez et al. JACS (2001) 123 p.7081-7093 and the first 
    three constants represent C0, C1, and C2.
 */
class Action_Jcoupling: public Action {
  public:
    Action_Jcoupling();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Jcoupling(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    /// Load Karplus parameters from a file
    int loadKarplus(std::string const&);

    /// Jcoupling calculation type.
    enum CalcType {CHOU = 0, PEREZ};

    // ----- jcoupleDihedral ---------------------
    /// Hold Karplus parameters for a dihedral composed of 4 atoms
    class jcoupleDihedral {
      public:
        /// CONSTRUCTOR
        jcoupleDihedral() : type_(CHOU)
        {
          offset_[0] = 0; offset_[1] = 0; offset_[2] = 0; offset_[3] = 0;
          C_[0] = 0; C_[1] = 0; C_[2] = 0; C_[3] = 0;
        }
        /// COPY CONSTRUCTOR
        jcoupleDihedral(jcoupleDihedral const& rhs) : type_(rhs.type_)
        {
          for (unsigned int i = 0; i != 4; i++) {
            atomName_[i] = rhs.atomName_[i];
            offset_[i] = rhs.offset_[i];
            C_[i] = rhs.C_[i];
          }
          type_ = rhs.type_;
        }
        /// ASSIGNMENT
        jcoupleDihedral& operator=(jcoupleDihedral const& rhs) {
          if (this == &rhs) return *this;
          for (unsigned int i = 0; i != 4; i++) {
            atomName_[i] = rhs.atomName_[i];
            offset_[i] = rhs.offset_[i];
            C_[i] = rhs.C_[i];
          }
          type_ = rhs.type_;
          return *this;
        }
        /// Set calculation type
        void SetType(CalcType t) { type_ = t; }
        /// Set atom res offset
        void SetOffset(int at, int off) { offset_[at] = off; }
        /// Set atom name
        void SetName(int at, const char* n) { atomName_[at] = NameType(n); }
        /// \return pointer to constants array, for reading in
        double* Cptr() { return C_; }
        /// \return reference to specified constant
        double& C(int at) { return C_[at]; }

        /// \return specified atom name
        NameType const& AtomName(int at) const { return atomName_[at]; }
        /// \return calculation type
        CalcType Type()                  const { return type_; }
        /// \return specified atom res offset
        int Offset(int at)               const { return offset_[at]; }
        /// \return specified constant
        double Constant(int i)           const { return C_[i]; }
        /// \return pointer to constants array, for assigning to jcouplingInfo
        const double* Carray()           const { return C_; }
      private:
        NameType atomName_[4]; ///< Name of each atom involved in the dihedral
        int offset_[4];        ///< Residue offset for each atom
        double C_[4];          ///< Constants for each atom
        CalcType type_;        ///< Calculation type (0=Chou, 1=Perez)
    };
    // -------------------------------------------

    /// Hold all dihedral constants for a given residue
    typedef std::vector<jcoupleDihedral> JcoupleDihedralArray;
    /// Pair a residue name with dihedral constants
    typedef std::pair<NameType, JcoupleDihedralArray> JcoupleDihedralPair;
    /// Map residue names to Karplus constants
    typedef std::map<NameType, JcoupleDihedralArray> JcoupleDihedralMap;
    /// Hold constants for all residues
    JcoupleDihedralMap JcoupleData_;
    /// Hold info for single j-coupling calculation TODO deprecate?
    struct jcouplingInfo {
      int residue; ///< Residue number
      int atom[4]; ///< Atom #s of the dihedral
      const double *C;  ///< Pointer to C in associated karplusConstant structure
      CalcType type;    ///< Calculation type (0=Chou, 1=Perez)
      DataSet* data_; ///< Hold Jcoupling vs frame
    };
    /// Hold info for all j-coupling calcs
    std::vector<jcouplingInfo> JcouplingInfo_;

    CharMask Mask1_;         ///< Mask to search for dihedrals in.
    int debug_;              ///< Debug level.
    int Nconstants_;
    Topology* CurrentParm_;
    CpptrajFile* outputfile_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    std::string setname_;
    int setcount_;
};
#endif  
