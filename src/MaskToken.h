#ifndef INC_MASKTOKEN_H
#define INC_MASKTOKEN_H
#include "Atom.h"
#include "Residue.h"
#include "Molecule.h"
#include "SymbolExporting.h"

/// Hold information used in mask selection. 
class MaskToken {
  public:
    enum MaskTokenType { 
      OP_NONE=0,
      ResNum, ResName, ResChain, OresNum,
      AtomNum, AtomName, AtomType, AtomElement,
      MolNum,
      SelectAll, OP_AND, OP_OR, OP_NEG, OP_DIST
    };
    enum DistOpType { BY_ATOM = 0, BY_RES, BY_MOL, BY_RESCENTER };
    MaskToken();
    const char *TypeName() const;
    void Print() const;
    int SetToken( MaskTokenType, std::string const& );
    int SetDistance( std::string const& );
    void SetOperator(MaskTokenType t) { type_ = t; onStack_ = false; }
    void SetOnStack()                 { onStack_ = true;             }

    inline MaskTokenType Type()   const { return type_;      }
    inline int Idx1()             const { return idx1_;      }
    inline int Idx2()             const { return idx2_;      }
    inline const NameType& Name() const { return name_;      }
    inline bool OnStack()         const { return onStack_;   }
    inline bool Within()          const { return d_within_;  }
    inline DistOpType DistOp()    const { return distOp_;    }
    inline double Distance2()     const { return distance2_; }
  private:
    static const char* MaskTypeString[];

    int MakeNameType();

    double distance2_;   ///< Distance cutoff squared
    NameType name_;      ///< Atom name/type/element, residue name, chain ID
    MaskTokenType type_; ///< Mask token type
    DistOpType distOp_;  ///< Distance selection type
    int idx1_;           ///< Begin atom/residue/molecule index
    int idx2_;           ///< End atom/residue/molecule index
    bool onStack_;       ///< True if resulting mask needs to go on stack
    bool d_within_;      ///< True if distance selection is within
};
// =============================================================================
/// Hold an array of MaskTokens. Basis of all Mask classes.
class MaskTokenArray {
  public:
    typedef std::vector<Atom> AtomArrayT;
    typedef std::vector<Residue> ResArrayT;
    typedef std::vector<Molecule> MolArrayT;
    MaskTokenArray();
    // Virtual destructor since this will be inherited
    virtual ~MaskTokenArray() {}
    /// \return original mask expression as std::string
    std::string const& MaskExpression() const { return maskExpression_;     }
    /// \return original mask expression as const char* 
    const char* MaskString()            const { return maskExpression_.c_str(); }
    /// \return true if mask expression has been set.
    bool MaskStringSet()                const { return !maskExpression_.empty(); }
    void SetDebug(int d)                      { debug_ = d;                 }
    /// Set mask expression and tokens
    int SetMaskString(std::string const&);
    /// Set mask expression and tokens
    int SetMaskString(const char*);
    /// Print mask string and number of selected atoms.
    void MaskInfo() const;
    /// Print [<expression>](<# selected>), no newline.
    void BriefMaskInfo() const;
    // -------------------------------------------
    /// Print selected atoms to screen.
    virtual void PrintMaskAtoms(const char*) const = 0;
    /// Select atoms based on current MaskTokens given atom/residue info
    virtual int SetupMask(AtomArrayT const&, ResArrayT const&, MolArrayT const&, const double*) = 0;
    /// Clear all mask information.
    virtual void ResetMask() = 0;
    /// Clear selected atoms only.
    virtual void ClearSelected() = 0;
    /// Invert mask selection.
    virtual void InvertMask() = 0;
    /// \return Number of selected atoms.
    virtual int Nselected() const = 0;
    // -------------------------------------------
    /// \return true if no atoms selected
    bool None() const { return Nselected() == 0; }
  protected:
    void ClearTokens() { maskTokens_.clear(); maskExpression_.clear(); }
    /// \return array of characters with selected atoms marked with SelectedChar_
    char* ParseMask(AtomArrayT const&, ResArrayT const&, MolArrayT const&, const double*) const;
    static CPPTRAJ_EXPORT char SelectedChar_;
    static CPPTRAJ_EXPORT char UnselectedChar_;
  private:
    typedef std::vector<MaskToken> MTarray;
    typedef MTarray::const_iterator token_iterator;

    static bool IsOperator(char);
    static bool IsOperand(char);
    static int OperatorPriority(char);
    /// Convert mask expression to MaskTokens
    int Tokenize();

    // Mask selection routines.
    int SelectDistance(const double*, char *, MaskToken const&,
                       AtomArrayT const&, ResArrayT const&, MolArrayT const&) const;
    void Mask_AND(char*, char*, unsigned int) const;
    void Mask_OR(char*, char*, unsigned int) const;
    void Mask_NEG(char *, unsigned int) const;
    void SelectResName(ResArrayT const&, NameType const&, char*) const;
    void SelectResNum(ResArrayT const&, int, int, char*) const;
    void SelectChainID(ResArrayT const&, NameType const&, char*) const;
    void SelectOriginalResNum(ResArrayT const&, int, int, char*) const;
    void SelectMolNum(MolArrayT const&, int, int, char*) const;
    void SelectElement(AtomArrayT const&, NameType const&, char*) const;
    void SelectAtomType(AtomArrayT const&, NameType const&, char*) const;
    void SelectAtomName(AtomArrayT const&, NameType const&, char*) const;
    void SelectAtomNum(AtomArrayT const&, int, int, char*) const;

    MTarray maskTokens_;
    std::string maskExpression_; ///< String specifying atom mask selection.
    int debug_;
};
#endif
