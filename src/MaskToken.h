#ifndef INC_MASKTOKEN_H
#define INC_MASKTOKEN_H
#include "Atom.h"
#include "Residue.h"
/// Hold information used in mask selection. 
class MaskToken {
  public:
    enum MaskTokenType { 
      OP_NONE=0, ResNum, ResName, AtomNum, AtomName, AtomType, AtomElement, MolNum,
      SelectAll, OP_AND, OP_OR, OP_NEG, OP_DIST
    };
    MaskToken();
    const char *TypeName() const;
    void Print() const;
    int SetToken( MaskTokenType, std::string const& );
    int SetDistance( std::string & );
    void SetOperator(MaskTokenType t) { type_ = t; onStack_ = false; }
    void SetOnStack()                 { onStack_ = true;             }

    inline MaskTokenType Type()   const { return type_;     }
    inline int Res1()             const { return res1_;     }
    inline int Res2()             const { return res2_;     }
    inline const NameType& Name() const { return name_;     }
    inline bool OnStack()         const { return onStack_;  }
    inline bool Within()          const { return d_within_; }
    inline bool ByAtom()          const { return d_atom_;   }
    inline double Distance()      const { return distance_; }
  private:
    static const char* MaskTypeString[];

    int MakeNameType();

    MaskTokenType type_;
    int res1_;
    int res2_;
    NameType name_;
    bool onStack_;
    // Distance criteria
    bool d_within_;
    bool d_atom_;
    double distance_;
};
// =============================================================================
/// Hold an array of MaskTokens. Basis of all Mask classes.
class MaskTokenArray {
  public:
    typedef std::vector<Atom> AtomArrayT;
    typedef std::vector<Residue> ResArrayT;
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
    virtual int SetupMask(AtomArrayT const&, ResArrayT const&, const double*) = 0;
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
    char* ParseMask(AtomArrayT const&, ResArrayT const&, const double*) const;
    static char SelectedChar_;
    static char UnselectedChar_; 
  private:
    typedef std::vector<MaskToken> MTarray;
    typedef MTarray::const_iterator token_iterator;

    static bool IsOperator(char);
    static bool IsOperand(char);
    static int OperatorPriority(char);
    /// Convert mask expression to MaskTokens
    int Tokenize();

    // Mask selection routines.
    int Mask_SelectDistance(const double*, char *, MaskToken const&,
                            AtomArrayT const&, ResArrayT const&) const;
    void Mask_AND(char*, char*, unsigned int) const;
    void Mask_OR(char*, char*, unsigned int) const;
    void Mask_NEG(char *, unsigned int) const;
    void MaskSelectResidues(ResArrayT const&, NameType const&, char*) const;
    void MaskSelectResidues(ResArrayT const&, int, int, char*) const;
    void MaskSelectElements(AtomArrayT const&, NameType const&, char*) const;
    void MaskSelectTypes(AtomArrayT const&, NameType const&, char*) const;
    void MaskSelectAtoms(AtomArrayT const&, NameType const&, char*) const;
    void MaskSelectAtoms(AtomArrayT const&, int, int, char*) const;

    MTarray maskTokens_;
    std::string maskExpression_; ///< String specifying atom mask selection.
    int debug_;
};
#endif
