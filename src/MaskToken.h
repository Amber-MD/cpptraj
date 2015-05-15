#ifndef INC_MASKTOKEN_H
#define INC_MASKTOKEN_H
#include <vector>
#include "Atom.h"
#include "Residue.h"
/// Hold information use in mask selection. 
class MaskToken {
  public:
    enum MaskTokenType { 
      OP_NONE=0, ResNum, ResName, AtomNum, AtomName, AtomType, AtomElement, SelectAll,
      OP_AND, OP_OR, OP_NEG, OP_DIST
    };
    MaskToken();
    const char *TypeName() const;
    void Print() const;
    int SetToken( MaskTokenType, std::string const& );
    int SetDistance( std::string & );
    void SetOperator(MaskTokenType);

    inline MaskTokenType Type()   const { return type_;     }
    inline int Res1()             const { return res1_;     }
    inline int Res2()             const { return res2_;     }
    inline const NameType& Name() const { return name_;     }
    inline bool OnStack()         const { return onStack_;  }
    inline bool Within()          const { return d_within_; }
    inline bool ByAtom()          const { return d_atom_;   }
    inline double Distance()      const { return distance_; }

    void SetOnStack();
  private:
    static const char* MaskTypeString[];

    MaskTokenType type_;
    int res1_;
    int res2_;
    NameType name_;
    bool onStack_;
    // Distance criteria
    bool d_within_;
    bool d_atom_;
    double distance_;

    void MakeNameType();
};

/// Hold an array of MaskTokens. Basis of all Mask classes.
class MaskTokenArray {
    typedef std::vector<MaskToken> MTarray;
    typedef MTarray::const_iterator token_iterator;
  public:
    typedef std::vector<Atom> AtomArrayT;
    typedef std::vector<Residue> ResArrayT;
    MaskTokenArray() : debug_(0) {}
    //inline token_iterator begintoken()  const { return maskTokens_.begin(); }
    //inline token_iterator endtoken()    const { return maskTokens_.end();   }
    std::string const& MaskExpression() const { return maskExpression_;     }
    const char* MaskString()            const { return maskExpression_.c_str(); }
    bool MaskStringSet()                const { return !maskExpression_.empty(); }
    void SetDebug(int d)                      { debug_ = d;                 }
    int SetMaskString(std::string const&);
    int SetMaskString(const char*);
    virtual int SetupMask(AtomArrayT const&, ResArrayT const&, const double*);
  protected:
    char* ParseMask(AtomArrayT const&, ResArrayT const&, const double*) const; 
  private:
    static bool IsOperator(char);
    static bool IsOperand(char);
    static int OperatorPriority(char);
    int Tokenize();


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
    std::string maskExpression_;
    int debug_;
};
#endif
