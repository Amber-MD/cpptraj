#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
#include <vector>
#include "NameType.h" // FIXME: Needed?
// ----- BOND/ANGLE/DIHEDRAL PARAMETERS ----------------------------------------
/// Hold bond parameters
class BondParmType {
  public:
    BondParmType() : rk_(0), req_(0) {}
    BondParmType(double rk, double req) : rk_(rk), req_(req) {}
    inline double Rk()  const { return rk_;  }
    inline double Req() const { return req_; }
  private:
    double rk_;
    double req_;
};
typedef std::vector<BondParmType> BondParmArray;
/// Hold bonded atom indices and parameter index
class BondType {
  public:
    BondType() : a1_(0), a2_(0), idx_(0) {}
    BondType(int a1, int a2, int idx) : a1_(a1), a2_(a2), idx_(idx) {}
    inline int A1()  const { return a1_;  }
    inline int A2()  const { return a2_;  }
    inline int Idx() const { return idx_; }
    void SetIdx(int i)     { idx_ = i;    }
  private:
    int a1_;
    int a2_;
    int idx_;
};
typedef std::vector<BondType> BondArray;
/// Hold angle parameters
class AngleParmType {
  public:
    AngleParmType() : tk_(0), teq_(0) {}
    AngleParmType(double tk, double teq) : tk_(tk), teq_(teq) {}
    inline double Tk()  const { return tk_;  }
    inline double Teq() const { return teq_; }
  private:
    double tk_;
    double teq_;
};
typedef std::vector<AngleParmType> AngleParmArray;
/// Hold angle atom indices and parameter index
class AngleType {
  public:
    AngleType() : a1_(0), a2_(0), a3_(0), idx_(0) {}
    AngleType(int a1, int a2, int a3, int idx) :
                  a1_(a1), a2_(a2), a3_(a3), idx_(idx) {}
    inline int A1()  const { return a1_;  }
    inline int A2()  const { return a2_;  }
    inline int A3()  const { return a3_;  }
    inline int Idx() const { return idx_; }
  private:
    int a1_;
    int a2_;
    int a3_;
    int idx_;
};
typedef std::vector<AngleType> AngleArray;
/// Hold dihedral parameters
class DihedralParmType {
  public:
    DihedralParmType() : pk_(0), pn_(0), phase_(0), scee_(0), scnb_(0) {}
    DihedralParmType(double k, double n, double p, double e, double b) :
                         pk_(k), pn_(n), phase_(p), scee_(e), scnb_(b) {}
    inline double Pk()    const { return pk_;    }
    inline double& Pk()         { return pk_;    }
    inline double Pn()    const { return pn_;    }
    inline double Phase() const { return phase_; }
    inline double SCEE()  const { return scee_;  }
    inline double SCNB()  const { return scnb_;  }
    void SetSCEE(double s)      { scee_ = s;     }
    void SetSCNB(double s)      { scnb_ = s;     }
  private:
    double pk_;
    double pn_;
    double phase_;
    double scee_;
    double scnb_;
};
typedef std::vector<DihedralParmType> DihedralParmArray;
/// Hold dihedral atom indices and parameter index
class DihedralType {
  public:
    /// Dihedral type; a3_ < 0 = E, a4_ < 0 = I
    enum Dtype { NORMAL=0, IMPROPER, END, BOTH };
    DihedralType() : a1_(0), a2_(0), a3_(0), a4_(0), type_(NORMAL), idx_(0) {}
    /// For use with Amber-style dihedral array
    DihedralType(int a1, int a2, int a3, int a4, int idx) :
                  a1_(a1), a2_(a2), a3_(a3), a4_(a4), idx_(idx)
    {
      if (a3_ < 0 && a4_ < 0) { a3_ = -a3_; a4_ = -a4_; type_ = BOTH;    }
      else if (a3_ < 0)       { a3_ = -a3;              type_ = END;     }
      else if (a4_ < 0)       { a4_ = -a4;              type_ = IMPROPER;}
      else                                              type_ = NORMAL;
    }
    DihedralType(int a1, int a2, int a3, int a4, Dtype t, int i) :
                 a1_(a1), a2_(a2), a3_(a3), a4_(a4), type_(t), idx_(i) {}
    inline int A1()     const { return a1_;   }
    inline int A2()     const { return a2_;   }
    inline int A3()     const { return a3_;   }
    inline int A4()     const { return a4_;   }
    inline Dtype Type() const { return type_; }
    inline int Idx()    const { return idx_;  }
  private:
    int a1_;
    int a2_;
    int a3_;
    int a4_;
    Dtype type_;
    int idx_;
};
typedef std::vector<DihedralType> DihedralArray;
// ----- NON-BONDED PARAMETERS -------------------------------------------------
/// Hold LJ 10-12 hbond params
class HB_ParmType {
  public:
    HB_ParmType() : asol_(0), bsol_(0), hbcut_(0) {}
    HB_ParmType(double a, double b, double c) :
                           asol_(a), bsol_(b), hbcut_(0) {}
    inline double Asol()  const { return asol_;  }
    inline double Bsol()  const { return bsol_;  }
    inline double HBcut() const { return hbcut_; }
  private:
    double asol_;
    double bsol_;
    double hbcut_;
};
typedef std::vector<HB_ParmType> HB_ParmArray;
/// Hold Lennard-Jones 6-12 interaction A and B parameters
class NonbondType {
  public:
    NonbondType() : A_(0), B_(0) {}
    NonbondType(double rk, double req) : A_(rk), B_(req) {}
    inline double A() const { return A_; }
    inline double B() const { return B_; }
  private:
    double A_;
    double B_;
};
typedef std::vector<NonbondType> NonbondArray;
/// Hold nonbonded interaction parameters
class NonbondParmType {
  public:
    NonbondParmType() : ntypes_(0) {}
    inline int Ntypes()           const { return ntypes_;  }
    std::vector<int> NBindex()    const { return nbindex_; }
    NonbondArray const& NBarray() const { return nbarray_; }
    /// In Amber, index < 0 means HB, otherwise LJ 6-12
    int GetLJindex(int type1, int type2) const {
      return nbindex_[ ntypes_ * type1 + type2 ];
    }
    NonbondType const& NBarray(int i) const { return nbarray_[i]; }
    HB_ParmType const& HBarray(int i) const { return hbarray_[i]; }
  private:
    int ntypes_;               ///< Number of unique atom types
    std::vector<int> nbindex_; ///< Hold indices into Lennard-Jones array nbarray
    NonbondArray nbarray_;     ///< Hold Lennard-Jones 6-12 A and B parameters for all pairs.
    HB_ParmArray hbarray_;     ///< Hold 10-12 Amber HBond params for all pairs.
};
// ----- LES PARAMETERS --------------------------------------------------------
/// Hold LES atom parameters
class LES_AtomType {
  public:
    LES_AtomType() : type_(0), cnum_(0), id_(0) {}
    LES_AtomType(int t, int c, int i) : type_(t), cnum_(c), id_(i) {}
    inline int Type() const { return type_; }
    inline int Copy() const { return cnum_; }
    inline int ID()   const { return id_;   }
  private:
    int type_;
    int cnum_;
    int id_;
};
typedef std::vector<LES_AtomType> LES_Array;
/// Hold LES parameters
class LES_ParmType {
  public:
    LES_ParmType() : ntypes_(0) {}
    LES_ParmType(int na, int nt, std::vector<double> const& fac) : ntypes_(nt), fac_(fac) {
      array_.reserve( na );
    }
    inline int Ntypes()                 const { return ntypes_;          }
    std::vector<double> const& FAC()    const { return fac_;             }
    LES_Array const& Array()            const { return array_;           }
    void AddLES_Atom(LES_AtomType const& lat) { array_.push_back( lat ); }
  private:
    int ntypes_;              ///< Number of LES regions
    LES_Array array_;         ///< LES parameters for each atom
    std::vector<double> fac_; ///< Scaling factor for typeA * typeB
};
// ----- PERTURBATION PARAMETERS -----------------------------------------------
/// Hold perturbed atom parameters
class PertAtom {
  public:
    PertAtom() : atname_(""), atsym_(""), lambda_(0), isPert_(0), atype_(0), charge_(0) {}
  private:
    NameType atname_; ///< atom name at lambda = 1 (igrper)
    NameType atsym_;  ///< atomic symbol at lambda = 1 (ismper)
    double lambda_;   ///< value of lambda for atom (almper)
    int isPert_;      ///< = 1 if atom is being perturbed (iaper)
    int atype_;       ///< atom type at lambda = 1 (iacper)
    double charge_;   ///< atom charge at lambda = 1 (cgper)
};
/// Hold perturbation parameters
class PertParmType {
  public:
    PertParmType() {}
  private:                           // Original: ixper, jxper, (kxper,( lpper,)) icxper{0,1}
    std::vector<int> pbond_;         ///< {AtomIdx1, AtomIdx2, ParmIdx0, ParmIdx1}
    std::vector<int> pangle_;        ///< {AtomIdx1, AtomIdx2, AtomIdx3, ParmIdx0, ParmIdx1}
    std::vector<int> pdih_;          ///< {AtomIdx1, AtomIdx2, AtomIdx3, AtomIdx4, ParmIdx0, ParmIdx1}
    std::vector<NameType> resnames_; ///< residue names at lambda = 1 (labper)
    std::vector<PertAtom> patoms_;   ///< Perturbed atom array
};
#endif
