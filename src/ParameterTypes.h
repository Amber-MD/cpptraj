#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
#include <vector>
#include "NameType.h" // For perturbed atom info
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
    enum Dtype { NORMAL=0, IMPROPER, END, BOTH };
    DihedralType() : a1_(0), a2_(0), a3_(0), a4_(0), type_(NORMAL), idx_(0) {}
    /// For use with Amber-style dihedral array; a3_ < 0 = E, a4_ < 0 = I
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
    NonbondType(double a, double b) : A_(a), B_(b) {}
    inline double A() const { return A_; }
    inline double B() const { return B_; }
  private:
    double A_;
    double B_;
};
typedef std::vector<NonbondType> NonbondArray;
/// Hold nonbonded interaction parameters
/** The nbindex array holds indices into nbarray (>=0) or hbarray (<0).
  * nbarray size should be (ntypes*(ntypes+1))/2 (half matrix).
  */
class NonbondParmType {
  public:
    NonbondParmType() : ntypes_(0) {}
    NonbondParmType(int n, std::vector<int> const& nbi, NonbondArray const& nba,
                    HB_ParmArray const& hba) :
                      ntypes_(n), nbindex_(nbi), nbarray_(nba), hbarray_(hba) {}
    inline int Ntypes()                  const { return ntypes_;     }
    std::vector<int> const& NBindex()    const { return nbindex_;    }
    NonbondArray     const& NBarray()    const { return nbarray_;    }
    HB_ParmArray     const& HBarray()    const { return hbarray_;    }
    NonbondType const& NBarray(int i)    const { return nbarray_[i]; }
    HB_ParmType const& HBarray(int i)    const { return hbarray_[i]; }
    /// In Amber, index < 0 means HB, otherwise LJ 6-12
    int GetLJindex(int type1, int type2) const {
      return nbindex_[ ntypes_ * type1 + type2 ];
    }
  private:
    int ntypes_;               ///< Number of unique atom types
    std::vector<int> nbindex_; ///< Hold indices into arrays nbarray/hbarray for atom type pairs
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
    int type_; ///< LES atom type
    int cnum_; ///< LES copy #
    int id_;   ///< LES region ID
};
typedef std::vector<LES_AtomType> LES_Array;
/// Hold LES parameters
class LES_ParmType {
  public:
    LES_ParmType() : ntypes_(0), ncopies_(0) {}
    LES_ParmType(int na, int nt, std::vector<double> const& fac) : 
                     ntypes_(nt), ncopies_(0), fac_(fac)
    {
      array_.reserve( na );
    }
    inline int Ntypes()                 const { return ntypes_;          }
    inline int Ncopies()                const { return ncopies_;         }
    std::vector<double> const& FAC()    const { return fac_;             }
    LES_Array           const& Array()  const { return array_;           }
    // FIXME: It seems that ncopies is not correctly reported in LES
    //        topology files. Do a manual count until this is fixed.
    void AddLES_Atom(LES_AtomType const& lat) {
      array_.push_back( lat );
      if (array_.back().Copy() > ncopies_ )
        ncopies_ = array_.back().Copy();
    }
  private:
    int ntypes_;              ///< Total number of LES types 
    int ncopies_;             ///< Total number of LES copies.
    LES_Array array_;         ///< LES parameters for each atom
    std::vector<double> fac_; ///< Scaling factor for typeA * typeB
};
// ----- CAP INFO --------------------------------------------------------------
class CapParmType {
  public:
    CapParmType() : natcap_(0), cutcap_(0), xcap_(0), ycap_(0), zcap_(0) {}
    CapParmType(int n, double c, double x, double y, double z) :
                    natcap_(n), cutcap_(c), xcap_(x), ycap_(y), zcap_(z) {}
    inline int NatCap()    const { return natcap_; }
    inline double CutCap() const { return cutcap_; }
    inline double xCap()   const { return xcap_;   }
    inline double yCap()   const { return ycap_;   }
    inline double zCap()   const { return zcap_;   }
  private:
    int natcap_;    ///< last atom before the start of the cap of waters
    double cutcap_; ///< the distance from the center of the cap to the outside
    double xcap_;   ///< X coordinate for the center of the cap
    double ycap_;   ///< Y coordinate for the center of the cap
    double zcap_;   ///< Z coordinate for the center of the cap
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
  private:                    // Original: ixper, jxper, (kxper,( lpper,)) icxper{1,0}
    std::vector<int> pbond_;  ///< {AtomIdx1, AtomIdx2, ParmIdx1, ParmIdx0}
    std::vector<int> pangle_; ///< {AtomIdx1, AtomIdx2, AtomIdx3, ParmIdx1, ParmIdx0}
    std::vector<int> pdih_;   ///< {AtomIdx1, AtomIdx2, AtomIdx3, AtomIdx4, ParmIdx1, ParmIdx0}
    std::vector<NameType> resnames_; ///< residue names at lambda = 0 (labper)
    std::vector<PertAtom> patoms_;   ///< Perturbed atom array
};
#endif
