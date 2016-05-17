#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
// ----- BOND/ANGLE/DIHEDRAL PARAMETERS ----------------------------------------
/// Hold bond parameters
class BondParmType {
  public:
    BondParmType() : rk_(0), req_(0) {}
    BondParmType(double rk, double req) : rk_(rk), req_(req) {}
    inline double Rk()  const { return rk_;  }
    inline double Req() const { return req_; }
    inline void SetRk(double rk)   { rk_ = rk;   }
    inline void SetReq(double req) { req_ = req; }
    bool operator<(const BondParmType& rhs) const {
      if (rk_ == rhs.rk_) {
        return (req_ < rhs.req_);
      } else return (rk_ < rhs.rk_);
    }
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
    bool operator<(const BondType& rhs) const {
      if (a1_ == rhs.a1_) {
        return (a2_ < rhs.a2_);
      } else return (a1_ < rhs.a1_);
    }
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
    inline void SetTk(double tk)   { tk_ = tk;   }
    inline void SetTeq(double teq) { teq_ = teq; } 
    bool operator<(const AngleParmType& rhs) const {
      if (tk_ == rhs.tk_) {
        return (teq_ < rhs.teq_);
      } else return (tk_ < rhs.tk_);
    }
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
    void SetIdx(int i)     { idx_ = i;    }
    bool operator<(AngleType const& rhs) const {
      if (a1_ == rhs.a1_) {
        if (a2_ == rhs.a2_) {
          return (a3_ < rhs.a3_);
        } else return (a2_ < rhs.a2_);
      } else return (a1_ < rhs.a1_);
    }
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
    DihedralParmType(double k, double p) :
                         pk_(k), pn_(0), phase_(p), scee_(0), scnb_(0) {}
    inline double Pk()    const { return pk_;    }
    inline double& Pk()         { return pk_;    }
    inline double Pn()    const { return pn_;    }
    inline double Phase() const { return phase_; }
    inline double SCEE()  const { return scee_;  }
    inline double SCNB()  const { return scnb_;  }
    void SetPk(double k)        { pk_ = k;       }
    void SetPn(double n)        { pn_ = n;       }
    void SetPhase(double p)     { phase_ = p;    }
    void SetSCEE(double s)      { scee_ = s;     }
    void SetSCNB(double s)      { scnb_ = s;     }
    bool operator<(DihedralParmType const& rhs) const {
      if (pk_ == rhs.pk_) {
        if (pn_ == rhs.pn_) {
          if (phase_ == rhs.phase_) {
            if (scee_ == rhs.scee_) {
              return (scnb_ < rhs.scnb_);
            } else return (scee_ < rhs.scee_);
          } else return (phase_ < rhs.phase_);
        } else return (pn_ < rhs.pn_);
      } else return (pk_ < rhs.pk_);
    }
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
    void SetIdx(int i)        { idx_ = i;     }
    bool operator<(DihedralType const& rhs) const {
      if (a1_ == rhs.a1_) {
        if (a2_ == rhs.a2_) {
          if (a3_ == rhs.a3_) {
            return (a4_ < rhs.a4_);
          } else return (a3_ < rhs.a3_);
        } else return (a2_ < rhs.a2_);
      } else return (a1_ < rhs.a1_);
    }
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
                           asol_(a), bsol_(b), hbcut_(c) {}
    inline double Asol()  const { return asol_;  }
    inline double Bsol()  const { return bsol_;  }
    inline double HBcut() const { return hbcut_; }
    void SetAsol(double a)  { asol_ = a;  }
    void SetBsol(double b)  { bsol_ = b;  }
    void SetHBcut(double h) { hbcut_ = h; }
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
    void SetA(double a) { A_ = a; }
    void SetB(double b) { B_ = b; }
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
//    NonbondParmType(int n, std::vector<int> const& nbi, NonbondArray const& nba,
//                    HB_ParmArray const& hba) :
//                      ntypes_(n), nbindex_(nbi), nbarray_(nba), hbarray_(hba) {}
    inline bool HasNonbond()             const { return ntypes_ > 0; }
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
    /// Set number of types and init nonbond index array.
    void SetNtypes(int n) {
      ntypes_ = n;
      nbindex_.assign(ntypes_ * ntypes_, -1); 
    }
    /// Set number of LJ terms and init LJ array TODO combine with SetNtypes?
    void SetNLJterms(int n)   { nbarray_.assign( n, NonbondType() ); }
    /// Set specified LJ term
    NonbondType& SetLJ(int i) { return nbarray_[i];                  }
    /// Set number of HB terms and init HB array TODO combine with SetNtypes?
    void SetNHBterms(int n)   { hbarray_.assign( n, HB_ParmType() ); }
    /// Set specified HB term
    HB_ParmType& SetHB(int i) { return hbarray_[i];                  }
    /// Set specified nbindex location to given value. FIXME no bounds check
    void SetNbIdx(int idx, int nbidx) { nbindex_[idx] = nbidx; }
    /// Add given LJ term to nonbond array and update nonbond index array.
    void AddLJterm(int ndx, int type1, int type2, NonbondType const& LJ) {
      nbindex_[ntypes_ * type1 + type2] = ndx;
      nbindex_[ntypes_ * type2 + type1] = ndx;
      if (ndx >= (int)nbarray_.size())
        nbarray_.resize(ndx+1);
      nbarray_[ndx] = LJ;
    }
    /// Add given HB term to HB array and update the nonbond index array.
    void AddHBterm(int type1, int type2, HB_ParmType const& HB) {
      int ndx = -((int)hbarray_.size())-1;
      nbindex_[ntypes_ * type1 + type2] = ndx;
      nbindex_[ntypes_ * type2 + type1] = ndx;
      hbarray_.push_back( HB );
    }
    void Clear() { ntypes_ = 0; nbindex_.clear(); nbarray_.clear(); hbarray_.clear(); }
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
    inline bool HasLES()                const { return ntypes_ > 0;      }
    inline int Ntypes()                 const { return ntypes_;          }
    inline int Ncopies()                const { return ncopies_;         }
    std::vector<double> const& FAC()    const { return fac_;             }
    LES_Array           const& Array()  const { return array_;           }
    void SetTypes(int n, std::vector<double> const& f) {
      ntypes_ = n;
      fac_ = f;
    }
    // FIXME: It seems that ncopies is not correctly reported in LES
    //        topology files. Do a manual count until this is fixed.
    void AddLES_Atom(LES_AtomType const& lat) {
      array_.push_back( lat );
      if (array_.back().Copy() > ncopies_ )
        ncopies_ = array_.back().Copy();
    }
    void Clear() { ntypes_ = 0; ncopies_ = 0; array_.clear(); fac_.clear(); }
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
    inline bool HasWaterCap() const { return cutcap_ > 0.0; }
    inline int NatCap()       const { return natcap_; }
    inline double CutCap()    const { return cutcap_; }
    inline double xCap()      const { return xcap_;   }
    inline double yCap()      const { return ycap_;   }
    inline double zCap()      const { return zcap_;   }
    void Clear() { natcap_ = 0; cutcap_ = 0.0; xcap_ = 0.0; ycap_ = 0.0; zcap_ = 0.0; }
    void SetNatcap(int n)    { natcap_ = n; }
    void SetCutCap(double c) { cutcap_ = c; }
    void SetXcap(double x)   { xcap_ = x;   }
    void SetYcap(double y)   { ycap_ = y;   }
    void SetZcap(double z)   { zcap_ = z;   }
  private:
    int natcap_;    ///< last atom before the start of the cap of waters
    double cutcap_; ///< the distance from the center of the cap to the outside
    double xcap_;   ///< X coordinate for the center of the cap
    double ycap_;   ///< Y coordinate for the center of the cap
    double zcap_;   ///< Z coordinate for the center of the cap
};
// ----- CHAMBER PARAMETERS ----------------------------------------------------
/// Hold CMAP grid parameters
class CmapGridType {
  public:
    CmapGridType() : resolution_(0) {}
    CmapGridType(int r, std::vector<double> const& g) :
                     resolution_(r), grid_(g) {}
    inline int Resolution()                  const { return resolution_; }
    inline std::vector<double> const& Grid() const { return grid_;       }
  private:
    int resolution_;           ///< Number of steps along each phi/psi CMAP axis
    std::vector<double> grid_; ///< CMAP grid
};
typedef std::vector<CmapGridType> CmapGridArray;
/// Hold CMAP atom indices and corresponding grid index
class CmapType {
  public:
    CmapType() : a1_(0), a2_(0), a3_(0), a4_(0), a5_(0), idx_(0) {}
    CmapType(int a1, int a2, int a3, int a4, int a5, int i) :
                 a1_(a1), a2_(a2), a3_(a3), a4_(a4), a5_(a5), idx_(i) {}
    inline int A1()     const { return a1_;   }
    inline int A2()     const { return a2_;   }
    inline int A3()     const { return a3_;   }
    inline int A4()     const { return a4_;   }
    inline int A5()     const { return a5_;   }
    inline int Idx()    const { return idx_;  }
  private:
    int a1_;
    int a2_;
    int a3_;
    int a4_;
    int a5_;
    int idx_;
};
typedef std::vector<CmapType> CmapArray;
/// Hold CHAMBER parameters
class ChamberParmType {
  public:
    ChamberParmType() : chmff_verno_(-1) {}
    bool                     HasChamber()   const { return chmff_verno_ > -1; }
    bool                     HasCmap()      const { return !cmapGrid_.empty(); }
    int                      FF_Version()   const { return chmff_verno_;  }
    std::string       const& FF_Type()      const { return chmff_type_;   }
    BondArray         const& UB()           const { return ub_;           }
    BondParmArray     const& UBparm()       const { return ubparm_;       }
    DihedralArray     const& Impropers()    const { return impropers_;    }
    DihedralParmArray const& ImproperParm() const { return improperparm_; }
    NonbondArray      const& LJ14()         const { return lj14_;         }
    CmapGridArray     const& CmapGrid()     const { return cmapGrid_;     }
    CmapArray         const& Cmap()         const { return cmap_;         }
    void SetLJ14(NonbondArray const& nb)          { lj14_ = nb;           }
    NonbondType SetLJ14(int idx)                  { return lj14_[idx];    }
    /// Set expected number of LJ14 terms TODO combine with SetVersion?
    void SetNLJ14terms(int n)                     { lj14_.assign( n, NonbondType() ); }
    void SetVersion(int i, std::string const& s)  { 
      chmff_verno_ = i;
      chmff_type_ = s;
    }
    void ReserveUBterms(unsigned int n)       { ub_.reserve( n );         }
    void AddUBterm(BondType const& bnd)       { ub_.push_back( bnd );     }
    void ResizeUBparm(unsigned int n)         { ubparm_.resize( n );      }
    BondParmType& SetUBparm(unsigned int idx) { return ubparm_[idx];      }
    void SetUB(BondArray const& ub, BondParmArray const& ubp) {
      ub_ = ub;
      ubparm_ = ubp;
    }
    void ReserveImproperTerms(unsigned int n)         { impropers_.reserve( n );     }
    void AddImproperTerm(DihedralType const& dih)     { impropers_.push_back( dih ); }
    void ResizeImproperParm(unsigned int n)           { improperparm_.resize( n );   }
    DihedralParmType& SetImproperParm(unsigned int i) { return improperparm_[i];     }
    void SetImproper(DihedralArray const& im, DihedralParmArray const& imp) {
      impropers_ = im;
      improperparm_ = imp;
    }
    void AddCmapGrid(CmapGridType const& g) { cmapGrid_.push_back(g); }
    void AddCmapTerm(CmapType const& c)     { cmap_.push_back(c);     }
    void Clear() {
      chmff_verno_=-1; chmff_type_.clear(); ub_.clear(); ubparm_.clear();
      impropers_.clear(); improperparm_.clear(); lj14_.clear();
      cmapGrid_.clear(); cmap_.clear();
    } 
  private:
    int chmff_verno_;                ///< CHARMM FF version number
    std::string chmff_type_;         ///< CHARMM FF type 
    BondArray ub_;                   ///< Urey-Bradley terms
    BondParmArray ubparm_;           ///< Urey-Bradley parameters
    DihedralArray impropers_;        ///< Improper terms
    DihedralParmArray improperparm_; ///< Improper parameters
    NonbondArray lj14_;              ///< Lennard-Jones 1-4 parameters
    CmapGridArray cmapGrid_;         ///< Hold CMAP grids
    CmapArray cmap_;                 ///< Hold atom indices and CMAP grid index
};
#endif
