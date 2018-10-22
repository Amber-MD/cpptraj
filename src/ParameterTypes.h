#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
#include <vector>
#include <cmath> // pow, sqrt, fabs
#include "NameType.h"
#include "Constants.h" // SMALL
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
    DihedralParmType(double k, double n, double p) :
                         pk_(k), pn_(n), phase_(p), scee_(0), scnb_(0) {}
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
/** Dihedrals can be marked normal (A1-A2-A3-A4), end (meaning 1-4 calc should
  * be skipped to avoid overcounting, e.g. for dihedrals with multiple 
  * multiplicities or certain ring dihedrals), improper, or both end and improper.
  */
class DihedralType {
  public:
    enum Dtype { NORMAL=0, IMPROPER, END, BOTH };
    /// Set skip 1-4 (end) and improper status
    inline void SetFromType(Dtype t) {
      switch (t) {
        case NORMAL   : skip14_ = false; improper_ = false; break;
        case IMPROPER : skip14_ = false; improper_ = true; break;
        case END      : skip14_ = true;  improper_ = false; break;
        case BOTH     : skip14_ = true;  improper_ = true; break;
      }
    }
    /// Default constructor
    DihedralType() : a1_(0), a2_(0), a3_(0), a4_(0), idx_(0), skip14_(false), improper_(false) {}
    /// For use with Amber-style dihedral array; a3_ < 0 = E, a4_ < 0 = I
    DihedralType(int a1, int a2, int a3, int a4, int idx) :
                  a1_(a1), a2_(a2), a3_(a3), a4_(a4), idx_(idx)
    {
      if (a3_ < 0 && a4_ < 0) { a3_ = -a3_; a4_ = -a4_; skip14_ = true;  improper_ = true;  }
      else if (a3_ < 0)       { a3_ = -a3;              skip14_ = true;  improper_ = false; }
      else if (a4_ < 0)       { a4_ = -a4;              skip14_ = false; improper_ = true;  }
      else                    {                         skip14_ = false; improper_ = false; }
    }
    /// Takes type, no index
    DihedralType(int a1, int a2, int a3, int a4, Dtype t) :
                 a1_(a1), a2_(a2), a3_(a3), a4_(a4), idx_(-1)
    {
      SetFromType(t);
    }
    /// Takes type and index
    DihedralType(int a1, int a2, int a3, int a4, Dtype t, int i) :
                 a1_(a1), a2_(a2), a3_(a3), a4_(a4), idx_(i)
    {
      SetFromType(t);
    }

    inline int A1()     const { return a1_;    }
    inline int A2()     const { return a2_;    }
    inline int A3()     const { return a3_;    }
    inline int A4()     const { return a4_;    }
    inline int Idx()    const { return idx_;   }  
    void SetIdx(int i)        { idx_ = i;      }
    void SetSkip14(bool b)    { skip14_ = b;   }
    void SetImproper(bool b)  { improper_ = b; }
    /// \return type based on skip 1-4 (end) and improper status
    inline Dtype Type() const {
      if (skip14_ && improper_) return BOTH;
      else if (skip14_)         return END;
      else if (improper_)       return IMPROPER;
      else                      return NORMAL;
    }
    /// Sort based on atom indices
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
    int idx_;
    bool skip14_; ///< If true the 1-4 interaction for this dihedral should be skipped.
    bool improper_; ///< If true this is an improper dihedral.
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
    double Radius() const {
      if (B_ > 0.0)
        return (0.5 * pow(2.0 * A_ / B_, (1.0/6.0)));
      else
        return 0.0;
    }
    double Depth() const {
      if (A_ > 0.0)
        return ( (B_ * B_) / (4.0 * A_) );
      else
        return 0.0;
    }
    /// \return True if A and B match
    bool operator==(NonbondType const& rhs) const {
      return ( (fabs(A_ - rhs.A_) < Constants::SMALL) &&
               (fabs(B_ - rhs.B_) < Constants::SMALL) );
    }
    /// \return True if A and B do not match
    bool operator!=(NonbondType const& rhs) const {
      return ( (fabs(A_ - rhs.A_) > Constants::SMALL) ||
               (fabs(B_ - rhs.B_) > Constants::SMALL) );
    }
  private:
    double A_;
    double B_;
};
typedef std::vector<NonbondType> NonbondArray;
/// Hold Lennard-Jones radius and well-depth
class LJparmType {
  public:
    LJparmType() : radius_(0.0), depth_(0.0) {}
    LJparmType(double r, double d) : radius_(r), depth_(d) {}
    double Radius() const { return radius_; }
    double Depth()  const { return depth_;  }
    void SetRadius(double r) { radius_ = r; }
    void SetDepth(double d)  { depth_ = d;  }
    /// \return True if radius and well depth match
    bool operator==(LJparmType const& rhs) const {
      return ( (fabs(radius_ - rhs.radius_) < Constants::SMALL) &&
               (fabs(depth_  - rhs.depth_ ) < Constants::SMALL) );
    }
    /// \return true if radius and well depth are less in that order
    bool operator<(LJparmType const& rhs) const {
      if (radius_ == rhs.radius_)
        return (depth_ < rhs.depth_);
      else
        return (radius_ < rhs.radius_);
    }
    /// Combine these LJ params with another using Lorentz-Berthelot rules.
    NonbondType Combine_LB(LJparmType const& rhs) const {
      double dR = radius_ + rhs.radius_;
      double dE = sqrt( depth_ * rhs.depth_ );
      double dR2 = dR * dR;
      double dR6 = dR2 * dR2 * dR2;
      double dER6 = dE * dR6;
      return NonbondType( dER6*dR6, 2.0*dER6 );
    }
  private:
    double radius_;
    double depth_;
};
typedef std::vector<LJparmType> LJparmArray;
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
    /// Set number of types, init NB index array, init LJ array.
    void SetupLJforNtypes(int n) { SetNtypes(n); nbarray_.assign((n*(n+1))/2, NonbondType()); }
    /// Set specified LJ term
    NonbondType& SetLJ(int i) { return nbarray_[i];                  }
    /// Set number of HB terms and init HB array TODO combine with SetNtypes?
    void SetNHBterms(int n)   { hbarray_.assign( n, HB_ParmType() ); }
    /// Set specified HB term
    HB_ParmType& SetHB(int i) { return hbarray_[i];                  }
    /// Set specified nbindex location to given value. FIXME no bounds check
    void SetNbIdx(int idx, int nbidx) { nbindex_[idx] = nbidx; }
    /// Add given LJ term to nonbond array and update nonbond index array.
    /** Certain routines in sander (like the 1-4 calcs) do NOT use the 
      * nonbond index array; instead they expect the nonbond arrays to be
      * indexed like '(ibig*(ibig-1)/2+isml)', where ibig is the larger atom
      * type index.
      */
    void AddLJterm(int type1, int type2, NonbondType const& LJ) {
      int ibig, isml;
      if (type1 > type2) {
        ibig = type1 + 1;
        isml = type2 + 1;
      } else {
        ibig = type2 + 1;
        isml = type1 + 1;
      }
      int ndx = (ibig*(ibig-1)/2+isml)-1;
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
    void SetType(int t) { type_ = t; }
    void SetCopy(int c) { cnum_ = c; }
    void SetID(int i)   { id_ = i;   }
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
    /// Prepare LES_ParmType to receive data based on given # atoms and # LES types.
    void Allocate(int natomsIn, int ntypesIn) {
      ntypes_ = ntypesIn;
      ncopies_ = 0;
      array_.clear();
      array_.resize( natomsIn );
      fac_.clear();
      fac_.resize( ntypes_ * ntypes_ );
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
    /// Set LES fac at given position.
    void SetFAC(int idx, double f)  { fac_[idx] = f; }
    /// Set given LES atom type
    void SetType(int idx, int type) { array_[idx].SetType( type ); }
    /// Set given LES atom copy number. FIXME see AddLES_Atom above.
    void SetCopy(int idx, int cnum) {
      array_[idx].SetCopy( cnum );
      if (cnum > ncopies_) ncopies_ = cnum;
    }
    /// Set given LES atom ID
    void SetID(int idx, int id)     { array_[idx].SetID( id ); }
    /// Clear all data
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
    CmapGridType(int r) : resolution_(r), grid_(r*r, 0.0) {}
    inline int Resolution()                  const { return resolution_;       }
    inline std::vector<double> const& Grid() const { return grid_;             }
    inline int Size()                        const { return (int)grid_.size(); }
    void SetGridPt(int idx, double d)              { grid_[idx] = d;           }
  private:
    int resolution_;           ///< Number of steps along each phi/psi CMAP axis
    std::vector<double> grid_; ///< CMAP grid (size is resolution_*resolution_)
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
    typedef std::vector<std::string> Sarray;
  public:
    ChamberParmType() {}
    bool                     HasCmap()      const { return !cmapGrid_.empty(); }
    Sarray            const& Description()  const { return chmff_desc_;   }
    BondArray         const& UB()           const { return ub_;           }
    BondParmArray     const& UBparm()       const { return ubparm_;       }
    DihedralArray     const& Impropers()    const { return impropers_;    }
    DihedralParmArray const& ImproperParm() const { return improperparm_; }
    NonbondArray      const& LJ14()         const { return lj14_;         }
    CmapGridArray     const& CmapGrid()     const { return cmapGrid_;     }
    CmapArray         const& Cmap()         const { return cmap_;         }
    NonbondType& SetLJ14(int idx)                 { return lj14_[idx];    }
    CmapGridType& SetCmapGrid(int idx)            { return cmapGrid_[idx];}
    /// Set expected number of LJ14 terms TODO combine with SetVersion?
    void SetNLJ14terms(int n)                 { lj14_.assign( n, NonbondType() ); }
    void AddDescription(std::string const& s) { chmff_desc_.push_back( s );       }
    void SetDescription(Sarray const& s)      { chmff_desc_ = s;                  }
    void ReserveUBterms(unsigned int n)       { ub_.reserve( n );                 }
    void AddUBterm(BondType const& bnd)       { ub_.push_back( bnd );             }
    BondArray& SetUB()                        { return ub_;                       }
    void ResizeUBparm(unsigned int n)         { ubparm_.resize( n );              }
    BondParmType& SetUBparm(unsigned int idx) { return ubparm_[idx];              }
    BondParmArray& SetUBparm()                { return ubparm_;                   }
    void ReserveImproperTerms(unsigned int n)         { impropers_.reserve( n );     }
    void AddImproperTerm(DihedralType const& dih)     { impropers_.push_back( dih ); }
    DihedralArray& SetImpropers()                     { return impropers_;           }
    void ResizeImproperParm(unsigned int n)           { improperparm_.resize( n );   }
    DihedralParmType& SetImproperParm(unsigned int i) { return improperparm_[i];     }
    DihedralParmArray& SetImproperParm()              { return improperparm_;        }
    void AddCmapGrid(CmapGridType const& g) { cmapGrid_.push_back(g); }
    void AddCmapTerm(CmapType const& c)     { cmap_.push_back(c);     }
    void Clear() {
      chmff_desc_.clear(); ub_.clear(); ubparm_.clear();
      impropers_.clear(); improperparm_.clear(); lj14_.clear();
      cmapGrid_.clear(); cmap_.clear();
    }
    /// \return true if any CHARMM parameters are set (based on indices).
    bool HasChamber() const {
      if (!ub_.empty()) return true;
      if (!impropers_.empty()) return true;
      if (!lj14_.empty()) return true;
      if (!cmap_.empty()) return true;
      return false;
    }
  private:
    Sarray chmff_desc_;              ///< CHARMM FF descriptions
    BondArray ub_;                   ///< Urey-Bradley terms
    BondParmArray ubparm_;           ///< Urey-Bradley parameters
    DihedralArray impropers_;        ///< Improper terms
    DihedralParmArray improperparm_; ///< Improper parameters
    NonbondArray lj14_;              ///< Lennard-Jones 1-4 parameters
    CmapArray cmap_;                 ///< Hold atom indices and CMAP grid index
    CmapGridArray cmapGrid_;         ///< Hold CMAP grids
};
#endif
