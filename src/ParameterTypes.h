#ifndef INC_PARAMETERTYPES_H
#define INC_PARAMETERTYPES_H
#include <cstddef> // size_t
#include <cmath> // pow, sqrt, fabs
#include <vector>
#include <string>
#include "Constants.h" // SMALL
// ----- Floating point comparison routines ------------------------------------
/// Floating point equals
static inline bool FEQ(double v1, double v2) {
  double delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta < Constants::SMALL);
}
/// Floating point not equals.
static inline bool FNE(double v1, double v2) {
  double delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta > Constants::SMALL);
}

/// Calculate relative error of A1 from A0 TODO check A0?
static inline double RELERR(double A0, double A1) {
  double delta = A0 - A1;
  if (delta < 0.0) delta = -delta;
  return (delta / A0);
}

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
    bool operator==(const BondParmType& rhs) const {
      return ( FEQ(rk_,  rhs.rk_ ) &&
               FEQ(req_, rhs.req_) );
    }
    bool operator!=(const BondParmType& rhs) const {
      return ( FNE(rk_,  rhs.rk_ ) ||
               FNE(req_, rhs.req_) );
    }
    bool operator<(const BondParmType& rhs) const {
      if (*this != rhs) {
        if (FEQ(rk_, rhs.rk_)) {
          return (req_ < rhs.req_);
        } else return (rk_ < rhs.rk_);
      } else
        return false;
    }
  private:
    double rk_;
    double req_;
};
/// Hold Array of bond parameters
class BondParmArray {
    typedef std::vector<BondParmType> BPArray;
  public:
    BondParmArray() {}
    /// Resize the bond parm array
    void resize(unsigned int n) { bondparm_.resize( n ); }
    /// Clear the bond parm array
    void clear() { bondparm_.clear(); }
    /// \return reference to specified bond parameter
    BondParmType& operator[](unsigned int idx) { return bondparm_[idx]; }
    /// Add bond parameter
    void push_back( BondParmType const& bp ) { bondparm_.push_back( bp ); }

    /// \return const reference to specified bond parameter
    BondParmType const& operator[](unsigned int idx) const { return bondparm_[idx]; }
    /// \return true if no bond parameters
    bool empty() const { return bondparm_.empty(); }
    /// \return number of bond parameters
    size_t size() const { return bondparm_.size(); }
    /// Const iterator
    typedef BPArray::const_iterator const_iterator;
    /// \return const iterator to beginning
    const_iterator begin() const { return bondparm_.begin(); }
    /// \return const iterator to end
    const_iterator end() const { return bondparm_.end(); }
    /// \return Underlying array
    std::vector<BondParmType> const& Array() const { return bondparm_; }
  private:
    BPArray bondparm_;
};
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
    bool operator==(const BondType& rhs) const {
      if ( (a1_ == rhs.a1_ && a2_ == rhs.a2_) ||
           (a2_ == rhs.a1_ && a1_ == rhs.a2_) ) return true;
      return false;
    }
  private:
    int a1_;
    int a2_;
    int idx_;
};
/// Hold array of bonds
class BondArray {
    typedef std::vector<BondType> BArray;
  public:
    /// CONSTRUCTOR
    BondArray() {}

    /// iterator
    typedef BArray::iterator iterator;
    /// begin
    iterator begin() { return bonds_.begin(); }
    /// end
    iterator end()   { return bonds_.end();   }
    /// const iterator
    typedef BArray::const_iterator const_iterator;
    /// const begin
    const_iterator begin() const { return bonds_.begin(); }
    /// const end
    const_iterator end()   const { return bonds_.end();   }

    /// Reserve space for # of bonds
    void reserve(size_t n) { bonds_.reserve(n); }
    /// Add bond
    void push_back(BondType const& b) { bonds_.push_back(b); }
    /// Clear bonds
    void clear() { bonds_.clear(); }
    /// Erase given bond from array
    void erase( iterator bnd ) { bonds_.erase( bnd ); }

    /// \return true if no bonds
    bool empty()  const { return bonds_.empty(); }
    /// \return number of bonds
    size_t size() const { return bonds_.size(); }
    /// \return specified bond
    BondType const& operator[](size_t idx) const { return bonds_[idx]; }
    /// \return underlying array
    std::vector<BondType> const& Array() const { return bonds_; }
  private:
    BArray bonds_;
};
/// Hold angle parameters
class AngleParmType {
  public:
    AngleParmType() : tk_(0), teq_(0) {}
    AngleParmType(double tk, double teq) : tk_(tk), teq_(teq) {}
    inline double Tk()  const { return tk_;  }
    inline double Teq() const { return teq_; }
    inline void SetTk(double tk)   { tk_ = tk;   }
    inline void SetTeq(double teq) { teq_ = teq; }
    bool operator==(const AngleParmType& rhs) const {
      return ( FEQ(tk_,  rhs.tk_ ) &&
               FEQ(teq_, rhs.teq_) );
    }
    bool operator!=(const AngleParmType& rhs) const {
      return ( FNE(tk_,  rhs.tk_ ) ||
               FNE(teq_, rhs.teq_) );
    }
    bool operator<(const AngleParmType& rhs) const {
      if (*this != rhs) {
        if (FEQ(tk_, rhs.tk_)) {
          return (teq_ < rhs.teq_);
        } else return (tk_ < rhs.tk_);
      } else
        return false;
    }
  private:
    double tk_;  ///< Angle force constant in kcal/mol*rad^2
    double teq_; ///< Angle equilibirum value in rad
};
/// Hold Array of angle parameters
class AngleParmArray {
    typedef std::vector<AngleParmType> APArray;
  public:
    AngleParmArray() {}
    /// Resize the angle parm array
    void resize(unsigned int n) { angleparm_.resize( n ); }
    /// Clear the angle parm array
    void clear() { angleparm_.clear(); }
    /// \return reference to specified angle parameter
    AngleParmType& operator[](unsigned int idx) { return angleparm_[idx]; }
    /// Add angle parameter
    void push_back( AngleParmType const& bp ) { angleparm_.push_back( bp ); }

    /// \return const reference to specified angle parameter
    AngleParmType const& operator[](unsigned int idx) const { return angleparm_[idx]; }
    /// \return true if no angle parameters
    bool empty() const { return angleparm_.empty(); }
    /// \return number of angle parameters
    size_t size() const { return angleparm_.size(); }
    /// Const iterator
    typedef APArray::const_iterator const_iterator;
    /// \return const iterator to beginning
    const_iterator begin() const { return angleparm_.begin(); }
    /// \return const iterator to end
    const_iterator end() const { return angleparm_.end(); }
    /// \return Underlying array
    std::vector<AngleParmType> const& Array() const { return angleparm_; }
  private:
    APArray angleparm_;
};
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
/// Hold array of angle parameters
class AngleArray {
    typedef std::vector<AngleType> AArray;
  public:
    /// CONSTRUCTOR
    AngleArray() {}

    /// iterator
    typedef AArray::iterator iterator;
    /// begin
    iterator begin() { return angles_.begin(); }
    /// end
    iterator end()   { return angles_.end();   }
    /// const iterator
    typedef AArray::const_iterator const_iterator;
    /// const begin
    const_iterator begin() const { return angles_.begin(); }
    /// const end
    const_iterator end()   const { return angles_.end();   }

    /// Reserve space for # of angles
    void reserve(size_t n) { angles_.reserve(n); }
    /// Add angle
    void push_back(AngleType const& b) { angles_.push_back(b); }
    /// Clear angles
    void clear() { angles_.clear(); }
    /// Erase given angle from array
    void erase( iterator bnd ) { angles_.erase( bnd ); }

    /// \return true if no angles
    bool empty()  const { return angles_.empty(); }
    /// \return number of angles
    size_t size() const { return angles_.size(); }
    /// \return specified angle
    AngleType const& operator[](size_t idx) const { return angles_[idx]; }
    /// \return underlying array
    std::vector<AngleType> const& Array() const { return angles_; }
    /// \return last angle added
    AngleType& back() { return angles_.back(); }
  private:
    AArray angles_;
};
/// Hold dihedral parameters
/** Note that for the '==' and '<' operators, direct comparisons are used
  * instead of the FEQ() function in order to match how LEaP would compare
  * dihedrals.
  */
class DihedralParmType {
  public:
    /// CONSTRUCTOR
    DihedralParmType() : pk_(0), pn_(0), phase_(0), scee_(0), scnb_(0) {}
    /// CONSTRUCTOR - PK, PN, Phase, SCEE, SCNB (Amber parameter file proper)
    DihedralParmType(double k, double n, double p, double e, double b) :
                         pk_(k), pn_(n), phase_(p), scee_(e), scnb_(b) {}
    /// CONSTRUCTOR - PK, PN, Phase (Amber parameter file improper)
    DihedralParmType(double k, double n, double p) :
                         pk_(k), pn_(n), phase_(p), scee_(0), scnb_(0) {}
    /// \return Dihedral force constant
    inline double Pk()    const { return pk_;    }
    /// \return Dihedral periodicity
    inline double Pn()    const { return pn_;    }
    /// \return Dihedral phase
    inline double Phase() const { return phase_; }
    /// \return 1-4 electrostatics scaling constant
    inline double SCEE()  const { return scee_;  }
    /// \return 1-4 vdW scaling constant
    inline double SCNB()  const { return scnb_;  }
    /// Set dihedral force constant
    void SetPk(double k)        { pk_ = k;       }
    /// Set dihedral periodicity
    void SetPn(double n)        { pn_ = n;       }
    /// Set dihedral phase
    void SetPhase(double p)     { phase_ = p;    }
    /// Set 1-4 electrostatics scaling constant
    void SetSCEE(double s)      { scee_ = s;     }
    /// Set 1-4 vdW scaling constant
    void SetSCNB(double s)      { scnb_ = s;     }
    /// \return True if two dihedral parameters are equal
    bool operator==(DihedralParmType const& rhs) const {
      return ( pn_    == rhs.pn_    &&
               pk_    == rhs.pk_    &&
               phase_ == rhs.phase_ &&
               scee_  == rhs.scee_  &&
               scnb_  == rhs.scnb_ );
    }
    /// \return True if this dihedral parameter should come before the given one
    bool operator<(DihedralParmType const& rhs) const {
      if ( pn_ == rhs.pn_ ) {
        if ( pk_ == rhs.pk_ ) {
          if ( phase_ == rhs.phase_ ) {
            if ( scee_ == rhs.scee_ ) {
              return ( scnb_ < rhs.scnb_ );
            } else return (scee_ < rhs.scee_);
          } else return (phase_ < rhs.phase_);
        } else return (pk_ < rhs.pk_);
      } else return (pn_ < rhs.pn_);
    }
/*    bool operator==(DihedralParmType const& rhs) const {
      return ( FEQ(pk_, rhs.pk_) &&
               FEQ(pn_, rhs.pn_) &&
               FEQ(phase_, rhs.phase_) &&
               FEQ(scee_, rhs.scee_) &&
               FEQ(scnb_, rhs.scnb_) );
    }
    bool operator<(DihedralParmType const& rhs) const {
      if (FEQ(pk_, rhs.pk_)) {
        if (FEQ(pn_, rhs.pn_)) {
          if (FEQ(phase_, rhs.phase_)) {
            if (FEQ(scee_, rhs.scee_)) {
              return (scnb_ < rhs.scnb_);
            } else return (scee_ < rhs.scee_);
          } else return (phase_ < rhs.phase_);
        } else return (pn_ < rhs.pn_);
      } else return (pk_ < rhs.pk_);
    }*/
  private:
    double pk_;    ///< Dihedral force constant
    double pn_;    ///< Dihedral periodicity
    double phase_; ///< Dihedral phase shift
    double scee_;  ///< 1-4 electrostatics scale factor
    double scnb_;  ///< 1-4 vdW scale factor
};
/// Hold Array of dihedral parameters
class DihedralParmArray {
    typedef std::vector<DihedralParmType> DPArray;
  public:
    /// CONSTRUCTOR
    DihedralParmArray() {}
    /// CONSTRUCTOR - Number of parameters, parameter
    DihedralParmArray(unsigned int n, DihedralParmType const& dp) :
      dihedralparm_(n, dp) {}
    /// Resize the dihedral parm array
    void resize(unsigned int n) { dihedralparm_.resize( n ); }
    /// Clear the dihedral parm array
    void clear() { dihedralparm_.clear(); }
    /// \return reference to specified dihedral parameter
    DihedralParmType& operator[](unsigned int idx) { return dihedralparm_[idx]; }
    /// Add dihedral parameter
    void push_back( DihedralParmType const& bp ) { dihedralparm_.push_back( bp ); }
    /// Iterator
    typedef DPArray::iterator iterator;
    /// \return iterator to beginning
    iterator begin() { return dihedralparm_.begin(); }
    /// \return iterator to end
    iterator end() { return dihedralparm_.end(); }

    /// \return const reference to specified dihedral parameter
    DihedralParmType const& operator[](unsigned int idx) const { return dihedralparm_[idx]; }
    /// \return true if no dihedral parameters
    bool empty() const { return dihedralparm_.empty(); }
    /// \return number of dihedral parameters
    size_t size() const { return dihedralparm_.size(); }
    /// Const iterator
    typedef DPArray::const_iterator const_iterator;
    /// \return const iterator to beginning
    const_iterator begin() const { return dihedralparm_.begin(); }
    /// \return const iterator to end
    const_iterator end() const { return dihedralparm_.end(); }
    /// \return Underlying array
    std::vector<DihedralParmType> const& Array() const { return dihedralparm_; }
    /// \return Last dihedral parameter added
    DihedralParmType const& back() const { return dihedralparm_.back(); }
    /// \return First dihedral parameter added
    DihedralParmType const& front() const { return dihedralparm_.front(); }
  private:
    DPArray dihedralparm_;
};
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
    int& ChangeA1() { return a1_; }
    int& ChangeA2() { return a2_; }
    int& ChangeA3() { return a3_; }
    int& ChangeA4() { return a4_; }
    /// \return type based on skip 1-4 (end) and improper status
    inline Dtype Type() const {
      if (skip14_ && improper_) return BOTH;
      else if (skip14_)         return END;
      else if (improper_)       return IMPROPER;
      else                      return NORMAL;
    }
    inline bool Skip14()     const { return skip14_; }
    inline bool IsImproper() const { return improper_; }
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
    /// \return true if any atom indices do not match
    bool operator!=(DihedralType const& rhs) const {
      return (a1_ != rhs.a1_ ||
              a2_ != rhs.a2_ ||
              a3_ != rhs.a3_ ||
              a4_ != rhs.a4_);
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
/// Hold array of dihedral parameters
class DihedralArray {
    typedef std::vector<DihedralType> DArray;
  public:
    /// CONSTRUCTOR
    DihedralArray() {}

    /// iterator
    typedef DArray::iterator iterator;
    /// begin
    iterator begin() { return dihedrals_.begin(); }
    /// end
    iterator end()   { return dihedrals_.end();   }
    /// const iterator
    typedef DArray::const_iterator const_iterator;
    /// const begin
    const_iterator begin() const { return dihedrals_.begin(); }
    /// const end
    const_iterator end()   const { return dihedrals_.end();   }

    /// Reserve space for # of dihedrals
    void reserve(size_t n) { dihedrals_.reserve(n); }
    /// Add dihedral
    void push_back(DihedralType const& b) { dihedrals_.push_back(b); }
    /// Clear dihedrals
    void clear() { dihedrals_.clear(); }
    /// Erase given dihedral from array
    void erase( iterator bnd ) { dihedrals_.erase( bnd ); }

    /// \return true if no dihedrals
    bool empty()  const { return dihedrals_.empty(); }
    /// \return number of dihedrals
    size_t size() const { return dihedrals_.size(); }
    /// \return specified dihedral
    DihedralType const& operator[](size_t idx) const { return dihedrals_[idx]; }
    /// \return last dihedral added
    DihedralType const& back() const { return dihedrals_.back(); }
    /// \return underlying array
    std::vector<DihedralType> const& Array() const { return dihedrals_; }
  private:
    DArray dihedrals_;
};
// ----- NON-BONDED PARAMETERS -------------------------------------------------
/// Hold LJ 10-12 hbond params
class HB_ParmType {
    /** Tolerance for comparison. A little larger than SMALL because A
      * and B tend to be large.
      */
    // NOTE: Probably should check __cpluscplus here instead of using a
    //       define, but this is guaranteed to be portable.
#     define tol_ 0.00000001
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
    /// \return True if A, B, and HBcut match
    bool operator==(HB_ParmType const& rhs) const {
      return ( (fabs(asol_ - rhs.asol_) < tol_) &&
               (fabs(bsol_ - rhs.bsol_) < tol_) &&
               (fabs(hbcut_ - rhs.hbcut_) < tol_) );
    }
    /// \return True if A, B, or HBcut do not match
    bool operator!=(HB_ParmType const& rhs) const {
      return ( (fabs(asol_ - rhs.asol_) > tol_) ||
               (fabs(bsol_ - rhs.bsol_) > tol_) ||
               (fabs(hbcut_ - rhs.hbcut_) > tol_) );
    }
    /// \return True if A less than zero, or B if A is equal. TODO add hbcut?
    bool operator<(HB_ParmType const& rhs) const {
      if (*this != rhs) {
        if ( (fabs(asol_ - rhs.asol_) < tol_) )
          return (bsol_ < rhs.bsol_);
        else
          return (asol_ < rhs.asol_);
      } else
        return false;
    }
  private:
    double asol_;
    double bsol_;
    double hbcut_;
#   undef tol_
};
typedef std::vector<HB_ParmType> HB_ParmArray;
/// Hold Lennard-Jones 6-12 interaction A and B parameters
class NonbondType {
    /** Tolerance for comparison. Using relative error here since A and
      * B tend to be large (particularly A) and are different magnitudes
      * from each other. Not using Constants::SMALL for the same reason.
      */
    // NOTE: Probably should check __cpluscplus here instead of using a
    //       define, but this is guaranteed to be portable.
#     define tol_ 0.00000001
      //static const double tol_ = 0.00000001;
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
      return ( RELERR(A_, rhs.A_) < tol_ &&
               RELERR(B_, rhs.B_) < tol_ );
    }
    /// \return True if A and B do not match
    bool operator!=(NonbondType const& rhs) const {
      return ( RELERR(A_, rhs.A_) > tol_ ||
               RELERR(B_, rhs.B_) > tol_ );
    }
    /// \return True if A less than zero, or B if A is equal.
    bool operator<(NonbondType const& rhs) const {
      if (*this != rhs) {
        if ( RELERR(A_, rhs.A_) < tol_ )
          return (B_ < rhs.B_);
        else
          return (A_ < rhs.A_);
      } else
        return false;
    }
  private:
    double A_; ///< The coefficient for the r^12 term
    double B_; ///< The coefficient for the r^6 term
#   undef tol_
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
      return ( FEQ(radius_, rhs.radius_) &&
               FEQ(depth_,  rhs.depth_) );
    }
    bool operator!=(LJparmType const& rhs) const {
      return ( FNE(radius_, rhs.radius_) ||
               FNE(depth_,  rhs.depth_) );
    }
    /// \return true if radius and well depth are less in that order
    bool operator<(LJparmType const& rhs) const {
      if (*this != rhs) {
        if (FEQ(radius_, rhs.radius_))
          return (depth_ < rhs.depth_);
        else
          return (radius_ < rhs.radius_);
      } else
        return false;
    }
    /// \return LJ A/B params using Lorentz-Berthelot rules.
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
    bool Has_C_Coeff()                   const { return !ccoef_.empty(); }
    std::vector<int> const& NBindex()    const { return nbindex_;    }
    /// \return Array of LJ 6-12 A and B parameters
    NonbondArray     const& NBarray()    const { return nbarray_;    }
    /// \return Array of LJ 6-12 1-4 A and B parameters
    NonbondArray     const& LJ14()       const { return lj14_;       }
    /// \return Array of LJ 10-12 (hbond) parameters
    HB_ParmArray     const& HBarray()    const { return hbarray_;    }
    /// \return Array of LJ 12-6-4 C parameters
    std::vector<double> const& LJC_Array() const { return ccoef_; }
    /// \return LJ 6-12 A and B parameter at specified index
    NonbondType const& NBarray(int i)    const { return nbarray_[i]; }
    /// \return LJ 6-12 1-4 A and B parameter at specified index
    NonbondType const& LJ14(int i)       const { return lj14_[i]; }
    /// \return LJC parameter at specified index
    double LJC_Array(int idx)            const { return ccoef_[idx]; }
    /// \return LJ 10-12 (hbond) parameter at specified index
    HB_ParmType const& HBarray(int i)    const { return hbarray_[i]; }
    /// In Amber, index < 0 means HB, otherwise LJ 6-12
    int GetLJindex(int type1, int type2) const {
      return nbindex_[ ntypes_ * type1 + type2 ];
    }
    /// Set number of types and init nonbond index array.
    void SetNtypes(unsigned int n) {
      ntypes_ = n;
      nbindex_.assign((size_t)ntypes_ * (size_t)ntypes_, -1); 
    }
    /// Set number of types, init NB index array, init LJ array.
    void SetupLJforNtypes(int n) { SetNtypes(n); nbarray_.assign((n*(n+1))/2, NonbondType()); }
    /// Set number of LJ 1-4 terms
    void SetNLJ14terms(int n)    { lj14_.assign( n, NonbondType() ); }
    /// Set number of LJC terms
    void SetNLJCterms(int n)     { ccoef_.assign( n, 0 ); }
    /// Create an array of LJC terms the same size as LJ A/B all set to zero
    void AllocateLJC() { ccoef_.assign( nbarray_.size(), 0.0 ); }
    /// Set specified LJ term
    NonbondType& SetLJ(int i)    { return nbarray_[i];                  }
    /// Set specified LJ 1-4 term.
    /** Reserve space if not yet allocated, allows it to be used in conjunction
      * with AddLJterm.
      */
    NonbondType& SetLJ14(int ndx)  {
      if (ndx >= (int)lj14_.size())
        lj14_.resize(ndx+1);
      return lj14_[ndx];
    }
    /// Set specified LJC term
    void SetLJC(int i, double ljc) { ccoef_[i] = ljc; }
    /// Set number of HB terms and init HB array TODO combine with SetNtypes?
    void SetNHBterms(int n)   { hbarray_.assign( n, HB_ParmType() ); }
    /// Set specified HB term
    HB_ParmType& SetHB(int i) { return hbarray_[i];                  }
    /// Add a LJ C parameter
    void AddLJC(double c) { ccoef_.push_back( c ); }
    /// Set specified nbindex location to given value.
    void SetNbIdx(int idx, int nbidx) { nbindex_[idx] = nbidx; }
    /// Add given LJ term to nonbond array and update nonbond index array.
    /** Certain routines in sander (like the 1-4 calcs) do NOT use the 
      * nonbond index array; instead they expect the nonbond arrays to be
      * indexed like '(ibig*(ibig-1)/2+isml)', where ibig is the larger atom
      * type index.
      * \return The index into nbarray_
      */
    int AddLJterm(int type1, int type2, NonbondType const& LJ) {
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
      return ndx;
    }
    /// Add given HB term to HB array and update the nonbond index array.
    void AddHBterm(int type1, int type2, HB_ParmType const& HB) {
      int ndx = -((int)hbarray_.size())-1;
      nbindex_[ntypes_ * type1 + type2] = ndx;
      nbindex_[ntypes_ * type2 + type1] = ndx;
      hbarray_.push_back( HB );
    }
    void Clear() { ntypes_ = 0; nbindex_.clear(); nbarray_.clear(); lj14_.clear(); hbarray_.clear(); }
  private:
    int ntypes_;               ///< Number of unique atom types
    std::vector<int> nbindex_; ///< Hold indices into arrays nbarray/hbarray for atom type pairs
    NonbondArray nbarray_;     ///< Hold Lennard-Jones 6-12 A and B parameters for all pairs.
    NonbondArray lj14_;        ///< Lennard-Jones 6-12 1-4 parameters
    HB_ParmArray hbarray_;     ///< Hold 10-12 Amber HBond params for all pairs.
    std::vector<double> ccoef_; ///< Hold Lennard-Jones C parameters for 12-6-4 LJ potential.
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
    void Allocate(unsigned int natomsIn, unsigned int ntypesIn) {
      ntypes_ = ntypesIn;
      ncopies_ = 0;
      array_.clear();
      array_.resize( natomsIn );
      fac_.clear();
      fac_.resize( (size_t)ntypes_ * (size_t)ntypes_ );
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
/// Hold CMAP grid parameters // FIXME use NameType
class CmapGridType {
  public:
    CmapGridType() : nCmapRes_(0), resolution_(0) {}
    CmapGridType(unsigned int r) : nCmapRes_(0), resolution_(r), grid_(r*r, 0.0) {}
    /// \return Grid resolution (in 1 dim, full res is resolutionXresolution)
    unsigned int Resolution()                const { return resolution_;       }
    /// \return CMAP grid
    std::vector<double> const& Grid()        const { return grid_;             }
    /// \return array of residue names this CMAP applies to
    std::vector<std::string> const& ResNames() const { return resNames_; }
    /// \return true if given name matches a residue name
/*    bool MatchesResName(std::string const& nameIn) const {
      for (std::vector<std::string>::const_iterator it = resNames_.begin(); it != resNames_.end(); ++it)
        if (nameIn == *it) return true;
      return false;
    }*/
    /// \return array of atom names this CMAP applies to
    std::vector<std::string> const& AtomNames() const { return atomNames_; }
    /// \return Expected number of CMAP residue names
    int NcmapResNames()                      const { return nCmapRes_; }
    /// \return Grid size as integer, used for topology write
    int Size()                               const { return (int)grid_.size(); }
    /// \return CMAP title
    std::string const& Title()               const { return title_; }
    /// Size of the CMAP grid in bytes
    unsigned int DataSize()                  const { return (sizeof(unsigned int) + (grid_.size()*sizeof(double))); }
    /// Set specified grid point
    void SetGridPt(int idx, double d)              { grid_[idx] = d;           }
    /// Prepare grid for given resolution
    void SetResolution(unsigned int r) {
      grid_.clear();
      grid_.reserve( r*r );
      resolution_ = r;
    }
    /// Add value to the grid
    void AddToGrid(double d) { grid_.push_back( d ); }
    /// Set CMAP title
    void SetTitle(std::string const& t) { title_ = t; }
    /// Number of residues this CMAP will be applied to
    void SetNumCmapRes(int n) {
      resNames_.clear();
      resNames_.reserve( n );
      nCmapRes_ = n;
    }
    /// Add residue name this CMAP will apply to
    void AddResName(std::string const& n) { resNames_.push_back( n ); }
    /// Add atom name this CMAP will apply to
    void AddAtomName(std::string const& n) { atomNames_.push_back( n ); }
    /// \return True if the CMAP is valid
    bool CmapIsValid() const {
      if (resolution_ == 0 || resolution_*resolution_ != grid_.size())
        return false;
      if (atomNames_.size() != 5)
        return false;
      if (resNames_.empty())
        return false;
      return true;
    }
    /// \return True if CMAP # res names matches expected # of names.
    bool CmapNresIsValid() const {
      if (nCmapRes_ == 0 || nCmapRes_ != (int)resNames_.size())
        return false;
      return true;
    }
    /// \return True if CMAP is completely empty
    bool empty() const {
      return (nCmapRes_ == 0 && resolution_ == 0 &&
              grid_.empty() && title_.empty() && resNames_.empty());
    }
    /// \return True if CMAP grid matches given grid. Expensive, so use sparingly.
    bool GridMatches(CmapGridType const& rhs) const {
      if (grid_.size() != rhs.grid_.size()) return false;
      for (unsigned int idx = 0; idx != grid_.size(); idx++) {
        if (FNE(grid_[idx], rhs.grid_[idx])) return false;
      }
      return true;
    }
  private:
    int nCmapRes_;             ///< Number of expected residues this CMAP will apply to
    unsigned int resolution_;  ///< Number of steps along each phi/psi CMAP axis
    std::vector<double> grid_; ///< CMAP grid (size is resolution_*resolution_)
    std::string title_;        ///< CMAP title (from parameter file)
    std::vector<std::string> resNames_;  ///< Residue name(s) this CMAP will apply to
    std::vector<std::string> atomNames_; ///< 5x atom names this CMAP will apply to
};
/// Hold an array of CMAP grids
class CmapGridArray {
    typedef std::vector<CmapGridType> CGarray;
  public:
    CmapGridArray() {}
    /// \return number of CMAP grids
    size_t size() const { return cmapgrids_.size(); }
    /// \return CMAP grid at specified index
    CmapGridType const& operator[](unsigned int idx) const { return cmapgrids_[idx]; }
    /// const iterator
    typedef CGarray::const_iterator const_iterator;
    /// const begin iterator
    const_iterator begin() const { return cmapgrids_.begin(); }
    /// const end iterator
    const_iterator end() const { return cmapgrids_.end(); }
    /// \return True if no CMAP grids
    bool empty() const { return cmapgrids_.empty(); }

    /// \return CMAP grid at specified index
    CmapGridType& operator[](unsigned int idx) { return cmapgrids_[idx]; }
    /// Add CMAP grid
    void push_back(CmapGridType const& cg) { cmapgrids_.push_back( cg ); }
    /// Clear all grids
    void clear() { cmapgrids_.clear(); }
    /// Reserve for specified # of grids
    void reserve(size_t n) { cmapgrids_.reserve(n); }
    /// Resize for specified # of grids
    void resize(size_t n) { cmapgrids_.resize(n); }
    /// iterator
    typedef CGarray::iterator iterator;
    /// begin iterator
    iterator begin() { return cmapgrids_.begin(); }
    /// end iterator
    iterator end() { return cmapgrids_.end(); }

  private:
    CGarray cmapgrids_;
};
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

    void SetIdx(int i) { idx_ = i; }
  private:
    int a1_;
    int a2_;
    int a3_;
    int a4_;
    int a5_;
    int idx_;
};
/// Hold an array of CMAP terms
class CmapArray {
    typedef std::vector<CmapType> CArray;
  public:
    CmapArray() {}
    /// \return number of CMAP terms
    size_t size() const { return cmap_.size(); }
    /// const iterator
    typedef CArray::const_iterator const_iterator;
    /// const begin iterator
    const_iterator begin() const { return cmap_.begin(); }
    /// const end iterator
    const_iterator end() const { return cmap_.end(); }
    /// \return True if no CMAP terms
    bool empty() const { return cmap_.empty(); }

    /// Add CMAP term
    void push_back(CmapType const& cm) { cmap_.push_back( cm ); }
    /// Clear all terms 
    void clear() { cmap_.clear(); }
    /// iterator
    typedef CArray::iterator iterator;
    /// begin iterator
    iterator begin() { return cmap_.begin(); }
    /// end iterator
    iterator end() { return cmap_.end(); }
  private:
    CArray cmap_;
};
/// Hold CHAMBER parameters
/*
class ChamberParmType {
    typedef std::vector<std::string> Sarray;
  public:
    ChamberParmType() {}

    Sarray            const& Description()  const { return chmff_desc_;   }
    BondArray         const& UB()           const { return ub_;           }
    BondParmArray     const& UBparm()       const { return ubparm_;       }
    DihedralArray     const& Impropers()    const { return impropers_;    }
    DihedralParmArray const& ImproperParm() const { return improperparm_; }
    NonbondArray      const& LJ14()         const { return lj14_;         }
    NonbondType& SetLJ14(int idx)                 { return lj14_[idx];    }
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
    void Clear() {
      chmff_desc_.clear(); ub_.clear(); ubparm_.clear();
      impropers_.clear(); improperparm_.clear(); lj14_.clear();
    }
    /// \return true if any CHARMM parameters are set (based on indices).
    bool HasChamber() const {
      if (!ub_.empty()) return true;
      if (!impropers_.empty()) return true;
      if (!lj14_.empty()) return true;
      return false;
    }
  private:
    Sarray chmff_desc_;              ///< CHARMM FF descriptions
    BondArray ub_;                   ///< Urey-Bradley terms
    BondParmArray ubparm_;           ///< Urey-Bradley parameters
    DihedralArray impropers_;        ///< Improper terms
    DihedralParmArray improperparm_; ///< Improper parameters
    NonbondArray lj14_;              ///< Lennard-Jones 1-4 parameters
};*/
#endif
