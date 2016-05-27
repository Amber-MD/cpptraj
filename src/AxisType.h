#ifndef INC_AXISTYPE_H
#define INC_AXISTYPE_H
#include "Topology.h"
#include "DataSetList.h"
#include "DataSet_1D.h"
/*! \file AxisType.h
    \brief Hold classes and functions used for NA structure analysis.
 */
/// Hold information for axis corresponding to base/base-pair.
class NA_Axis {
  public:
    NA_Axis();
    /// Used to set rotation matrix/origin for base pair axis
    void StoreRotMatrix(Matrix_3x3 const&, Vec3 const&);
    void PrintAxisInfo(const char*) const;
    void FlipYZ();
    void FlipXY();
    Matrix_3x3 const& Rot() const { return R_;      }
    Vec3 const& Oxyz()      const { return origin_; }
    Vec3 const& Rx()        const { return RX_;     }
    Vec3 const& Ry()        const { return RY_;     }
    Vec3 const& Rz()        const { return RZ_;     }
  private:
    Matrix_3x3 R_;       ///< Rotation matrix for this axis
    Vec3 origin_;        ///< Origin of this axis
    // Rotation X|Y|Z vecs are stored to avoid constant calls to R.
    Vec3 RX_;            ///< Rotation X vector (col 1)
    Vec3 RY_;            ///< Rotation Y vector (col 2)
    Vec3 RZ_;            ///< Rotation Z vector (col 3)
};
// -----------------------------------------------------------------------------
// Forward declaration for RefBase so NA_Base::Setup_Base() can work
class RefBase;
/// Hold information for NA base.
class NA_Base {
    /// Type for phosphate/sugar atoms (index into atomIdx_).
    enum PSType { PHOS, O4p, C1p, C2p, C3p, C4p };
  public:
    enum PmethodType { ALTONA=0, CREMER };
    /// Type for each standard NA base.
    enum NAType { UNKNOWN_BASE = 0, ADE, CYT, GUA, THY, URA };
    /// Type of hydrogen bond atom.
    enum HBType { NONE = 0, DONOR, ACCEPTOR };
    NA_Base();
    NA_Base(const NA_Base&);
    NA_Base& operator=(const NA_Base&);
    int Setup_Base(RefBase const&, Residue const&, int,
                   std::vector<Atom> const&, DataSetList&, std::string const&);
    void CalcPucker(int, PmethodType);
    void SetInputFrame(Frame const&);
    void SetC3Idx(int i)                 { c3idx_ = i;             }
    void SetC5Idx(int i)                 { c5idx_ = i;             }
    void PrintAtomNames() const;
    NA_Axis&       Axis()                { return axis_;           } //TODO: Remove?
    NA_Axis const& Axis()          const { return axis_;           }
    NAType Type()                  const { return type_;           }
    int ResNum()                   const { return rnum_;           }
    int C3resIdx()                 const { return c3idx_;          }
    int C5resIdx()                 const { return c5idx_;          }
    char BaseChar()                const { return bchar_;          }
    Frame const& Ref()             const { return Ref_;            }
    Frame const& Input()           const { return Inp_;            }
    AtomMask const& InputFitMask() const { return inpFitMask_;     }
    AtomMask const& RefFitMask()   const { return refFitMask_;     }
    const char* atomName(int i)    const { return *(anames_[i]);   }
    NameType const& AtomName(int i)const { return anames_[i];      }
    std::string const& BaseName()  const { return basename_;       }
    bool HasPatom()                const { return atomIdx_[PHOS] != -1; }
    bool HasO4atom()               const { return atomIdx_[O4p] != -1;  }
#   ifdef NASTRUCTDEBUG
    const char* ResName()       const { return *rname_;         }
    const char* RefName(int i)  const { return *(refnames_[i]); }
#   endif
    int Natom()                const { return Inp_.Natom();             }
    HBType HbondType(int i)    const { return hb_[i];                   }
    const double* HBxyz(int i) const { return Inp_.XYZ(i);              }
    const double* Pxyz()       const { return Inp_.XYZ(atomIdx_[PHOS]); }
    const double* O4xyz()      const { return Inp_.XYZ(atomIdx_[O4p]);  }
    DataSet_1D* Pucker()       const { return pucker_;                  }
  private:
    const double* C1xyz()      const { return Inp_.XYZ(atomIdx_[C1p]);  }
    const double* C2xyz()      const { return Inp_.XYZ(atomIdx_[C2p]);  }
    const double* C3xyz()      const { return Inp_.XYZ(atomIdx_[C3p]);  }
    const double* C4xyz()      const { return Inp_.XYZ(atomIdx_[C4p]);  }
    /// Find index in Input corresponding to atom name.
    int FindAtom(NameType const&) const;

    NA_Axis axis_;                  ///< Reference axis.
    DataSet_1D* pucker_;            ///< Hold sugar pucker data.
    int rnum_;                      ///< Original residue number
    int c3idx_;                     ///< Index of c3' neighbor res.
    int c5idx_;                     ///< Index of c5' neighbor res.
    char bchar_;                    ///< 1 char base name.
    NAType type_;                   ///< Base type.
    Frame Ref_;                     ///< Reference coords.
    std::vector<NameType> anames_;  ///< Atom names (Input)
    std::string basename_;          ///< Base short name
#   ifdef NASTRUCTDEBUG
    NameType rname_;                 ///< Residue name
    std::vector<NameType> refnames_; ///< Atom names (Ref)
#   endif  
    Frame Inp_;                     ///< Input coords.
    std::vector<HBType> hb_;        ///< Hydrogen bond type of each Input atom.
    int atomIdx_[6];                ///< Indices of Input phosphate/sugar atoms.
    AtomMask parmMask_;             ///< Mask corresponding to atoms in parm.
    AtomMask inpFitMask_;           ///< Mask of input atoms to be used in RMS fit.
    AtomMask refFitMask_;           ///< Mask of ref atoms to be used in RMS fit.
};
// -----------------------------------------------------------------------------
/// Hold info for a NA base reference atom.
class NA_Atom {
  public:
    NA_Atom() : x_(0.0), y_(0.0), z_(0.0), hb_type_(NA_Base::NONE), rms_fit_(0) {}
    NA_Atom(double x, double y, double z, NA_Base::HBType t, int r, const char* n) :
      x_(x), y_(y), z_(z), hb_type_(t), rms_fit_(r), aname_(n) {}
    NameType const& Name()    const { return aname_; }
    const char* name()        const { return *aname_; }
    double X()                const { return x_; }
    double Y()                const { return y_; }
    double Z()                const { return z_; }
    NA_Base::HBType HB_type() const { return hb_type_; }
    int RmsFit()              const { return rms_fit_; }
  private:
    double x_, y_, z_;
    NA_Base::HBType hb_type_;
    int rms_fit_;
    NameType aname_;
};
// -----------------------------------------------------------------------------
/// Hold information for a NA base reference.
class RefBase {
    typedef std::vector<NA_Atom> NA_Array;
  public:
    typedef std::vector<NameType> NameArray;
    RefBase() {}
    RefBase(char b, NameType const& n, NA_Base::NAType t) : names_(1,n), baseChar_(b), type_(t) {}
    RefBase(char b, NameArray const& a,NA_Base::NAType t) : names_(a), baseChar_(b), type_(t) {}
    void AddName( NameType const& n ) { names_.push_back( n ); }
    void AddAtom( NA_Atom const& a )  { atoms_.push_back( a ); }
    void PrintInfo() const;
    bool NameMatches(NameType const&) const;
    typedef NA_Array::const_iterator const_iterator;
    const_iterator begin()             const { return atoms_.begin(); }
    const_iterator end()               const { return atoms_.end();   }
    NA_Atom const& operator[](int idx) const { return atoms_[idx];    }
    char BaseChar()                    const { return baseChar_;      }
    NA_Base::NAType Type()             const { return type_;          }
    bool empty()                       const { return atoms_.empty(); }
  private:
    NA_Array atoms_;
    NameArray names_; ///< List of names corresponding to this reference.
    char baseChar_;   ///< Base 1 char name
    NA_Base::NAType type_; ///< Base type
};
// -----------------------------------------------------------------------------
/// Hold information for all defined NA base references
class NA_Reference {
  public:
    enum RetType { BASE_OK = 0, BASE_ERROR, NOT_FOUND };
    NA_Reference();
    /// \return NA_Base set up with correct reference.
    RetType SetupBaseRef(NA_Base&, Topology const&, int, DataSetList&, std::string const&);
    /// Add given name to first reference base of the specified type.
    void AddNameToBaseType(NameType const&, NA_Base::NAType);
    /// Load a reference from a file
    int LoadFromFile(FileName const&);
  private:
    typedef std::vector<RefBase> BaseArray;
    BaseArray bases_;
};
#endif
