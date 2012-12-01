#ifndef INC_AXISTYPE_H
#define INC_AXISTYPE_H
#include "Topology.h"
/*! \file AxisType.h
    \brief Hold classes and functions used for NA structure analysis.
 */
/// Hold information for NA base.
class NA_Base {
  public:
     /// Type for each standard NA base.
    enum NAType { UNKNOWN_BASE, ADE, CYT, GUA, THY, URA };
    NA_Base();
    NA_Base(const NA_Base&);
    NA_Base& operator=(const NA_Base&);
    static NAType ID_BaseFromName(NameType const&);
    NA_Base(Topology const&, int, NAType);
    void SetInputFrame(Frame const&);
    void PrintAtomNames();
    NAType Type()               const { return type_;           }
    const char* ResName()       const { return *rname_;         }
    int ResNum()                      { return rnum_;           }
    Frame const& Ref()                { return Ref_;            }
    Frame const& Input()              { return Inp_;            }
    AtomMask const& InputFitMask()    { return inpFitMask_;     }
    AtomMask const& RefFitMask()      { return refFitMask_;     }
    const char* AtomName(int i) const { return *(anames_[i]);   }
    bool HasPatom()             const { return patomidx_ != -1; }
    bool HasO4atom()            const { return o4atomidx_ != -1;}
#   ifdef NASTRUCTDEBUG
    const char* RefName(int i)     { return *(refnames_[i]); }
    int HBidx(int i)         const { return hbidx_[i];       }
#   endif
    const double* HBxyz(int i) const { return Inp_.XYZ(hbidx_[i]); }
    const double* Pxyz()       const { return Inp_.XYZ(patomidx_); }
    const double* O4xyz()      const { return Inp_.XYZ(o4atomidx_);}
  private:
    NameType rname_;                ///< Residue name
    int rnum_;                      ///< Original residue number
    NAType type_;                   ///< Base type.
    Frame Ref_;                     ///< Reference coords.
    std::vector<NameType> anames_;  ///< Atom names (Input)
#   ifdef NASTRUCTDEBUG
    std::vector<NameType> refnames_; ///< Atom names (Ref)
#   endif  
    Frame Inp_;                     ///< Input coords.
    int hbidx_[3];                  ///< Indices of h-bonding atoms
    int patomidx_;                  ///< Index of phosphorous atom if present
    int o4atomidx_;                 ///< Index of O4' atom if present
    AtomMask parmMask_;             ///< Mask corresponding to atoms in parm.
    AtomMask inpFitMask_;           ///< Mask of input atoms to be used in RMS fit.
    AtomMask refFitMask_;           ///< Mask of ref atoms to be used in RMS fit.

    int FindAtom(NameType const&);
};

/// Hold information for axis corresponding to base/base-pair.
class NA_Axis {
  public:
    NA_Axis();
    NA_Axis(Matrix_3x3 const&, Vec3 const&, int);
    NA_Axis(int,int,bool);
    void StoreRotMatrix(Matrix_3x3 const&, Vec3 const&);
    void PrintAxisInfo(const char*);
    void FlipYZ();
    void FlipXY();
    Matrix_3x3 const& Rot() const { return R_;      }
    Vec3 const& Oxyz()      const { return origin_; }
    Vec3 const& Rx()        const { return RX_;     }
    Vec3 const& Ry()        const { return RY_;     }
    Vec3 const& Rz()        const { return RZ_;     }
    int Res1()              const { return residue_number_; }
    int Res2()              const { return second_resnum_;  }
    bool IsAnti()           const { return isAnti_;         }
  private:
    Matrix_3x3 R_;       ///< Rotation matrix for this axis
    Vec3 origin_;        ///< Origin of this axis
    // Rotation X|Y|Z vecs are stored to avoid constant calls to R.
    Vec3 RX_;            ///< Rotation X vector (col 1)
    Vec3 RY_;            ///< Rotation Y vector (col 2)
    Vec3 RZ_;            ///< Rotation Z vector (col 3)
    int residue_number_; ///< Residue number
    int second_resnum_;  ///< Second residue if this is a base pair
    bool isAnti_;        ///< If basepair, true if pair is anti-parallel
};
#endif  
