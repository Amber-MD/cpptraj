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
/// Hold information for NA base.
class NA_Base {
    /// Type for phosphate/sugar atoms (index into atomIdx_).
    enum PSType { PHOS, O4p, C1p, C2p, C3p, C4p };
  public:
    enum PmethodType { ALTONA=0, CREMER };
    /// Type for each standard NA base.
    enum NAType { UNKNOWN_BASE, ADE, CYT, GUA, THY, URA };
    /// Type of hydrogen bond atom.
    enum HBType { NONE = 0, DONOR, ACCEPTOR };
    NA_Base();
    NA_Base(const NA_Base&);
    NA_Base& operator=(const NA_Base&);
    static NAType ID_BaseFromName(NameType const&);
    int Setup_Base(Topology const&, int, NAType, DataSetList&, std::string const&);
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
#endif  
