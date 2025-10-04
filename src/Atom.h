#ifndef INC_ATOM_H
#define INC_ATOM_H
#include <vector>
#include "NameType.h"
#include "SymbolExporting.h"
/// Hold information for an atom
class Atom {
  public:
    /** Recognized atomic elements. This must be kept in sync with:
      *  AtomicElementNum_
      *  AtomicElementName_
      *  AtomicElementMass_
      *  AtomicElementRadius_
      *  NUMELEMENTS_
      *  MapAtom::AtomicElementChar
      */
    enum AtomicElementType { UNKNOWN_ELEMENT = 0,
      HYDROGEN,   BORON,      CARBON,   NITROGEN,  OXYGEN,     FLUORINE,   //  6
      PHOSPHORUS, SULFUR,     CHLORINE, BROMINE,   IRON,       CALCIUM,    // 12
      IODINE,     MAGNESIUM,  COPPER,   LITHIUM,   POTASSIUM,  RUBIDIUM,   // 18
      CESIUM,     ZINC,       SODIUM,   ALUMINUM,  ARGON,      ARSENIC,    // 24
      SILVER,     GOLD,       ASTATINE, BERYLLIUM, BARIUM,     BISMUTH,    // 30
      CHROMIUM,   COBALT,     CADMIUM,  FRANCIUM,  GALLIUM,    GERMANIUM,  // 36
      HELIUM,     HAFNIUM,    MERCURY,  INDIUM,    IRIDIUM,    KRYPTON,    // 42
      MANGANESE,  MOLYBDENUM, NEON,     NICKEL,    NIOBIUM,    OSMIUM,     // 48
      PALLADIUM,  PLATINUM,   LEAD,     POLONIUM,  RUTHENIUM,  RHODIUM,    // 54
      RHENIUM,    RADON,      RADIUM,   SILICON,   SCANDIUM,   SELENIUM,   // 60
      STRONTIUM,  TIN,        ANTIMONY, TITANIUM,  TECHNETIUM, TELLURIUM,  // 66
      TANTALUM,   THALLIUM,   VANADIUM, TUNGSTEN,  XENON,      ZIRCONIUM,  // 72
      YTTRIUM,    LUTETIUM,   SAMARIUM,
      EXTRAPT,    DRUDE
    };
    // Constructors and assignment ---------------
    Atom();
    virtual ~Atom() {}
    /// Take atom name, and (optional) 2 character element name.
    Atom(NameType const&, const char*);
    /// Take atom name, type name, and charge.
    Atom(NameType const&, NameType const&, double);
    /// Take atom name, type name, and type index.
    Atom(NameType const&, NameType const&, int);
    /// Take atom name, charge, mass, and type name (TODO deprecate)
    Atom(NameType const&, double, double, NameType const&);
    /// Take atom name, charge, mass, type name, and element
    Atom(NameType const&, double, double, NameType const&, AtomicElementType);
    /// Atom name, charge, polarizability, atomic num, mass, type index, type name, gb radius and screen parameters.
    Atom(NameType const&, double, double, int, double, int, NameType const&, double, double);
    Atom(const Atom &);
    void swap(Atom &, Atom &);
    Atom &operator=(Atom);
    // Functions related to bonded atoms --------- 
    typedef std::vector<int>::const_iterator bond_iterator;
    inline bond_iterator bondbegin() const { return bonds_.begin();     }
    inline bond_iterator bondend()   const { return bonds_.end();       }
    inline int Bond(int idx)         const { return bonds_[idx];        }
    std::vector<int> const& BondIdxArray() const { return bonds_; }
    /// Add atom index # to this atoms list of bonded atoms.
    void AddBondToIdx(int idxIn)           { bonds_.push_back( idxIn ); }
    void ClearBonds()                      { bonds_.clear() ;           }
    void RemoveBondToIdx(int);
    void SortBonds();
    // TODO: Use this routine in AtomMap etc
    /// \return true if this atom is bonded to given atom index 
    bool IsBondedTo(int) const;
    // Functions that set internal vars ----------
    void SetResNum(int resnumIn)             { resnum_ = resnumIn;  }
    void SetMol(int molIn)                   { mol_ = molIn;        }
    void SetCharge(double qin)               { charge_ = qin;       }
    void SetPolar(double pin)                { polar_ = pin;        }
    void SetMass(double min)                 { mass_ = min;         }
    void SetGBradius(double rin)             { gb_radius_ = rin;    }
    void SetGBscreen(double sin)             { gb_screen_ = sin;    }
    void SetTypeIndex(int tin)               { atype_index_ = tin;  }
    void SetName(NameType const& nin)        { aname_ = nin;        }
    void SetTypeName(NameType const& tin)    { atype_ = tin;        }
    // Internal vars -----------------------------
    inline bool NoMol()                const { return ( mol_ < 0 ); }
    inline const char *c_str()         const { return *aname_; }
    inline int ResNum()                const { return resnum_; }
    inline AtomicElementType Element() const { return element_; }
    inline int AtomicNumber()          const { return AtomicElementNum_[element_];  }
    inline const char* ElementName()   const { return AtomicElementName_[element_]; }
    inline double ElementRadius()      const { return AtomicElementRadius_[element_]; }
    double PaulingElectroNeg()         const { return PaulingElectroNeg_[element_]; }
    inline const NameType& Name()      const { return aname_; }
    inline const NameType& Type()      const { return atype_; }
    inline int TypeIndex()             const { return atype_index_; }
    inline int MolNum()                const { return mol_; }
    inline int Nbonds()                const { return (int)bonds_.size(); }
    inline double Mass()               const { return mass_; }
    inline double Charge()             const { return charge_; }
    inline double Polar()              const { return polar_; }
    inline double GBRadius()           const { return gb_radius_; }
    inline double Screen()             const { return gb_screen_; }
    /// \return Optimal bond length in Ang. based on given element types.
    static double GetBondLength(AtomicElementType, AtomicElementType);
    /// \return Element name of element
    static const char* ElementName(AtomicElementType e) { return AtomicElementName_[e]; }
    /// \return PARSE radius in Ang. based on element.
    double ParseRadius() const;
    /// Determine element from given atomic number. Use mass/name if number < 1.
    void DetermineElement(int);
    /// Set atomic mass from current element
    void SetMassFromElement();
  protected:
    static const size_t NUMELEMENTS_ = 78;
  private:
    static CPPTRAJ_EXPORT const int AtomicElementNum_[];
    static CPPTRAJ_EXPORT const char* AtomicElementName_[];
    static CPPTRAJ_EXPORT const double AtomicElementMass_[];
    static CPPTRAJ_EXPORT const double AtomicElementRadius_[];
    static CPPTRAJ_EXPORT const double PaulingElectroNeg_[];
    double charge_;    ///< Charge in e-
    double polar_;     ///< Atomic polarizability in Ang^3
    double mass_;      ///< mass in amu
    double gb_radius_; ///< GB radius in Ang
    double gb_screen_; ///< GB screening parameter
    NameType aname_;   ///< Atom name
    NameType atype_;   ///< Atom type name
    int atype_index_;  ///< Atom type index for nonbond params
    AtomicElementType element_;
    int resnum_;       ///< Index into residues array.
    int mol_;          ///< Index into molecules array.
    std::vector<int> bonds_; ///< Indices of atoms bonded to this one.

    static void WarnBondLengthDefault(AtomicElementType, AtomicElementType, double);
    void SetElementFromName();
    void SetElementFromSymbol(char,char);
    void SetElementFromMass();
};
#endif
