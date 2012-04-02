#ifndef INC_ATOM_H
#define INC_ATOM_H
#include <vector>
#include "NameType.h"
// Class: AtomType
/// Hold information for an atom
class Atom {
  public:
    enum AtomicElementType { UNKNOWN_ELEMENT = 0,
      HYDROGEN,   BORON,     CARBON,   NITROGEN, OXYGEN,    FLUORINE,
      PHOSPHORUS, SULFUR,    CHLORINE, BROMINE,  IRON,      CALCIUM,
      IODINE,     MAGNESIUM, COPPER,   LITHIUM,  POTASSIUM, RUBIDIUM,
      CESIUM,     ZINC,      SODIUM
    };
    static const char AtomicElementName[][3];
    // Constructors and assignment
    Atom();
    Atom(int, NameType, double (&)[3]);
    Atom( NameType, double (&)[3], NameType, double );
    Atom( int, NameType, double, int, double, int, NameType, double, double,int );
    Atom(const Atom &);
    void swap(Atom &, Atom &);
    Atom &operator=(Atom);
    // Iterator over bonded atom #s
    typedef std::vector<int>::const_iterator bond_iterator;
    inline bond_iterator bondbegin() const {
       return bonds_.begin();
    }
    inline bond_iterator bondend() const {
      return bonds_.end();
    }
    // Info functions
    void PrintXYZ();
    void Info();
    // Functions that set internal vars
    inline void SetNum(int num) {
      anum_ = num;
    }
    void SetName(NameType);
    void SetResNum(int);
    void SetMol(int);
    bool NoMol();
    // Inline functions returning internal vars
    inline const char *c_str() const {
      return *aname_;
    }
    inline int Num() const {
      return anum_;
    }
    inline int ResNum() const {
      return resnum_;
    }
    inline AtomicElementType Element() const {
      return element_;
    }
    inline NameType Name() const { // NOTE: return reference?
      return aname_;
    }
    inline void ReplaceAsterisk() {
      aname_.ReplaceAsterisk();
    }
    inline NameType Type() const {
      return atype_;
    }
    inline int Mol() const {
      return mol_;
    }
    inline int Nbonds() const {
      return (int)bonds_.size();
    }
    inline const double *XYZ() const {
      return coords_;
    }
    inline double Mass() const {
      return mass_;
    }
    inline double Charge() const {
      return charge_;
    }
    inline double Radius() const {
      return gb_radius_;
    }
    /// Add atom # to this atoms list of bonded atoms.
    void AddBond(int);
    void ClearBonds();

  private:
    static const int AtomicElementNum[];

    double coords_[3];
    double charge_;
    double mass_;
    double gb_radius_;
    double gb_screen_;
    NameType aname_;
    NameType atype_;
    NameType itree_;
    int irotat_;
    int atype_index_;
    int join_;
    AtomicElementType element_;
    int anum_;
    int resnum_;
    int mol_;
    std::vector<int> bonds_;
    // Store bond indices?

    void SetElementFromName();
};
#endif
