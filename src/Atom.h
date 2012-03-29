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

    Atom();
    Atom(int, NameType, double (&)[3]);
    Atom( NameType, double, int, double, int, double, double );
    Atom(const Atom &);
    void swap(Atom &, Atom &);
    Atom &operator=(Atom);

    typedef std::vector<int>::const_iterator bond_iterator;
    inline bond_iterator bondbegin() const {
       return bonds_.begin();
    }
    inline bond_iterator bondend() const {
      return bonds_.end();
    }

    void PrintXYZ();

    void SetName(NameType);
    void SetResNum(int);
    void SetMol(int);
    bool NoMol();

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
    inline NameType Name() {
      return aname_;
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

    void AddBond(int);

    void SetElementFromName();
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
};
#endif
