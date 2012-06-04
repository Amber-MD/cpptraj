#ifndef INC_MOLECULE_H
#define INC_MOLECULE_H
// Class: Molecule
/// Hold information for a molecule
class Molecule {
  public:
    Molecule();
    Molecule(int, int);

    void SetFirst(int begin) { beginAtom_ = begin; }
    void SetLast(int last)   { endAtom_ = last;    }
    void SetSolvent()        { isSolvent_ = true;  }

    inline int BeginAtom() const  { return beginAtom_;   }
    inline int EndAtom() const    { return endAtom_;     } 
    inline bool IsSolvent() const { return isSolvent_;   }
    inline int NumAtoms() const   { return (endAtom_ - beginAtom_); }

  private:
    int beginAtom_;
    int endAtom_;
    bool isSolvent_;
};
#endif
