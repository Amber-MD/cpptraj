#ifndef INC_MOLECULE_H
#define INC_MOLECULE_H
// Class: Molecule
/// Hold information for a molecule
class Molecule {
  public:
    Molecule();
    Molecule(int, int, int);

    void SetFirst(int,int);
    void SetLast(int);
    void SetSolvent();

    inline int BeginAtom() const  { return beginAtom_;   }
    inline int EndAtom() const    { return endAtom_;     } 
    inline int FirstRes()         { return firstResNum_; }
    inline bool IsSolvent() const { return isSolvent_;   }
    inline int NumAtoms() const {
      return (endAtom_ - beginAtom_);
    }

  private:
    int firstResNum_;
    int beginAtom_;
    int endAtom_;
    bool isSolvent_;
};
#endif
