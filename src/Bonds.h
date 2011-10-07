#ifndef INC_BONDS_H
#define INC_BONDS_H
/// Class: BondInfo
/// Used to determine how many unique non-covalently bonded molecules
/// exist within a set of atoms.
class BondInfo {
    // Mol maxbonds #bonds bond1 bond2 bond3 bond4 bond5 bond6 bond7 
    // 0   1        2      3     4     5     6     7     8     9
    typedef int bondinfo[10];
    bondinfo *Molecule;
    int natom;

    int BondAtoms(int,int);
    void VisitAtom(int, int);
  public:
    BondInfo();
    ~BondInfo();
    int Setup(int);
    void SetValence(int,char*);
    int CreateBond(int,int);
    void PrintBonds();
    int *DetermineMolecules(int*);
};
// Other functions
int AtomicNumberFromName(char *);
char ElementFromName(char *);
int MaxValence(char*);
double GetBondedCut(char *, char *);
double GetBondedCut(char, char);
#endif
