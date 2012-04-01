#ifndef INC_MOL2FILE_H
#define INC_MOL2FILE_H
#include "Topology.h"
#include "FileIO.h"
/// Used to access mol2 files.
class Mol2File {
  public: 
    Mol2File();

    enum TRIPOSTAG { MOLECULE=0, ATOM, BOND, SUBSTRUCT };

    bool IsMol2Keyword();
    bool GetLine(FileIO*);
    int ScanTo( FileIO *, TRIPOSTAG );
    bool ReadMolecule( FileIO*);
    int NextMolecule( FileIO * );
    void Mol2Bond(int &, int &);
    Atom Mol2Atom();
    Residue Mol2Residue();
    void Mol2XYZ(double *);

    void SetMol2Natoms(int);
    void SetMol2Nbonds(int);
    int Mol2Natoms();
    int Mol2Nbonds();
    std::string &Mol2Title();
  private:
    static const char TRIPOSTAGTEXT[][22];
    static const size_t BUF_SIZE_ = 256;
    char buffer_[BUF_SIZE_];
    int mol2debug_;
    int mol2atoms_;
    int mol2bonds_;
    std::string mol2title_;
};
#endif  
