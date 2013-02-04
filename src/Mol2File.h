#ifndef INC_MOL2FILE_H
#define INC_MOL2FILE_H
#include "CpptrajFile.h"
#include "Atom.h"
/// Used to access mol2 files.
class Mol2File : public CpptrajFile {
  public: 
    Mol2File();

    enum TRIPOSTAG { MOLECULE=0, ATOM, BOND, SUBSTRUCT };

    static bool IsMol2Keyword(const char*);
    static bool ID_Mol2(CpptrajFile&);
    /// Scan to the specified TRIPOS section of file.
    int ScanTo( TRIPOSTAG );
    void WriteHeader( TRIPOSTAG );
    /// Read in MOLECULE section of mol2file.
    bool ReadMolecule();
    bool WriteMolecule(bool,int);
    //// Used to only read # atoms in next MOLECULE record.
    int NextMolecule();
    /// Read in the next Mol2 BOND line. Get the indices of the bonded atoms.
    int Mol2Bond(int &, int &);
    /// Read in the next Mol2 ATOM line. Get the X Y and Z coords.
    int Mol2XYZ(double *);
    /// Convert current line to Atom 
    Atom Mol2Atom();
    /// Convert current line to Residue
    NameType Mol2Residue(int&);

    void SetMol2Natoms(int nIn)               { mol2atoms_ = nIn; }
    void SetMol2Nbonds(int nIn)               { mol2bonds_ = nIn; }
    void SetMol2Title(std::string const& tIn) { mol2title_ = tIn; }
    int Mol2Natoms()               { return mol2atoms_; }
    int Mol2Nbonds()               { return mol2bonds_; }
    std::string const& Mol2Title() { return mol2title_; }
  private:
    static const char* TRIPOSTAGTEXT[];
    int mol2debug_;
    int mol2atoms_;
    int mol2bonds_;
    std::string mol2title_;
};
#endif  
