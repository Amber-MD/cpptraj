#ifndef INC_MOL2FILE_H
#define INC_MOL2FILE_H
#include "CpptrajFile.h"
#include "Atom.h"
#include "Residue.h"
#include <map>
/// Used to access mol2 files.
class Mol2File : private CpptrajFile {
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
    Residue Mol2Residue();
    /// Write mol2 atom line: at#, atom, res#, res, coords
    void WriteMol2Atom(int, Atom const&, int, const char*, const double*);
    /// Write mol2 bond line; bond#, atom1, atom2
    void WriteMol2Bond(int, int, int, NameType const&, NameType const&);
    /// Write mol2 substructure line; res#, resname, firstatom
    void WriteMol2Substructure(int, const char*, int);
    /// Clear any existing atom mapping
    void ClearAmberMapping();
    /// Read in mapping from Amber to SYBYL
    int ReadAmberMapping(FileName const&, FileName const&, int);

    void SetMol2Natoms(int nIn)               { mol2atoms_ = nIn;  }
    void SetMol2Nbonds(int nIn)               { mol2bonds_ = nIn;  }
    void SetMol2Title(std::string const& tIn) { mol2title_ = tIn;  }
    int Mol2Natoms()                    const { return mol2atoms_; }
    int Mol2Nbonds()                    const { return mol2bonds_; }
    std::string const& Mol2Title()      const { return mol2title_; }
    // CpptrajFile functions that should be accessible.
    using CpptrajFile::SetupRead;
    using CpptrajFile::SetupWrite;
    using CpptrajFile::SetupAppend;
    using CpptrajFile::OpenFile;
    using CpptrajFile::OpenRead;
    using CpptrajFile::OpenWriteNumbered;
    using CpptrajFile::CloseFile;
    using CpptrajFile::Filename;
    using CpptrajFile::Rewind;
  private:
    static const char* TRIPOSTAGTEXT[];
    int mol2debug_;
    int mol2atoms_;
    int mol2bonds_;
    std::string mol2title_;
    //typedef std::map<NameType, NameType, NameType::Compare> AtypeMap;
    typedef std::map<NameType, NameType> AtypeMap;
    AtypeMap Atype_to_Sybyl_; ///< Map Amber type name to SYBYL type name. 
    typedef std::pair<NameType,NameType> AtomPair;
    typedef std::map<AtomPair, int> BndMap;
    BndMap Apair_to_Bond_; ///< Map pairs of Amber atom type names to SYBYL bond type
    static const char* SYBYL_BOND_[]; ///< Hold SYBYL bond type text: 1, 2, 3, am, ar
};
#endif
