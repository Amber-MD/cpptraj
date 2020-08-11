#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
#include "CpptrajFile.h"
#include "Atom.h"
#include "Residue.h"
/// Used to access PDB files
class PDBfile : public CpptrajFile {
  public:
    class SSBOND;
    class Link;
    // NOTE: PDB_RECNAME must correspond with this.
    enum PDB_RECTYPE {ATOM=0, HETATM, CRYST1, TER, END, ANISOU, END_OF_FILE, 
                      CONECT, LINK, MISSING_RES, UNKNOWN};
    /// CONSTRUCTOR
    PDBfile();
    /// Check if either of the first two lines contain valid PDB records.
    static bool ID_PDB(CpptrajFile&);
    /// \return the type of the next PDB record read.
    PDB_RECTYPE NextRecord();
    /// \return ATOM/HETATM alt. loc ID
    char pdb_AltLoc() const;
    /// \return Atom info with name and element for ATOM/HETATM; set altLoc and #.
    Atom pdb_Atom(char&, int&);
    /// \return Atom info with name and element for ATOM/HETATM.
    Atom pdb_Atom() { char al; int n; return pdb_Atom(al,n); }
    /// \return Residue info with name, number, icode, and chainID for ATOM/HETATM.
    Residue pdb_Residue();
    /// Set given XYZ array with coords from ATOM/HETATM record.
    void pdb_XYZ(double*);
    /// Get occupancy and B-factor from ATOM/HETATM record.
    void pdb_OccupancyAndBfactor(float&, float&);
    /// Get charge and radius from PQR ATOM/HETATM record.
    void pdb_ChargeAndRadius(float&, float&);
    /// Set given XYZ array with A/B/C/alpha/beta/gamma from CRYST1 record, verbose.
    void pdb_Box_verbose(double*);
    /// Set given XYZ array with A/B/C/alpha/beta/gamma from CRYST1 record, terse.
    void pdb_Box_terse(double*);
    /// Set given array with atom and #s of bonded atoms from CONECT record.
    int pdb_Bonds(int*);
    /// \return Link record.
    Link pdb_Link();
    /// \return current record type.
    PDB_RECTYPE RecType()         const { return recType_; }
    // -------------------------------------------
    /// Set whether column 21 can be used for 4-letter residue names.
    void SetUseCol21(bool b) { useCol21_ = b; }
    /// \return true if using column 21 for 4-letter residue names.
    bool UseCol21() const { return useCol21_; }
    /// Write PDB record header.
    void WriteRecordHeader(PDB_RECTYPE, int, NameType const&, char,
                           NameType const&, char, int, char, const char*);
    /// Write HETATM record using internal atom numbering
    void WriteHET(int, double, double, double);
    /// Write no-name ATOM record using internal atom numbering
    void WriteATOM(int, double, double, double, const char*, double);
    /// Write ATOM record with given name using internal atom numbering
    void WriteATOM(const char*, int, double, double, double, const char*, double);
    /// Write PDB ATOM/HETATM record, no B-factor, occ, elt, or charge.
    void WriteCoord(PDB_RECTYPE, int, NameType const&, NameType const&, char, int,
                    double, double, double);
    /// Write PDB ATOM/HETATM record, no alt loc, chain ID, icode.
    void WriteCoord(PDB_RECTYPE, int, NameType const&, NameType const&, int,
                         double, double, double, float, float, const char*, int);
    /// Write complete PDB ATOM/HETATM record
    void WriteCoord(PDB_RECTYPE, int, NameType const&, char, NameType const&, char, int,
                    char, double, double, double, float, float, const char *, int, bool);
    /// \return True if coordinate write has overflowed; reset overflow status.
    bool CoordOverflow() { bool stat = coordOverflow_; coordOverflow_ = false; return stat; }
    /// Write ANISOU record.
    void WriteANISOU(int, NameType const&, NameType const&, char, int,
                     const double*, const char *, int);
    /// Write TITLE
    void WriteTITLE(std::string const&);
    /// Write CRYST1
    void WriteCRYST1(const double*, const char*);
    /// Write MODEL
    void WriteMODEL(int);
    /// Write CONECT for an atom
    void WriteCONECT(int, std::vector<int> const&, Atom const&);
    /// Write single CONECT
    void WriteCONECT(int, int);
    /// Write SSBOND record
    void WriteSSBOND(int, SSBOND const&, float);
    /// Write ENDMDL
    void WriteENDMDL();
    /// Write END
    void WriteEND();
  private:
    /// \return true if the first 6 chars of buffer match a PDB keyword
    static bool IsPDBkeyword(std::string const&);
    /// Read box info from CRYST1 record
    void readCRYST1(double*);

    int anum_;               ///< Atom number for writing.
    PDB_RECTYPE recType_;    ///< Current record type.
    bool lineLengthWarning_; ///< True if any read line is shorter than 80 char
    bool coordOverflow_;     ///< True if coords on write exceed field width
    bool useCol21_;          ///< If true, use column 21 for 4-char res name
    /// Recognized PDB record types; corresponds to PDB_RECTYPE
    static const char* PDB_RECNAME[];
};
/// Hold information for an SSBOND record.
class PDBfile::SSBOND {
  public:
    SSBOND();
    SSBOND(int, int, Residue const&, Residue const&);
    SSBOND(SSBOND const&);
    SSBOND& operator=(SSBOND const&);
    const char* name1() const { return name1_; }
    const char* name2() const { return name2_; }
    int Idx1()  const { return idx1_;  }
    int Idx2()  const { return idx2_;  }
    int Rnum1() const { return rnum1_; }
    int Rnum2() const { return rnum2_; }
    char Chain1() const { return chain1_; }
    char Chain2() const { return chain2_; }
    char Icode1() const { return icode1_; }
    char Icode2() const { return icode2_; }
  private:
    // TODO SymOP
    int idx1_;      ///< Index into Topology for first SG
    int idx2_;      ///< Index into Topology for second SG
    int rnum1_;     ///< Residue sequence number 1
    int rnum2_;     ///< Residue sequence number 2
    char chain1_;   ///< Chain ID 1
    char chain2_;   ///< Chain ID 2
    char icode1_;   ///< Residue insertion code 1
    char icode2_;   ///< Residue insertion code 2
    char name1_[4]; ///< Residue name 1
    char name2_[4]; ///< Residue name 2
};
/// Hold information for a LINK record.
class PDBfile::Link {
  public:
    Link();
    Link(const char*, char, const char*, char, int, char,
         const char*, char, const char*, char, int, char);
    Link(Link const&);
    Link& operator=(Link const&);
    int Rnum1() const { return rnum1_; }
    int Rnum2() const { return rnum2_; }
    const char* aname1() const { return aname1_; }
    const char* aname2() const { return aname2_; }
    char AltLoc1()       const { return altloc1_;}
    char AltLoc2()       const { return altloc2_;}
    const char* rname1() const { return rname1_; }
    const char* rname2() const { return rname2_; }
    char Chain1()        const { return chain1_; }
    char Chain2()        const { return chain2_; }
    char Icode1()        const { return icode1_; }
    char Icode2()        const { return icode2_; }
  private:
    // TODO symmop
    int rnum1_;
    int rnum2_;
    char aname1_[5];
    char aname2_[5];
    char altloc1_;
    char altloc2_;
    char rname1_[4];
    char rname2_[4];
    char chain1_;
    char chain2_;
    char icode1_;
    char icode2_;
};
#endif
