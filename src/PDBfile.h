#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
#include "CpptrajFile.h"
#include "Atom.h"
#include "Residue.h"
/// Used to access PDB files
class PDBfile : public CpptrajFile {
  public:
    // NOTE: PDB_RECNAME must correspond with this.
    enum PDB_RECTYPE {ATOM=0, HETATM, CRYST1, TER, END, ANISOU, END_OF_FILE, 
                      CONECT, UNKNOWN};
    PDBfile() : anum_(1), recType_(UNKNOWN), lineLengthWarning_(false) {}
    /// Check if either of the first two lines contain valid PDB records.
    static bool ID_PDB(CpptrajFile&);
    /// \return the type of the next PDB record read.
    PDB_RECTYPE NextRecord();
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
    /// Set given XYZ array with A/B/C/alpha/beta/gamma from CRYST1 record.
    void pdb_Box(double*);
    /// Set given array with atom and #s of bonded atoms from CONECT record.
    int pdb_Bonds(int*);
    /// \return current record type.
    PDB_RECTYPE RecType()         const { return recType_; }
    // -------------------------------------------
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
    /// Write ANISOU record.
    void WriteANISOU(int, NameType const&, NameType const&, char, int,
                     int, int, int, int, int, int, const char *, int);
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
    /// Write ENDMDL
    void WriteENDMDL();
    /// Write END
    void WriteEND();
  private:
    /// \return true if the first 6 chars of buffer match a PDB keyword
    static bool IsPDBkeyword(std::string const&);

    int anum_;            ///< Atom number for writing.
    PDB_RECTYPE recType_; ///< Current record type.
    bool lineLengthWarning_; ///< True if any read line is shorter than 80 char
    static const char* PDB_RECNAME[];
};
#endif
