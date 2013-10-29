#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
#include "CpptrajFile.h"
#include "Atom.h"
/// Used to access PDB files
class PDBfile : public CpptrajFile {
  public:
    PDBfile() : anum_(1) {}
    // NOTE: PDB_RECNAME must correspond with this.
    enum PDB_RECTYPE {ATOM=0, HETATM, TER, ANISOU};
    /// \return true if the first 6 chars of buffer match a PDB keyword
    static bool IsPDBkeyword(std::string const&);
    /// Check if either of the first two lines contain valid PDB records.
    static bool ID_PDB(CpptrajFile&);
    /// \return true if current line has an ATOM/HETATM record.
    bool IsPDBatomKeyword();
    /// \return true if current line has TER keyword.
    bool IsPDB_TER();
    /// \return true if current line has END keyword.
    bool IsPDB_END();
    /// \return Atom based on current line.
    Atom pdb_Atom();
    /// \return Residue based on current line.
    NameType pdb_Residue(int&);
    /// Set XYZ based on current line.
    void pdb_XYZ(double*);

    /// Write TER record
    void WriteTER(int, NameType const&, char, int);
    /// Write HETATM record using internal atom numbering
    void WriteHET(int, double, double, double);
    /// Write no-name ATOM record using internal atom numbering
    void WriteATOM(int, double, double, double, const char*, double);
    /// Write ATOM record with given name using internal atom numbering
    void WriteATOM(const char*, int, double, double, double, const char*, double);
    /// Write PDB ATOM/HETATM record, no B-factor, occ, elt, or charge.
    void WriteCoord(PDB_RECTYPE, int, NameType const&, NameType const&, char, int,
                    double, double, double);
    /// Write complete PDB ATOM/HETATM record
    void WriteCoord(PDB_RECTYPE, int, NameType const&, NameType const&, char, int,
                    double, double, double, float, float, const char *, int, bool);
    /// Write ANISOU record.
    void WriteANISOU(int, NameType const&, NameType const&, char, int,
                     int, int, int, int, int, int, const char *, int);
  private:
    /// Write PDB record header.
    void WriteRecordHeader(PDB_RECTYPE, int, NameType const&,
                           NameType const&, char, int);

    int anum_;
    static const char* PDB_RECNAME[];
};
#endif
