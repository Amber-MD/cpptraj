#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
#include "Topology.h" // Atom, Residue
#include "CpptrajFile.h"
/// Used to access PDB files
class PDBfile : public CpptrajFile {
  public:
    PDBfile() : anum_(1) {}
    // NOTE: PDB_RECNAME must correspond with this.
    enum PDB_RECTYPE {ATOM=0, HETATM, TER/*, END, HEADER, TITLE, COMPND, AUTHOR,
                      CRYST1, REMARK, MODEL, JRNL, SEQRES, BLANK*/};
    /// \return true if the first 6 chars of buffer match a PDB keyword
    static bool IsPDBkeyword(const char*);
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
    Residue pdb_Residue();
    /// Set XYZ based on current line.
    void pdb_XYZ(double*);
    /// Write HETATM record using internal atom numbering
    void WriteHET(int, double, double, double);
    /// Write ATOM record using internal atom numbering
    void WriteATOM(int, double, double, double, const char*, double);
    /// Write TER record
    void WriteTER(int, NameType const&, char, int);
    /// Write PDB ATOM/HETATM record
    void WriteRec(PDB_RECTYPE, int, NameType const&, NameType const&, char, int,
                  double, double, double, float, float, const char *, bool);
    /// Write PDB ATOM/HETATM record
    void WriteRec(PDB_RECTYPE, int, NameType const&, NameType const&, char, int,
                  double, double, double);
/*    int pdb_atomNumber(char*);
    NameType pdb_atomName(char*);
    NameType pdb_resName(char*);
    char pdb_chainID(char*);
    int pdb_resNum(char*);*/
  private:
    int anum_;
    static const char PDB_RECNAME[][7];
};
#endif
