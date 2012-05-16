#ifndef INC_PDBTYPE_H
#define INC_PDBTYPE_H
#include "Topology.h" // Atom, Residue
#include "FileIO.h"
/// Used to access PDB files
class PDBtype {
  public:
    PDBtype();

    enum PDB_RECTYPE {PDBATOM=0, PDBHET, PDBTER};
  
    bool PDB_GetNextRecord(FileIO *);

    bool IsPDBkeyword();
    bool ID(FileIO *);
    bool IsPDBatomKeyword();
    bool IsPDB_TER();
    bool IsPDB_END();
    Atom pdb_Atom();
    Residue pdb_Residue();
    void pdb_XYZ(double*);
    const double* XYZ();
    void pdb_write_ATOM(FileIO*,PDB_RECTYPE,int,NameType,NameType,char,int,
                        double,double,double,float,float,char *,bool);
    void pdb_write_ATOM(FileIO*,PDB_RECTYPE,int,NameType,NameType,char,int,
                        double,double,double);
/*    int pdb_atomNumber(char*);
    NameType pdb_atomName(char*);
    NameType pdb_resName(char*);
    char pdb_chainID(char*);
    int pdb_resNum(char*);*/
  private: 
    static const char PDB_RECNAME[][7];
    static const size_t BUF_SIZE_ = 83;
    char buffer_[BUF_SIZE_];
    double XYZ_[3];
};
#endif
