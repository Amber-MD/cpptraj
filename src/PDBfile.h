#ifndef INC_PDBFILE_H
#define INC_PDBFILE_H
#include "Topology.h" // Atom, Residue
#include "FileIO.h"
/// Used to access PDB files
class PDBfile {
  public:
    PDBfile();

    enum PDB_RECTYPE {PDBATOM=0, PDBHET, PDBTER};
  
    bool PDB_GetNextRecord(FileIO *);

    bool IsPDBkeyword();
    bool IsPDBatomKeyword();
    bool IsPDB_TER();
    bool IsPDB_END();
    Atom pdb_Atom();
    Residue pdb_Residue();
    void pdb_XYZ(double*);
    void pdb_write_ATOM(FileIO*,PDB_RECTYPE,int,NameType,NameType,char,int,
                        double,double,double,float,float,char *,bool);
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
/*
bool isPDBkeyword(char *);
bool isPDBatomKeyword(char *);
char *pdb_title(char *);
int pdb_atom(char *);
NameType pdb_name(char *);
int pdb_resname(char *, char*);
char pdb_chain(char *);
int pdb_resnum(char *);
int pdb_xyz(char *, double *);
double pdb_occ(char *);
double pdb_Bfactor(char *);
char *pdb_lastChar(char *);
char *pdb_elt(char *);
char *pdb_charge(char*);
int pdb_write_ATOM(char *, PDB_RECTYPE, int, char *, char *, char, int,
                    double, double, double, float, float, char *,bool);*/
#endif
