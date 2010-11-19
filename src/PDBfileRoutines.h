#ifndef INC_PDBFILEROUTINES_H
#define INC_PDBFILEROUTINES_H

char *pdb_title(char *);
int pdb_atom(char *);
char *pdb_name(char *);
char *pdb_resname(char *);
char pdb_chain(char *);
int pdb_resnum(char *);
int pdb_xyz(char *, double *);
double pdb_occ(char *);
double pdb_Bfactor(char *);
char *pdb_lastChar(char *);
char *pdb_elt(char *);
char *pdb_charge(char*);

#endif
