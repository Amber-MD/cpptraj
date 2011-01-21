#ifndef INC_PDBFILEROUTINES_H
#define INC_PDBFILEROUTINES_H

char *pdb_title(char *);
int pdb_atom(char *);
int pdb_name(char *, char *);
int pdb_resname(char *, char *);
char pdb_chain(char *);
int pdb_resnum(char *);
int pdb_xyz(char *, double *);
double pdb_occ(char *);
double pdb_Bfactor(char *);
char *pdb_lastChar(char *);
char *pdb_elt(char *);
char *pdb_charge(char*);
void pdb_write_ATOM(char *, const char *, int, char *, char *, char, int,
                    double, double, double, float, float, char *);

#endif
