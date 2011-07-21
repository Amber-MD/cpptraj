#ifndef INC_NAME_H
#define INC_NAME_H
// NAMESIZE: Default size for atom and residue names, 5 + NULL.
// Amber atom/residue names are 4, but some mol2 atom types are larger.
#define NAMESIZE 6
#define NAME_DEFAULT "     "
// Default type for atom names, res names, atom types etc
typedef char NAME[NAMESIZE];
#endif
