#ifdef __cplusplus
extern "C" {
#endif

/*  ________________________________________________________________________
 */

/* From ptraj_local.h */
/* The parameter file assumes the atom, residue, symbol, etc. names to be *
 * four characters (we will store them as strings, requiring a newline,   *
 * hence the size is 5).                                                  *
 * NOTE: This is also defined in AmberParm.h.                             */
#define NAME_SIZE 6
#define NAME_DEFAULT "    "
typedef char Name[NAME_SIZE];


/*
 *  Enhanced atom selection parser: Viktor Hornak, Stony Brook University
 */
#define  MAXSELE 1000
#define  ALL      0
#define  NUMLIST  1
#define  NAMELIST 2
#define  TYPELIST 3
#define  ELEMLIST 4

extern char * parseMaskString(char *, int, int, Name *, Name *, int *, void *, char, Name *, int);


/* DAN ROE: NOTE: Function prototypes are here because they are out of
 * order in mask.c. Sloppy.
 */
int isOperator(char);
int isOperand(char);
int priority(char);
int tokenize(char *, char *);
int torpn(char *, char *);

char * eval(char *, int, int, Name *, Name *, int *, void *, char, Name *);
char * selectDist(char *, char *, int, int, int *, void *, char);
char * binop(char, char *, char *, int);
char * neg(char *, int);
int isElemMatch(char *, char *);
int isNameMatch(char *, char *);
void resnum_select(int, int, char *, int, int *);
void resname_select(char *, char *, int, Name *, int *);
void all_select(char *, int);
void atnum_select(int, int, char *, int);
void atname_select(char *, char *, int, Name *);
void attype_select(char *, char *, int, Name *);
void atelem_select(char *, char *, int, Name *);
void residue_numlist(char *, char *, int, int *);
void residue_namelist(char *, char *, int, Name *, int *);
void atom_numlist(char *, char *, int);
void atom_namelist(char *, char *, int, Name *);
void atom_typelist(char *, char *, int, Name *);
void atom_elemlist(char *, char *, int, Name *);
char * selectElemMask(char *, int, int, Name *, Name *, int *, Name *);
char * parseMaskString(char *, int, int, Name *, Name *, int *, void *, char, Name *, int);

#ifdef __cplusplus
}
#endif
