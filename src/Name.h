#ifndef INC_NAME_H
#define INC_NAME_H
/*! \file Name.h
    \brief Definition of type used to hold atom/residue names

  The NAME parameters control how things like AmberParm and PtrajMask 
  (the mask parser) expect strings to look. In the Amber-format parameter
  file the atom, residue, symbol, etc. names are four characters. When stored 
  as a string, the NULL character is required, requiring a size of 5. The size
  has been increased to 6 in cpptraj to accomodate the slightly larger atom
  type names that can be found in MOL2 files.
 */
/// Default size for atom and residue names, 5 + NULL.
#define NAMESIZE 6
#define NAME_DEFAULT "     "
/// Default type for atom names, res names, atom types etc
typedef char NAME[NAMESIZE];

void PadWithSpaces(NAME);
void TrimName(NAME);
void WrapName(NAME);
void ReplaceAsterisk(NAME);

#  ifdef __cplusplus
// Test class
class NAME_ {
  public:
    NAME_();
    NAME_(const NAME_&);
    NAME_& operator=(const NAME_&);

    char *Assign(char*,int);

    const char * operator*() const;

  private:
    static const int NAME_SIZE;
    char c_array[6];
};
#  endif
#endif
