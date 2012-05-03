#ifndef INC_ATOMMASK_H
#define INC_ATOMMASK_H
#include <string>
#include <vector>
#include "MaskToken.h"
// Class: AtomMask
/// Hold info on selected atoms based on mask expression.
/** AtomMask is used to hold an array of integers that represent atom numbers
  * of atoms selected based on a mask string. 
  * First the mask string is set via SetMaskString, where it is converted into
  * mask tokens. Then the actual mask can be set up using a Topology in one of
  * two ways:
  * - Integer mask
  *    Although an array of ints becomes larger than a simple character mask 
  *    once more than 25% of the system is selected, it tends to be faster 
  *    than the character array up until about 80% of the system is selected, 
  *    at which point the speed is comparable. This is the default way to use
  *    AtomMask and is how most of the routines in the Frame class have been
  *    written to use AtomMask.
  * - Character mask
  *    This is the original way to use the atom mask, useful e.g. when 
  *    you need to know what atoms are not selected as well as what atoms
  *    are selected. Unlike the integer mask, the character mask is not
  *    directly accessible by outside routines, only by the AtomInCharMask
  *    routine.
  */
// *** NOTE ***
// invertMask, AddAtom, AddAtoms, and AddAtomRange currently only apply
// to the Selected array.
class AtomMask {
  public:
    AtomMask();
    AtomMask(char*);
    AtomMask(const AtomMask &);
    AtomMask & operator=(const AtomMask&);

    /// AtomMask default iterator
    typedef std::vector<int>::const_iterator const_iterator;
    /// Iterator to the beginning of Selected
    const_iterator begin() const;
    /// Iterator at end of Selected
    const_iterator end() const;
    /// Return number of selected atoms
    int Nselected();
    /// Return selected atom at idx
    const int &operator[](int);
    /// Return original mask expression as char*
    const char *MaskString();
    /// Return original mask expression as std::string
    std::string MaskExpression();
    /// Reset atom mask
    void ResetMask();
    /// Switch char used to denote selected atoms (T->F, F->T)
    void InvertMask();
    /// Return the number of atoms mask has in common with another mask
    int NumAtomsInCommon(AtomMask &);
    /// Add given atom to Selected array 
    void AddAtom(int);
    /// Add a list of atoms to mask
    void AddAtoms(std::vector<int>&);
    /// Add minAtom <= atom < maxAtom to mask
    void AddAtomRange(int,int);
    /// Add atoms in given mask to this mask at positon, update position
    void AddMaskAtPosition(AtomMask &, int &);
    /// Print all mask atoms in to a line
    void PrintMaskAtoms(const char*);
    /// Return true if no atoms selected. 
    bool None();
    /// Set the mask string. If NULL, set * (all)
    int SetMaskString(const char*);
    /// Set up Selected based on given char mask 
    void SetupMask(char*,int,int);
    /// Set up CharMask based on given char mask 
    void SetupCharMask(char*, int, int);
    /// True if given atom is T in CharMask
    bool AtomInCharMask(int atom);
    /// True if mask has no expression
    bool NoMaskString() { return (maskString_.empty()); }
    /// Convert mask type (char->int, int->char)
    int ConvertMaskType();

    typedef std::vector<MaskToken>::const_iterator token_iterator;
    inline token_iterator begintoken() const {
      return maskTokens_.begin();
    }
    inline token_iterator endtoken() const {
      return maskTokens_.end();
    }
  private:
    int debug_;
    std::vector<char> CharMask_; ///< Char array of atoms, T if selected, F if not.
    char maskChar_;              ///< The character used to denote a selected atom (default 'T')
    std::string maskString_;     ///< String specifying atom selection
    /** Number of atoms mask was set-up with. Needed when converting from
      * integer mask to Character mask. */
    int Natom_;
    int nselected_;              ///< Number of selected atoms in mask
    std::vector<int> Selected_;  ///< Int array of selected atom numbers, 1 for each selected atom
    std::vector<MaskToken> maskTokens_;

    int Tokenize();
};
#endif
