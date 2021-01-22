#ifndef INC_CHARMASK_H
#define INC_CHARMASK_H
#include "MaskToken.h"
class Unit;
/// Atom mask using character array.
/** This is the original way to use the atom mask, useful e.g. when 
  * you need to know what atoms are not selected as well as what atoms
  * are selected.
  * CharMasks currently must be set up using the SetupCharMask() routines
  * in the Topology class.
  */
class CharMask : public MaskTokenArray {
  public:
    CharMask() : nselected_(0) {}
    /// Construct from mask expression
    CharMask(std::string const& e) : nselected_(0) { SetMaskString(e); }
    /// Construct from given char array and # selected atoms.
    CharMask(std::vector<char> const& c, int n) : CharMask_(c), nselected_(n) {}
    /// Construct array with nothing selected
    CharMask(int n) : CharMask_(std::vector<char>(n,UnselectedChar_)), nselected_(0) {}
    /// Initialize mask with SelectedChar_ if true or UnselectedChar_ if false
    void InitCharMask(int, bool);
    /// \return true if given atom is selected. 
    bool AtomInCharMask(int) const;
    /// \return true if any atoms within given range are selected.
    bool AtomsInCharMask(int,int) const;
    /// \return true if any atoms of the unit are selected.
    bool AtomsInCharMask(Unit const&) const;
    /// \return Integer array of selected atoms.
    std::vector<int> ConvertToIntMask() const;
    /// \return size of char mask
    int Natom() const { return (int)CharMask_.size(); }
    /// Add an atom to the mask as selected (true) or not selected (false)
    void AddAtom(bool); 
    /// Assign selected status to atom
    void SelectAtom(int, bool);
    // -------------------------------------------
    /// Print all mask atoms in to a line
    void PrintMaskAtoms(const char*) const;
    /// Set up character mask based on current mask expression. 
    int SetupMask(AtomArrayT const&, ResArrayT const&, MolArrayT const&, const double*);
    /// Reset atom mask
    void ResetMask();
    /// Clear any selected atoms in mask.
    void ClearSelected();
    /// Invert mask selection
    void InvertMask();
    /// \return number of selected atoms.
    int Nselected() const { return nselected_; }
  private:
    std::vector<char> CharMask_; ///< Array of size Natoms, 'T' selected, 'F' otherwise.
    int nselected_; ///< Total # of selected atoms.
};
#endif
