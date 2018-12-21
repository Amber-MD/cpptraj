#ifndef INC_ATOMMASK_H
#define INC_ATOMMASK_H
#include "MaskToken.h"
/// Atom mask using integer array. 
/** AtomMask is used to hold an array of integers that represent atom numbers
  * of atoms selected based on a mask expression set via SetMaskString. 
  * Although an array of ints becomes larger than a simple character mask 
  * once more than 25% of the system is selected, it tends to be faster 
  * than the character array up until about 80% of the system is selected, 
  * at which point the speed is comparable. This is the most common way to use
  * atom masks in cpptraj and is what most of the routines in the Frame class
  * have been written to use.
  * AtomMasks must currently be set up using the SetupIntMask() routines
  * in the Topology class.
  */
class AtomMask : public MaskTokenArray {
  public:
    AtomMask() : Natom_(0), maskChar_(SelectedChar_) {}
    /// Construct using mask expression
    AtomMask(std::string const& e) : Natom_(0), maskChar_(SelectedChar_) { SetMaskString(e); }
    /// Construct from given integer array and total # atoms.
    AtomMask(std::vector<int> const& s, int n) : Selected_(s), Natom_(n),
                                                 maskChar_(SelectedChar_) {}
    /// Create mask selecting atoms from begin up to (not including) end.
    AtomMask(int,int);
    /// Create mask with single atom selected.
    AtomMask(int);
    /// \return Internal selected atom array.
    std::vector<int> const& Selected()  const { return Selected_;             }
    /// \return true if masks select the same atoms
    bool operator==(AtomMask const&) const;
    /// AtomMask default iterator
    typedef std::vector<int>::const_iterator const_iterator;
    /// \return const iterator to the beginning of Selected
    const_iterator begin()              const { return Selected_.begin();     }
    /// \return const iterator at end of Selected
    const_iterator end()                const { return Selected_.end();       }
    /// \return last selected atom
    int back()                          const { return Selected_.back();      }
    /// \return selected atom at idx
    const int& operator[](int idx)      const { return Selected_[idx];        }
    /// \return Number of atoms from corresponding topology.
    int NmaskAtoms()                    const { return Natom_;                }
    /// Flip current mask expression.
    void InvertMaskExpression();
    /// Invert current mask
    void InvertMask();
    /// \return the number of atoms mask has in common with another mask
    int NumAtomsInCommon(AtomMask const&) const;
    /// Add atom to Selected array; assumes atoms will be in order.
    void AddSelectedAtom(int i) { Selected_.push_back( i ); }
    /// Add given atom to Selected array 
    void AddAtom(int);
    /// Add a list of atoms to mask
    void AddAtoms(std::vector<int> const&);
    /// Add minAtom <= atom < maxAtom to mask
    void AddAtomRange(int,int);
    /// Add atoms in given mask to this mask at positon, update position
    void AddMaskAtPosition(AtomMask const&, int);
    /// Convert from integer mask to char mask.
    std::vector<char> ConvertToCharMask() const;
    /// Set total number of atoms; needed for conversion to CharMask
    void SetNatoms(int n) { Natom_ = n; }
    // -------------------------------------------
    /// Print all mask atoms in to a line
    void PrintMaskAtoms(const char*) const;
    /// Set up integer mask based on current mask expression. 
    int SetupMask(AtomArrayT const&, ResArrayT const&, MolArrayT const&, const double*);
    /// Reset atom mask
    void ResetMask();
    /// Clear any selected atoms in mask.
    void ClearSelected()  { Selected_.clear();            }
    /// \return number of selected atoms
    int Nselected() const { return (int)Selected_.size(); }
  private:
    std::vector<int> Selected_; ///< Int array of selected atom numbers
    int Natom_;     ///< Total # of atoms mask was set up with, used to convert to CharMask
    char maskChar_; ///< Current character used for selected atoms.
};
#endif
