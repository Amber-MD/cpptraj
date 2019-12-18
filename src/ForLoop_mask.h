#ifndef INC_FORLOOP_MASK_H
#define INC_FORLOOP_MASK_H
#include <vector>
#include "ForLoop.h"
/// For loop over mask expression
class ForLoop_mask : public ForLoop {
  public:
    ForLoop_mask() {}

    int SetupFor(CpptrajState&, std::string const&, ArgList&);
  private:
    enum MaskType {ATOMS=0, RESIDUES, MOLECULES, MOLFIRSTRES, MOLLASTRES, NTYPES};
    typedef std::vector<int> Iarray;

    Iarray Idxs_;                ///< (MASK only) Selected atom/residue/molecule indices
    Iarray::const_iterator idx_; ///< (MASK only) Current atom/residue/molecule index

};
#endif
