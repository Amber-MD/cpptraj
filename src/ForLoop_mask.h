#ifndef INC_FORLOOP_MASK_H
#define INC_FORLOOP_MASK_H
#include <vector>
#include "ForLoop.h"
class Topology;
/// For loop over mask expression
class ForLoop_mask : public ForLoop {
    enum MaskType {ATOMS=0, RESIDUES, MOLECULES, MOLFIRSTRES, MOLLASTRES, NTYPES};
  public:
    ForLoop_mask();

    static void helpText();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    bool EndFor(DataSetList&);
  private:
    typedef std::vector<int> Iarray;

    Iarray Idxs_;                ///< (MASK only) Selected atom/residue/molecule indices
    Iarray::const_iterator idx_; ///< (MASK only) Current atom/residue/molecule index
    MaskType mtype_;             ///< Hold specific mask type
    std::string maskExpr_;       ///< The atom mask expression
    Topology* currentTop_;       ///< Topology to use when setting up mask.
};
#endif
