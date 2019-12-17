#ifndef INC_CONTROLBLOCK_FOR_H
#define INC_CONTROLBLOCK_FOR_H
#include "ControlBlock.h"
/// Loop over mask expression or integer 
class ControlBlock_For : public ControlBlock {
  public:
    ControlBlock_For() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new ControlBlock_For(); }

    int SetupBlock(CpptrajState&, ArgList&);
    bool EndBlock(ArgList const&) const;
    void AddCommand(ArgList const& c) { commands_.push_back(c); }
    const_iterator begin() const { return commands_.begin(); }
    const_iterator end()   const { return commands_.end();   }
    void Start();
    DoneType CheckDone(VariableArray&);
  private:
    enum ForType {ATOMS=0, RESIDUES, MOLECULES, MOLFIRSTRES, MOLLASTRES, INTEGER, LIST, UNKNOWN};
    enum OpType { INCREMENT=0, DECREMENT, LESS_THAN, GREATER_THAN, NO_OP };
    typedef std::vector<int> Iarray;
    typedef std::vector<std::string> Sarray;
    class LoopVar {
      public:
      LoopVar() : varType_(UNKNOWN) {}
      Iarray Idxs_;                ///< (MASK only) Selected atom/residue/molecule indices
      Iarray::const_iterator idx_; ///< (MASK only) Current atom/residue/molecule index
      Sarray List_;                ///< (LIST only) List of strings to iterate over.
      Sarray::const_iterator sdx_; ///< (LIST only) Iterator to current list item.
      std::string varname_;        ///< Loop variable name
      ForType varType_;            ///< Loop variable type
      OpType endOp_;               ///< (INTEGER only) end operator
      OpType incOp_;               ///< (INTEGER only) increment operator
      int start_;                  ///< (INTEGER only) initial value
      int end_;                    ///< (INTEGER only) end value
      int inc_;                    ///< (INTEGER only) increment value
      int currentVal_;             ///< (INTEGER only) current value
    };
    typedef std::vector<LoopVar> Marray;
    Marray Vars_;
    ArgArray commands_;
};
#endif
