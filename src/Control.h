#ifndef INC_CONTROL_H
#define INC_CONTROL_H
#include "CpptrajState.h"
#include "VariableArray.h"
/// Control block structures.
class ControlBlock : public DispatchObject {
  public:
    typedef std::vector<ArgList> ArgArray;
    typedef ArgArray::const_iterator const_iterator;
    /// Hold variable/value pairs
    typedef VariableArray Varray;
    /// Control block states
    enum DoneType { DONE = 0, NOT_DONE, ERROR };

    ControlBlock() : DispatchObject(BLOCK) {}
    virtual ~ControlBlock() {}
    /// \return Description of control block.
    std::string const& Description() const { return description_; }
    /// Set up control block.
    virtual int SetupBlock(CpptrajState&, ArgList&) = 0;
    /// Check for control block end command.
    virtual bool EndBlock(ArgList const&) const = 0;
    /// Add command to control block.
    virtual void AddCommand(ArgList const&) = 0;
    /// \return iterator to first command in the block.
    virtual const_iterator begin() const = 0;
    /// \return iterator to last command in the block.
    virtual const_iterator end() const = 0;
    /// Start control block. Init internal variables if necessary.
    virtual void Start() = 0;
    /// Add/update variables and increment, check block state.
    virtual DoneType CheckDone(Varray&) = 0;
  protected:
    std::string description_; ///< Describe control TODO private?
};

/// Loop over mask expression or integer 
class ControlBlock_For : public ControlBlock {
  public:
    ControlBlock_For() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new ControlBlock_For(); }

    int SetupBlock(CpptrajState&, ArgList&);
    bool EndBlock(ArgList const& a) const { return (a.CommandIs("done")); }
    void AddCommand(ArgList const& c) { commands_.push_back(c); }
    const_iterator begin() const { return commands_.begin(); }
    const_iterator end()   const { return commands_.end();   }
    void Start();
    DoneType CheckDone(Varray&);
  private:
    enum ForType {ATOMS=0, RESIDUES, MOLECULES, MOLFIRSTRES, MOLLASTRES, INTEGER, UNKNOWN};
    enum OpType { INCREMENT, DECREMENT, LESS_THAN, GREATER_THAN, NO_OP };
    typedef std::vector<int> Iarray;
    class LoopVar {
      public:
      LoopVar() : varType_(UNKNOWN) {}
      Iarray Idxs_;                ///< (MASK only) Selected atom/residue/molecule indices
      Iarray::const_iterator idx_; ///< (MASK only) Current atom/residue/molecule index
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
// =============================================================================
/// Work with script variables
class Control : public DispatchObject {
  public:
    /// Hold variable/value pairs
    typedef VariableArray Varray;

    Control() : DispatchObject(CONTROL) {}
    virtual ~Control() {}
    virtual CpptrajState::RetType SetupControl(CpptrajState&, ArgList&, Varray&) = 0;
};

/// Create/update script variables
class Control_Set : public Control {
  public:
    Control_Set() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_Set(); }

    CpptrajState::RetType SetupControl(CpptrajState&, ArgList&, Varray&);
};

/// List all variables and values.
class Control_Show : public Control {
  public:
    Control_Show() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_Show(); }

    CpptrajState::RetType SetupControl(CpptrajState&, ArgList&, Varray&);
};
#endif
