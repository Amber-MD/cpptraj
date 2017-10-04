#ifndef INC_CONTROL_H
#define INC_CONTROL_H
#include "CpptrajState.h"
#include "VariableArray.h"
/// Control structures.
class Control : public DispatchObject {
  public:
    typedef std::vector<ArgList> ArgArray;
    typedef ArgArray::const_iterator const_iterator;
    /// Hold variable/value pairs
    typedef VariableArray Varray;

    Control() : DispatchObject(CONTROL) {}
    /// Set up control structure.
    virtual int SetupControl(CpptrajState&, ArgList&, Varray&) = 0;
    /// \return true if this is a control block, false otherwise.
    virtual bool IsBlock() const = 0;

    // ----- BLOCK -------------------------------
    /// \return Description of control structure.
    std::string const& Description() const { return description_; }
    /// Check for control structure end command.
    virtual bool EndControl(ArgList const&) const = 0;
    /// Add command to control structure.
    virtual void AddCommand(ArgList const&) = 0;
    /// Control states
    enum DoneType { DONE = 0, NOT_DONE, ERROR };
    /// \return iterator to first command in the block.
    virtual const_iterator begin() const = 0;
    /// \return iterator to last command in the block.
    virtual const_iterator end() const = 0;
    /// Start control block. Init internal variables if necessary.
    virtual void Start() = 0;
    /// Add/update variables and increment, check control state.
    virtual DoneType CheckDone(Varray&) = 0;
  protected:
    std::string description_; ///< Describe control TODO private?
};

/// Loop over mask expression or integer 
class Control_For : public Control {
  public:
    Control_For() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_For(); }

    int SetupControl(CpptrajState&, ArgList&, Varray&);
    bool EndControl(ArgList const& a) const { return (a.CommandIs("done")); }
    void AddCommand(ArgList const& c) { commands_.push_back(c); }
    bool IsBlock() const { return true; }

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

/// Create/update script variables
class Control_Set : public Control {
  public:
    Control_Set() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_Set(); }

    int SetupControl(CpptrajState&, ArgList&, Varray&);
    bool IsBlock() const { return false; }

    bool EndControl(ArgList const& a) const { return false; }
    void AddCommand(ArgList const& c) {}
    const_iterator begin() const { return const_iterator(); }
    const_iterator end()   const { return const_iterator(); }
    void Start() {}
    DoneType CheckDone(Varray&) { return ERROR; }
};
#endif
