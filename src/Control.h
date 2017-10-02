#ifndef INC_CONTROL_H
#define INC_CONTROL_H
#include "CpptrajState.h"
/// Control structures.
class Control : public DispatchObject {
  public:
    typedef std::vector<ArgList> ArgArray;
    typedef ArgArray::const_iterator const_iterator;
    Control() : DispatchObject(CONTROL) {}
    /// Set up control structure.
    virtual int SetupControl(CpptrajState&, ArgList&) = 0;
    /// Check for control structure end command.
    virtual bool EndControl(ArgList const&) const = 0;
    /// Add command to control structure.
    virtual void AddCommand(ArgList const&) = 0;
    /// \return Description of control structure.
    std::string const& Description() const { return description_; }

    enum DoneType { DONE = 0, NOT_DONE, ERROR };

    /// Hold variable, value
    typedef std::pair<std::string,std::string> VarPair;
    /// Hold variable/value pairs
    typedef std::vector<VarPair> Varray;
    /// \return Number of commands in the block
    virtual unsigned int Ncommands() const = 0;
    /// \return iterator to first command in the block.
    virtual const_iterator begin() const = 0;
    /// \return iterator to last command in the block.
    virtual const_iterator end() const = 0;
    /// Start control block. Check/init variables.
    virtual void Start(Varray&) = 0;
    /// Update variables and increment, check control state.
    virtual DoneType CheckDone(Varray&) = 0;
  protected:
    std::string description_; ///< Describe control TODO private?
};

/// Loop over mask expression etc
class Control_For : public Control {
  public:
    Control_For() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_For(); }

    int SetupControl(CpptrajState&, ArgList&);
    bool EndControl(ArgList const& a) const { return (a.CommandIs("done")); }
    void AddCommand(ArgList const& c) { commands_.push_back(c); }

    unsigned int Ncommands() const { return commands_.size(); }
    const_iterator begin() const { return commands_.begin(); }
    const_iterator end()   const { return commands_.end();   }
    void Start(Varray&);
    DoneType CheckDone(Varray&);
  private:
    enum ForType {ATOMS=0, RESIDUES, MOLECULES, UNKNOWN};
    typedef std::vector<int> Iarray;
    class MaskHolder {
      public:
      MaskHolder() : varType_(UNKNOWN) {}
      Iarray Idxs_;                ///< Hold selected atom/residue/molecule indices
      Iarray::const_iterator idx_; ///< Hold current atom/residue/molecule index
      std::string varname_;
      ForType varType_;
    };
    typedef std::vector<MaskHolder> Marray;
    Marray Masks_;
    ArgArray commands_;
};
#endif
