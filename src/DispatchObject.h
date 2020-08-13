#ifndef INC_DISPATCHOBJECT_H
#define INC_DISPATCHOBJECT_H
/// Abstract base class that all dispatchable objects will inherit.
/** A DispatchObject is the most basic unit of command in cpptraj. This
  * base class only contains a category used to sort the command - all
  * other logic is implemented by child classes.
  */
class ArgList;
class DispatchObject {
  public:
    /// Object categories. DEPRECATED should always be last.
    enum Otype { NONE=0, GENERAL,  SYSTEM,  COORDS, TRAJ, PARM,
                 ACTION, ANALYSIS, CONTROL, DEPRECATED };
    /// CONSTRUCTOR
    DispatchObject() : type_(NONE), hidden_(false) {}
    /// CONSTRUCTOR - take object type
    DispatchObject(Otype o) : type_(o), hidden_(false) {} 
    /// DESTRUCTOR - virtual since this will be inherited
    virtual ~DispatchObject() {}
    /// Print help for this object to screen.
    virtual void Help() const = 0;
    /// Print help - takes arguments.
    virtual void Help(ArgList&) const { Help(); }
    /// \return Pointer to new instance of this object.
    virtual DispatchObject* Alloc() const = 0;
    /// \return Object category
    Otype Type() const { return type_; }
    /// \return True if object should be hidden.
    bool Hidden() const { return hidden_; }
    /// Set object hidden status.
    void SetHidden(bool h) { hidden_ = h; }
  private:
    Otype type_;  ///< The object type.
    bool hidden_; ///< True if object should be hidden, i.e. should not show up in 'help'.
};
#endif
