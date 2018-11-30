#ifndef INC_DISPATCHOBJECT_H
#define INC_DISPATCHOBJECT_H
#include "ArgList.h"
/// Abstract base class that all dispatchable objects will inherit.
/** A DispatchObject is the most basic unit of command in cpptraj. This
  * base class only contains a category used to sort the command - all
  * other logic is implemented by child classes.
  */
class DispatchObject {
  public:
    /// Object categories. HIDDEN and DEPRECATED should always be last.
    enum Otype { NONE=0, GENERAL,  SYSTEM, COORDS, TRAJ, PARM, ACTION, ANALYSIS,
                 CONTROL,
                 HIDDEN, DEPRECATED };
    /// CONSTRUCTOR
    DispatchObject() : type_(NONE) {}
    /// CONSTRUCTOR - take object type
    DispatchObject(Otype o) : type_(o) {} 
    /// DESTRUCTOR - virtual since this will be inherited
    virtual ~DispatchObject() {}
    /// Print help for this object to screen.
    virtual void Help() const {}
    /// Print help for this object, optionally with subtopics.
    // TODO make this the only version and virtual
    virtual void Help(ArgList&) const { return Help(); }
    /// \return Pointer to new instance of this object.
    virtual DispatchObject* Alloc() const = 0;
    /// \return Keyword for given object category.
    static const char* ObjKeyword(Otype);
    /// \return Object category
    Otype Type() const { return type_; }
  private:
    Otype type_; ///< The object type.
};
#endif
