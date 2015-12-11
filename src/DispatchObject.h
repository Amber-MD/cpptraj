#ifndef INC_DISPATCHOBJECT_H
#define INC_DISPATCHOBJECT_H
/// Abstract base class that all dispatchable objects will inherit.
class DispatchObject {
  public:
    /// Object categories. HIDDEN and DEPRECATED should always be last.
    enum Otype { NONE=0, GENERAL,  SYSTEM, COORDS, TRAJ, PARM, ACTION, ANALYSIS,
                 HIDDEN, DEPRECATED };
    /// CONSTRUCTOR
    DispatchObject() : type_(NONE) {}
    /// CONSTRUCTOR - take object type
    DispatchObject(Otype o) : type_(o) {} 
    /// DESTRUCTOR - virtual since this will be inherited
    virtual ~DispatchObject() {}
    /// Help function for object.
    virtual void Help() const = 0;
    /// \return Pointer to this object
    virtual DispatchObject* Alloc() const = 0;
    /// \return Keyword for given object category.
    static const char* ObjKeyword(Otype);
    /// \return Object type
    Otype Type() const { return type_; }
  private:
    Otype type_; ///< The object type.
};
#endif
