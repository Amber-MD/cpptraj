#ifndef INC_CONTROL_H
#define INC_CONTROL_H
#include "DispatchObject.h"
#include "CpptrajState.h" // TODO do we need to return CpptrajState::RetType?
class CpptrajState;
class ArgList;
class VariableArray;
/// Work with script variables
class Control : public DispatchObject {
  public:
    Control() : DispatchObject(CONTROL) {}
    virtual ~Control() {}
    virtual CpptrajState::RetType SetupControl(CpptrajState&, ArgList&, VariableArray&) = 0;
};

/// Create/update script variables
class Control_Set : public Control {
  public:
    Control_Set() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_Set(); }

    CpptrajState::RetType SetupControl(CpptrajState&, ArgList&, VariableArray&);
};

/// List all variables and values.
class Control_Show : public Control {
  public:
    Control_Show() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Control_Show(); }

    CpptrajState::RetType SetupControl(CpptrajState&, ArgList&, VariableArray&);
};
#endif
