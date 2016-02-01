#ifndef INC_DEPRECATED_H
#define INC_DEPRECATED_H
#include "DispatchObject.h"
/// Deprecated are for commands which are no longer used. Only Help() is needed.
class Deprecated : public DispatchObject {
  public:
    Deprecated() : DispatchObject(DEPRECATED) {}
    DispatchObject* Alloc() const { return 0; }
};

class Deprecated_MinDist      : public Deprecated { public: void Help() const; };
class Deprecated_Hbond        : public Deprecated { public: void Help() const; };
class Deprecated_TopSearch    : public Deprecated { public: void Help() const; };
class Deprecated_ParmBondInfo : public Deprecated { public: void Help() const; };
class Deprecated_ParmResInfo  : public Deprecated { public: void Help() const; };
class Deprecated_ParmMolInfo  : public Deprecated { public: void Help() const; };
class Deprecated_AvgCoord     : public Deprecated { public: void Help() const; };
class Deprecated_DihScan      : public Deprecated { public: void Help() const; };
#endif
