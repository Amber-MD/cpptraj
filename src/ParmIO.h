#ifndef INC_PARMIO_H
#define INC_PARMIO_H
#include "CpptrajFile.h"
#include "AmberParm.h"
// Class: ParmIO
/// Base class that all ParmIO objects inherit from
class ParmIO : public CpptrajFile {
  public:
    ParmIO();
    virtual ~ParmIO() { }
    ParmIO(const ParmIO&);
    ParmIO &operator=(const ParmIO&);
    
    void SetDebug(int);

    virtual int ReadParm(AmberParm &) { return 1; }
    virtual int WriteParm(AmberParm &) { return 1; }
};
#endif
