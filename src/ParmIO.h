#ifndef INC_PARMIO_H
#define INC_PARMIO_H
#include "CpptrajFile.h"
#include "AmberParm.h"
// Class: ParmIO
/// Base class that all ParmIO objects inherit from
class ParmIO {
  public:
    // NOTE: Move to file eventually
    ParmIO() { debug=0; }
    virtual ~ParmIO() { }
    
    void SetDebug(int debugIn) { debug = debugIn; }

    virtual int ReadParm(AmberParm &, CpptrajFile &) { return 1; }
    virtual int WriteParm(AmberParm &, CpptrajFile &) { return 1; }
  protected:
    int debug;
};
#endif
