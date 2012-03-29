#ifndef INC_PARMIO_H
#define INC_PARMIO_H
#include "CpptrajFile.h"
#include "Topology.h"
// Class: ParmIO
/// Base class that all ParmIO objects inherit from
class ParmIO : public CpptrajFile {
  public:
    ParmIO();
    virtual ~ParmIO() { }
    ParmIO(const ParmIO&);
    ParmIO &operator=(const ParmIO&);

    void SetDebug(int);

    virtual int ReadParm(Topology &) { return 1; }
    virtual int WriteParm(Topology &) { return 1; }
    virtual bool ID_ParmFormat() { return false; }
};
#endif
