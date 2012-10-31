#ifndef INC_PARMIO_H
#define INC_PARMIO_H
#include "Topology.h"
#include "CpptrajFile.h"
// Class: ParmIO
/// Abstract base class that all ParmIO objects inherit from
class ParmIO {
  public:
    virtual ~ParmIO() { }
    virtual bool ID_ParmFormat(CpptrajFile&) = 0; 
    virtual int ReadParm(std::string const&, Topology&) = 0;
    virtual int WriteParm(std::string const&, Topology const&) = 0;
    virtual void SetDebug(int) = 0;
    typedef ParmIO* (*AllocatorType)();
};
#endif
