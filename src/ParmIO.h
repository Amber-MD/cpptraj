#ifndef INC_PARMIO_H
#define INC_PARMIO_H
#include "ArgList.h"
#include "Topology.h"
#include "CpptrajFile.h"
#include "BaseIOtype.h"
#include "BondSearch.h"
/// Abstract base class that all ParmIO objects inherit from
class ParmIO : public BaseIOtype {
  public:
    ParmIO() : debug_(0), Offset_(0.20) {}
    virtual ~ParmIO() { }
    virtual bool ID_ParmFormat(CpptrajFile&) = 0;
    virtual int processReadArgs(ArgList&) = 0; 
    virtual int ReadParm(FileName const&, Topology&) = 0;
    virtual int processWriteArgs(ArgList&) = 0;
    virtual int WriteParm(FileName const&, Topology const&) = 0;
    void SetDebug(int d)       { debug_ = d;                   }
    void SetOffset(double oIn) { if (oIn > 0.0) Offset_ = oIn; }
    void SetBondSearchType(BondSearchType t) { searchType_ = t; }
  protected:
    int debug_;
    double Offset_;             ///< Distance offset for use in bond search
    BondSearchType searchType_; ///< If bond search needed, which algorithm to use.
};
#endif
