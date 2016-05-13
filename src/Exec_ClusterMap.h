#ifndef INC_EXEC_CLUSTERMAP_H
#define INC_EXEC_CLUSTERMAP_H
#include "Exec.h"
#include "DataSet_MatrixFlt.h"
#ifdef TIMER
#include "Timer.h"
#endif
// EXPERIMENTAL ALPHA CODE
class Exec_ClusterMap : public Exec {
  public:
    Exec_ClusterMap();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ClusterMap(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
