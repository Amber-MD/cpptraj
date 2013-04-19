#ifndef INC_ACTION_CREATERESERVOIR_H
#define INC_ACTION_CREATERESERVOIR_H
#include "Action.h"
#include "Traj_AmberNetcdf.h"
// Class: Action_CreateReservoir
/// Create a RREMD structure reservoir.
class Action_CreateReservoir : public Action {
  public:
    Action_CreateReservoir();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CreateReservoir(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    Traj_AmberNetcdf reservoir_;
    Topology* original_trajparm_;
    DataSet* ene_;
    DataSet* bin_;
    std::string filename_;
    bool trajIsOpen_;
    size_t nframes_;
};
#endif
