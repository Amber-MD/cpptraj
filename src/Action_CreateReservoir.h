#ifndef INC_ACTION_CREATERESERVOIR_H
#define INC_ACTION_CREATERESERVOIR_H
#include "Action.h"
#include "RemdReservoirNC.h"
#include "DataSet_1D.h"
/// Create a NetCDF structure reservoir for REMD.
class Action_CreateReservoir : public Action {
  public:
    Action_CreateReservoir();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_CreateReservoir(); }
    void Help() const;
  private:
    // ----- Inherited functions -----------------
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
    // -------------------------------------------
    RemdReservoirNC reservoir_;   ///< Output structure reservoir file
    Topology* original_trajparm_; ///< Topology corresponding to frames to be used. 
    DataSet_1D* ene_;             ///< Data set containing potential energies
    DataSet_1D* bin_;             ///< Data set containing cluster bins
    std::string title_;           ///< Hold user-specified title
    double reservoirT_;           ///< The reservoir temperature
    int iseed_;                   ///< Reservoir seed
    FileName filename_;           ///< Reservoir file name
    bool trajIsOpen_;             ///< True if reservoir has been initialized and opened
    bool useVelocity_;            ///< If true, include velocities if present.
    bool useForce_;               ///< If true, include forces if present.
    size_t nframes_;              ///< Use to keep track of number of frames written.
};
#endif
