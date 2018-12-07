#include "Action_CreateReservoir.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_CreateReservoir::Action_CreateReservoir() :
  original_trajparm_(0),
  ene_(0),
  bin_(0),
  reservoirT_(0.0),
  iseed_(0),
  trajIsOpen_(false),
  useVelocity_(false),
  useForce_(false),
  nframes_(0)
{}

void Action_CreateReservoir::Help() const {
  mprintf("\t<filename> ene <energy data set> [bin <cluster bin data set>]\n"
          "\ttemp0 <temp0> iseed <iseed> [novelocity] [noforce]\n"
          "\t[%s] [title <title>]\n", DataSetList::TopArgs);
  mprintf("  Create structure reservoir for use with reservoir REMD simulations using\n"
          "  energies in <energy data set>, temperature <temp0> and random seed <iseed>\n"
          "  Do not include velocities/forces if 'novelocity'/'noforce' specified.\n"
          "  If <cluster bin data set> is specified from e.g. a previous\n"
          "  'clusterdihedral' command, the reservoir can be used for non-Boltzmann\n"
          "  reservoir REMD (rremd==3).\n");
}

// Action_CreateReservoir::Init()
Action::RetType Action_CreateReservoir::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifndef BINTRAJ
  mprinterr("Error: NetCDF reservoir requires NetCDF support. Recompile with -DBINTRAJ.\n");
  return Action::ERR;
# endif
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'createreservoir' action does not work with > 1 thread (%i threads currently).\n", init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  // Get keywords
  filename_.SetFileName( actionArgs.GetStringNext() );
  if (filename_.empty()) {
    mprinterr("Error: createreservoir: No filename specified.\n");
    return Action::ERR;
  }
  reservoirT_ = actionArgs.getKeyDouble("temp0", -1.0);
  if (reservoirT_ < 0.0) {
    mprinterr("Error: Reservoir temperature must be specified and cannot be < 0.0\n");
    return Action::ERR;
  }
  iseed_ = actionArgs.getKeyInt("iseed", 0);
  if (iseed_ < 1) {
    mprinterr("Error: Reservoir random seed must be specified and > 0\n");
    return Action::ERR;
  }
  useVelocity_ = !actionArgs.hasKey("novelocity");
  useForce_ = !actionArgs.hasKey("noforce");
  // Get parm for reservoir traj
  original_trajparm_ = init.DSL().GetTopology( actionArgs );
  if (original_trajparm_ == 0) {
    mprinterr("Error: createreservoir: no topology.\n");
    return Action::ERR;
  }
  // Get energy data set
  std::string eneDsname = actionArgs.GetStringKey("ene");
  DataSet* dstmp = init.DSL().GetDataSet( eneDsname );
  if (dstmp == 0) {
    mprinterr("Error: could not get energy data set %s\n", eneDsname.c_str());
    return Action::ERR;
  }
  if (dstmp->Type() != DataSet::FLOAT &&
      dstmp->Type() != DataSet::DOUBLE &&
      dstmp->Type() != DataSet::XYMESH)
  {
    mprinterr("Error: energy data set %s must be type FLOAT, DOUBLE, or XYMESH.\n",
              dstmp->legend());
    return Action::ERR;
  }
  if (dstmp->Ndim() != 1) {
    mprinterr("Error: energy data set is not 1D (%u)\n", dstmp->Ndim());
    return Action::ERR;
  }
  ene_ = static_cast<DataSet_1D*>( dstmp );
  // Get bin data set
  std::string binDSname = actionArgs.GetStringKey("bin");
  if (!binDSname.empty()) {
    dstmp = init.DSL().GetDataSet( binDSname );
    if (dstmp == 0) {
      mprinterr("Error: could not get bin data set %s\n", binDSname.c_str());
      return Action::ERR;
    } else if (dstmp->Ndim() != 1) {
      mprinterr("Error: bin data set must be one dimensional.\n");
      return Action::ERR;
    }
    bin_ = static_cast<DataSet_1D*>( dstmp );
  }
  trajIsOpen_ = false;
  nframes_ = 0;
  // Setup output reservoir file
  reservoir_.SetDebug( debugIn );
  // Set title
  title_ = actionArgs.GetStringKey("title");
  if (title_.empty())
    title_.assign("Cpptraj Generated structure reservoir");

  mprintf("    CREATERESERVOIR: '%s', energy data '%s'", filename_.full(), ene_->legend());
  if (bin_ != 0)
    mprintf(", bin data '%s'", bin_->legend());
  mprintf("\n\tTitle: %s\n", title_.c_str());
  mprintf("\tReservoir temperature= %.2f, random seed= %i\n", reservoirT_, iseed_);
  if (useVelocity_)
    mprintf("\tVelocities will be written to reservoir if present.\n");
  else
    mprintf("\tVelocities will not be written to reservoir.\n");
  if (useForce_)
    mprintf("\tForces will be written to reservoir if present.\n");
  else
    mprintf("\tForces will not be written to reservoir.\n");
  mprintf("\tTopology: %s\n", original_trajparm_->c_str());
  return Action::OK;
}

// Action_CreateReservoir::Setup()
Action::RetType Action_CreateReservoir::Setup(ActionSetup& setup) {
  // Check that input parm matches current parm
  if (original_trajparm_->Pindex() != setup.Top().Pindex()) {
    mprintf("Info: createreservoir was set up for topology %s\n", original_trajparm_->c_str());
    mprintf("Info: skipping topology %s\n", setup.Top().c_str());
    return Action::SKIP;
  }
  if (!trajIsOpen_) {
    mprintf("\tCreating reservoir file %s\n", filename_.full());
    CoordinateInfo cinfo = setup.CoordInfo();
    if (!useVelocity_) cinfo.SetVelocity(false);
    if (!useForce_)    cinfo.SetForce(false);
    if (reservoir_.InitReservoir(filename_, title_, cinfo, setup.Top().Natom(),
                                 (bin_ != 0), reservoirT_, iseed_))
    {
      mprinterr("Error: Could not set up NetCDF reservoir.\n");
      return Action::ERR;
    }
    trajIsOpen_ = true;
    nframes_ = 0;
  }
  return Action::OK;
}

// Action_CreateReservoir::DoAction()
Action::RetType Action_CreateReservoir::DoAction(int frameNum, ActionFrame& frm) {
  int bin = -1;
  if (bin_ != 0) bin = (int)bin_->Dval(frm.TrajoutNum());
  if (reservoir_.WriteReservoir(nframes_++, frm.Frm(), ene_->Dval(frm.TrajoutNum()), bin))
    return Action::ERR;
  return Action::OK;
}

// Action_CreateReservoir::Print()
void Action_CreateReservoir::Print() {
  mprintf("\tReservoir %s: %u frames.\n", filename_.base(), nframes_);
  reservoir_.CloseReservoir();
  trajIsOpen_ = false;
}
