#include "Action_GIST.h"
#include "CpptrajStdio.h"

Action_GIST::Action_GIST() :
  gO_(0),
  gH_(0),
  Esw_(0),
  Eww_(0),
  dTStrans_(0),
  dTSorient_(0),
  dTSsix_(0),
  neighbor_norm_(0),
  dipole_(0),
  order_norm_(0),
  dipolex_(0),
  dipoley_(0),
  dipolez_(0),
  datafile_(0),
  BULK_DENS_(0.0),
  temperature_(0.0),
  NFRAME_(0),
  doOrder_(false),
  doEij_(false),
  skipE_(false)
{}

void Action_GIST::Help() const {
  mprintf("\t[doorder] [doeij] [skipE] [refdens <rdval>] [Temp <tval>]\n"
          "\t[gridcntr <xval> <yval> <zval>]\n"
          "\t[griddim <xval> <yval> <zval>] [gridspacn <spaceval>]\n"
          "\t[out <filename>] [noimage]\n");
}

Action::RetType Action_GIST::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  std::string gistout = actionArgs.GetStringKey("out");
  if (gistout.empty()) gistout.assign("gist-output.dat");
  datafile_ = init.DFL().AddCpptrajFile( gistout, "GIST output" );
  if (datafile_ == 0) return Action::ERR;
  doOrder_ = actionArgs.hasKey("doorder");
  doEij_ = actionArgs.hasKey("doeij");
  skipE_ = actionArgs.hasKey("skipE");
  // Set Bulk Density 55.5M
  BULK_DENS_ = actionArgs.getKeyDouble("refdens", 0.0334);
  if ( BULK_DENS_ > (0.0334*1.2) )
    mprintf("Warning: water reference density is high, consider using 0.0334 for 1g/cc water density\n");
  else if ( BULK_DENS_ < (0.0334*0.8) )
    mprintf("Warning: water reference density is low, consider using 0.0334 for 1g/cc water density\n");
  temperature_ = actionArgs.getKeyDouble("temp", 300.0);
  if (temperature_ < 0.0) {
    mprinterr("Error: Negative temperature specified.\n");
    return Action::ERR;
  }
  // Grid spacing
  double gridspacing = actionArgs.getKeyDouble("gridspacn", 0.50);
  // Grid center
  Vec3 gridcntr(0.0);
  if ( actionArgs.hasKey("gridcntr") ) {
    gridcntr[0] = actionArgs.getNextDouble(-1);
    gridcntr[1] = actionArgs.getNextDouble(-1);
    gridcntr[2] = actionArgs.getNextDouble(-1);
  } else
    mprintf("Warning: No grid center values specified, using default (origin)\n");
  // Grid dimensions
  int nx = 40;
  int ny = 40;
  int nz = 40;
  if ( actionArgs.hasKey("griddim") ) {
    nx = actionArgs.getNextInteger(-1);
    ny = actionArgs.getNextInteger(-1);
    nz = actionArgs.getNextInteger(-1);
  } else
    mprintf("Warning: No grid dimension values specified, using default (40,40,40)\n");
  // Data set name
  std::string dsname = actionArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = init.DSL().GenerateDefaultName("GIST");

  // Set up grid data sets
  gO_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "gO"));
  gH_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "gH"));
  Esw_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "Esw"));
  Eww_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "Eww"));
  dTStrans_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "dTStrans"));
  dTSorient_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "dTSorient"));
  dTSsix_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "dTSsix"));
  neighbor_norm_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "neighbor"));
  dipole_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname, "dipole"));

  order_norm_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_DBL, MetaData(dsname, "order"));
  dipolex_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_DBL, MetaData(dsname, "dipolex"));
  dipoley_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_DBL, MetaData(dsname, "dipoley"));
  dipolez_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_DBL, MetaData(dsname, "dipolez"));
 
  // Set up grids. TODO non-orthogonal as well
  Vec3 v_spacing( gridspacing );
  gO_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  gH_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  Esw_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  Eww_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dTStrans_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dTSorient_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dTSsix_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  neighbor_norm_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dipole_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);

  order_norm_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dipolex_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dipoley_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);
  dipolez_->Allocate_N_C_D(nx, ny, nz, gridcntr, v_spacing);

  // Set up grid params
  Box gbox;
  gbox.SetBetaLengths( 90.0, (double)nx * gridspacing,
                             (double)ny * gridspacing,
                             (double)nz * gridspacing );
  grid_.Setup_O_Box( nx, ny, nz, gO_->GridOrigin(), gbox );

  mprintf("    GIST:\n");
  if(doOrder_)
    mprintf("\tDo Order calculation\n");
  else
    mprintf("\tSkip Order calculation\n");
  if(doEij_)
    mprintf("\tCompute and print water-water Eij matrix\n");
  else
    mprintf("\tSkip water-water Eij matrix\n");
  mprintf("\tWater reference density: %6.4f\n", BULK_DENS_); // TODO units
  mprintf("\tSimulation temperature: %6.4f K\n", temperature_);
  if (image_.UseImage())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tDistances will not be imaged.\n");
  gO_->GridInfo();
  mprintf("\t#Please cite these papers if you use GIST results in a publication:\n"
          "\t#    Steven Ramsey, Crystal Nguyen, Romelia Salomon-Ferrer, Ross C. Walker, Michael K. Gilson, and Tom Kurtzman J. Comp. Chem. 37 (21) 2016\n"
          "\t#    Crystal Nguyen, Michael K. Gilson, and Tom Young, arXiv:1108.4876v1 (2011)\n"
          "\t#    Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson,\n"
          "\t#      J. Chem. Phys. 137, 044101 (2012)\n"
          "\t#    Lazaridis, J. Phys. Chem. B 102, 3531–3541 (1998)\n");

  return Action::OK;
}

// Action_GIST::Setup()
Action::RetType Action_GIST::Setup(ActionSetup& setup) {
  // We need box info
  if (setup.CoordInfo().TrajBox().Type() == Box::NOBOX) {
    mprinterr("Error: Must have explicit solvent with periodic boundaries!");
    return Action::ERR;
  }
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );

  // Get molecule number for each solvent molecule
  mol_nums_.clear();
  unsigned int midx = 0;
  for (Topology::mol_iterator mol = setup.Top().MolStart();
                              mol != setup.Top().MolEnd(); ++mol, ++midx)
  {
    if (mol->IsSolvent())
      mol_nums_.push_back( midx );
  }

  water_voxel_.assign( mol_nums_.size(), -1 );

  return Action::OK;
}

// Action_GIST::DoAction()
Action::RetType Action_GIST::DoAction(int frameNum, ActionFrame& frm) {
  NFRAME_++;

  // Loop over each solvent molecule
  for (Iarray::const_iterator smol = mol_nums
