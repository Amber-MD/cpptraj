#include <cfloat>
#include <cmath> // ceil
#include "Action_Bounds.h"
#include "CpptrajStdio.h"
#include "DataSet_3D.h" // For allocating grid

// CONSTRUCTOR
Action_Bounds::Action_Bounds() :
  outfile_(0),
  offset_(1),
  grid_(0),
  ds_xmin_(0),
  ds_ymin_(0),
  ds_zmin_(0),
  ds_xmax_(0),
  ds_ymax_(0),
  ds_zmax_(0)
{}

void Action_Bounds::Help() const {
  mprintf("\t[<mask>] [out <filename>] [name <setname>]\n"
          "\t[dx <dx> [dy <dy>] [dz <dz>] [offset <offset>]]\n"
          "  Calcuate the max/min coordinates (X,Y,Z) of atoms in <mask>.\n"
          "    <mask>         : Atoms to calculate boundaries for.\n"
          "    out <filename> : Write boundaries to <filename>.\n"
          "    name <setname> : Output data set name.\n"
          "  If dx/dy/dz are specified, a grid will be created from the calculated\n"
          "  boundaries after processing is complete. In this case 'name' must be\n"
          "  specified.\n"
          "    dx <dx>          : Grid spacing in the X direction.\n"
          "    dy <dy>          : Grid spacing in the Y direction. Use <dx> if not given.\n"
          "    dz <dz>          : Grid spacing in the Z direction. Use <dy> if not given.\n"
          "    offset <offset>] : Number of bins to add in each direction to grid.\n");
}

// Action_Bounds::Init()
Action::RetType Action_Bounds::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  outfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("out"), "Bounds",
                                 DataFileList::TEXT, true);
  dxyz_[0] = actionArgs.getKeyDouble("dx", -1.0);
  dxyz_[1] = actionArgs.getKeyDouble("dy", -1.0);
  dxyz_[2] = actionArgs.getKeyDouble("dz", -1.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  std::string dsname = actionArgs.GetStringKey("name");
  offset_ = actionArgs.getKeyInt("offset", 1);
  if (dxyz_[0] > -1.0) {
    // Set up grid
    if (dsname.empty()) {
      mprinterr("Error: Grid name must be specified if spacing specified.\n");
      return Action::ERR;
    } 
    // Set default spacings if necessary. Y uses X, Z uses Y.
    if (dxyz_[1] < 0.0) dxyz_[1] = dxyz_[0];
    if (dxyz_[2] < 0.0) dxyz_[2] = dxyz_[1];
    grid_ = init.DSL().AddSet(DataSet::GRID_FLT, dsname, "Bounds");
    if (grid_ == 0) return Action::ERR;
  } else {
    if (dsname.empty())
      dsname = init.DSL().GenerateDefaultName("BOUNDS");
  }
  // Set up data sets
  ds_xmin_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "xmin"));
  ds_ymin_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "ymin"));
  ds_zmin_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "zmin"));
  ds_xmax_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "xmax"));
  ds_ymax_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "ymax"));
  ds_zmax_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "zmax"));
  if (ds_xmin_ == 0 || ds_ymin_ == 0 || ds_zmin_ == 0 ||
      ds_xmax_ == 0 || ds_ymax_ == 0 || ds_zmax_ == 0)
    return Action::ERR;

  min_[0] = DBL_MAX;
  min_[1] = min_[0];
  min_[2] = min_[0];
  max_[0] = -DBL_MAX;
  max_[1] = max_[0];
  max_[2] = max_[0];

  mprintf("    BOUNDS: Calculating bounds for atoms in mask [%s]\n", mask_.MaskString());
  mprintf("\tOutput to '%s'\n", outfile_->Filename().full());
  if (grid_ != 0) {
    mprintf("\tGrid %s will be created after processing using\n"
            "\t  spacings dX= %g  dY= %g  dZ= %g  offset= %i Bins.\n",
            grid_->legend(), dxyz_[0], dxyz_[1], dxyz_[2], offset_);
  }
# ifdef MPI
  trajComm_ = init.TrajComm();
  // Since grid is not allocated until Print(), no sync needed.
  if (grid_ != 0) grid_->SetNeedsSync(false);
# endif
  return Action::OK;
}

// Action_Bounds::Setup()
Action::RetType Action_Bounds::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: bounds: No atoms selected in mask.\n");
    return Action::SKIP;
  }
  return Action::OK;
}

// Action_Bounds::DoAction()
Action::RetType Action_Bounds::DoAction(int frameNum, ActionFrame& frm) {
  for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
  {
    const double* xyz = frm.Frm().XYZ( *atom );
    if (xyz[0] < min_[0]) min_[0] = xyz[0];
    if (xyz[0] > max_[0]) max_[0] = xyz[0];
    if (xyz[1] < min_[1]) min_[1] = xyz[1];
    if (xyz[1] > max_[1]) max_[1] = xyz[1];
    if (xyz[2] < min_[2]) min_[2] = xyz[2];
    if (xyz[2] > max_[2]) max_[2] = xyz[2];
  }
  return Action::OK;
}

#ifdef MPI
int Action_Bounds::SyncAction() {
  double buf[3];
  trajComm_.Reduce( &buf, min_, 3, MPI_DOUBLE, MPI_MIN );
  if (trajComm_.Master())
    std::copy( buf, buf+3, min_ );
  trajComm_.Reduce( &buf, max_, 3, MPI_DOUBLE, MPI_MAX );
  if (trajComm_.Master())
    std::copy( buf, buf+3, max_ );
  return 0;
}
#endif

void Action_Bounds::Print() {
  static const char cXYZ[3] = {'X', 'Y', 'Z'};
  Vec3 center;
  size_t nxyz[3];

  mprintf("    BOUNDS: Output to %s\n", outfile_->Filename().full());
  ds_xmin_->Add(0, min_  );
  ds_ymin_->Add(0, min_+1);
  ds_zmin_->Add(0, min_+2);
  ds_xmax_->Add(0, max_  );
  ds_ymax_->Add(0, max_+1);
  ds_zmax_->Add(0, max_+2);
  for (int i = 0; i < 3; i++) {
    outfile_->Printf("%f < %c < %f", min_[i], cXYZ[i], max_[i]);
    if (dxyz_[i] > 0.0) {
      center[i] = (max_[i] + min_[i]) / 2.0;
      long int nbins = (long int)ceil( (max_[i] - min_[i]) / dxyz_[i] ) + (long int)offset_;
      nxyz[i] = (size_t)nbins;
      outfile_->Printf("\tCenter= %f  Bins=%zu", center[i], nxyz[i]);
    }
    outfile_->Printf("\n");
  }
  if (grid_ != 0) {
    DataSet_3D& grid3d = static_cast<DataSet_3D&>( *grid_ );
    if (grid3d.Allocate_N_C_D( nxyz[0], nxyz[1], nxyz[2], center, dxyz_ ))
      mprinterr("Error: Could not allocate grid %s\n", grid_->legend());
  }
}
