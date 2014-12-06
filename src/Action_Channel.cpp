#include "Action_Channel.h"
#include "CpptrajStdio.h"
#include "DataSet_GridFlt.h"

void Action_Channel::Help() {
  mprintf("\t<solute mask> [<solvent mask>] [out <file>] [dx <dx> [dy <dy>] [dz <dz>]]\n");
}

Action::RetType Action_Channel::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Keywords.
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  dxyz_[0] = actionArgs.getKeyDouble("dx", 0.35);
  dxyz_[1] = actionArgs.getKeyDouble("dy", dxyz_[0]);
  dxyz_[2] = actionArgs.getKeyDouble("dz", dxyz_[1]);
  // solute mask
  std::string sMask = actionArgs.GetMaskNext();
  if (sMask.empty()) {
    mprinterr("Error: No solute mask specified.\n");
    return Action::ERR;
  }
  soluteMask_.SetMaskString( sMask );
  // solvent mask
  sMask = actionArgs.GetMaskNext();
  if (sMask.empty())
    sMask.assign(":WAT@O");
  solventMask_.SetMaskString( sMask );

  // Grid Data Set
  grid_ = DSL->AddSet(DataSet::GRID_FLT, actionArgs.GetStringNext(), "Channel");
  if (grid_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddSet( grid_ );

  mprintf("    CHANNEL: Solute mask [%s], solvent mask [%s]\n",
          soluteMask_.MaskString(), solventMask_.MaskString());
  mprintf("\tSpacing: XYZ={ %g %g %g }\n", dxyz_[0], dxyz_[1], dxyz_[2]);
  return Action::OK;
}

Action::RetType Action_Channel::Setup(Topology* currentParm, Topology** parmAddress) {
  // Initial grid setup
  if (grid_->Size() == 0) {
    DataSet_3D& GRID = static_cast<DataSet_3D&>( *grid_ );
    Box const& box = currentParm->ParmBox();
    if (box.Type() == Box::NOBOX) {
      mprinterr("Error: No box information to set up grid.\n");
      return Action::ERR;
    } else if (box.Type() == Box::ORTHO) {
      // FIXME: May need to update parm box info or set up on first frame.
      if (GRID.Allocate_X_C_D(box.Lengths(), box.Center(), dxyz_)) return Action::ERR; 
    } else {
      Vec3 nxyz = box.Lengths() / dxyz_;
      if (GRID.Allocate_N_O_Box((size_t)nxyz[0], (size_t)nxyz[1], (size_t)nxyz[2],
                                Vec3(0.0), box)) return Action::ERR;
    }
    GRID.GridInfo();
  }
  // Set up masks
  if (currentParm->SetupIntegerMask( soluteMask_ ) ||
      currentParm->SetupIntegerMask( solventMask_)   )
    return Action::ERR;
  soluteMask_.MaskInfo();
  if (soluteMask_.None()) {
    mprinterr("Error: channel: No solute atoms selected.\n");
    return Action::ERR;
  }
  solventMask_.MaskInfo();
  if (solventMask_.None()) {
    mprinterr("Error: channel: No solvent atoms selected.\n");
    return Action::ERR;
  }
  // Set up solute van der Waals. FIXME: Handle case where no LJ params
  const double one_over_6 = 1.0 / 6.0;
  radii_.clear();
  for (AtomMask::const_iterator uAtom = soluteMask_.begin();
                                uAtom != soluteMask_.end(); ++uAtom)
  {
    NonbondType const& LJ = currentParm->GetLJparam(*uAtom, *uAtom);
    radii_.push_back( 0.5 * pow(2 * LJ.A() / LJ.B(), one_over_6);
  }
  return Action::OK;
}

Action::RetType Action_Channel::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress)
{
  // TODO: Gridding should be a DataSet_3d routine.
  DataSet_GridFlt& GRID = static_cast<DataSet_GridFlt&>( *grid_ ); 
  const float SOLUTE = 1.0;
  std::vector<double>::const_iterator radius = radii_.begin();
  for (AtomMask::const_iterator uAtom = soluteMask_.begin();
                                uAtom != soluteMask_.end(); ++uAtom, ++radius)
  {
    // Super naive approach. 
