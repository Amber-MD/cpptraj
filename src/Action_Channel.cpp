#include <cmath> // pow
#include "Action_Channel.h"
#include "CpptrajStdio.h"
#include "DataSet_GridFlt.h"

// CONSTRUCTOR
Action_Channel::Action_Channel() : Action(HIDDEN),
  grid_(0), dxyz_(-1.0) {}

void Action_Channel::Help() const {
  mprintf("\t<solute mask> [<solvent mask>] [out <file>] [dx <dx> [dy <dy>] [dz <dz>]]\n");
}

Action::RetType Action_Channel::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Keywords.
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  dxyz_[0] = actionArgs.getKeyDouble("dx", 0.35);
  dxyz_[1] = actionArgs.getKeyDouble("dy", dxyz_[0]);
  dxyz_[2] = actionArgs.getKeyDouble("dz", dxyz_[1]);
  // solute mask
  std::string sMask = actionArgs.GetMaskNext();
  if (sMask.empty()) {
    mprinterr("Error: No solute mask specified.\n");
    return Action::ERR;
  }
  if (soluteMask_.SetMaskString( sMask )) return Action::ERR;
  // solvent mask
  sMask = actionArgs.GetMaskNext();
  if (sMask.empty())
    sMask.assign(":WAT@O");
  if (solventMask_.SetMaskString( sMask )) return Action::ERR;

  // Grid Data Set
  grid_ = init.DSL().AddSet(DataSet::GRID_FLT, actionArgs.GetStringNext(), "Channel");
  if (grid_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( grid_ );

  mprintf("Warning: *** THIS ACTION IS EXPERIMENTAL AND NOT FULLY IMPLEMENTED. ***\n");
  mprintf("    CHANNEL: Solute mask [%s], solvent mask [%s]\n",
          soluteMask_.MaskString(), solventMask_.MaskString());
  mprintf("\tSpacing: XYZ={ %g %g %g }\n", dxyz_[0], dxyz_[1], dxyz_[2]);
  return Action::OK;
}

Action::RetType Action_Channel::Setup(ActionSetup& setup) {
  // Initial grid setup
  if (grid_->Size() == 0) {
    DataSet_3D& GRID = static_cast<DataSet_3D&>( *grid_ );
    Box const& box = setup.CoordInfo().TrajBox();
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
  if (setup.Top().SetupIntegerMask( soluteMask_ ) ||
      setup.Top().SetupIntegerMask( solventMask_)   )
    return Action::ERR;
  soluteMask_.MaskInfo();
  if (soluteMask_.None()) {
    mprintf("Warning: No solute atoms selected.\n");
    return Action::SKIP;
  }
  solventMask_.MaskInfo();
  if (solventMask_.None()) {
    mprintf("Warning: No solvent atoms selected.\n");
    return Action::SKIP;
  }
  // Set up solute van der Waals. FIXME: Handle case where no LJ params
  radii_.clear();
  for (AtomMask::const_iterator uAtom = soluteMask_.begin();
                                uAtom != soluteMask_.end(); ++uAtom)
    radii_.push_back( setup.Top().GetVDWradius( *uAtom ) );
  return Action::OK;
}

Action::RetType Action_Channel::DoAction(int frameNum, ActionFrame& frm) {
  // TODO: Gridding should be a DataSet_3d routine.
  DataSet_GridFlt& GRID = static_cast<DataSet_GridFlt&>( *grid_ );
  const float SOLUTE = 1.0;
  const float BULK = 0.0;
  long int nx = (long int)GRID.NX();
  long int ny = (long int)GRID.NY();
  long int nz = (long int)GRID.NZ();
  // Reset grid
  for (DataSet_GridFlt::iterator gval = GRID.begin(); gval != GRID.end(); ++gval)
    *gval = BULK;
  // Loop over solute atoms
  std::vector<double>::const_iterator radius = radii_.begin();
  for (AtomMask::const_iterator uAtom = soluteMask_.begin();
                                uAtom != soluteMask_.end(); ++uAtom, ++radius)
  {
    // Super naive approach.
    //double r2 = (*radius) * (*radius);
    Vec3 pt(frm.Frm().XYZ(*uAtom));
    mprintf("\nAtom %i  radius= %g Ang.\n", *uAtom + 1, *radius);
    pt.Print("   coords");
    // Minimum and maxmimum indicies
    Vec3 minPt = pt - *radius;
    Vec3 maxPt = pt + *radius;
    minPt.Print("min point");
    maxPt.Print("max point");
    long int min_i, min_j, min_k;
    GRID.Bin().Indices(minPt[0],minPt[1],minPt[2],min_i,min_j,min_k);
    long int max_i, max_j, max_k;
    GRID.Bin().Indices(maxPt[0],maxPt[1],maxPt[2],max_i,max_j,max_k);
    mprintf("\tGrid dims: %li <= i < %li\n", std::max(0L,min_i), std::min(max_i,nx));
    mprintf("\tGrid dims: %li <= j < %li\n", std::max(0L,min_j), std::min(max_j,ny));
    mprintf("\tGrid dims: %li <= k < %li\n", std::max(0L,min_k), std::min(max_k,nz));
    // TODO: Check spherical distance.
    // TODO cast down to ints?
    for (long int i = std::max(0L,min_i); i <= std::min(max_i,nx); i++)
      for (long int j = std::max(0L,min_j); j <= std::min(max_j,ny); j++)
        for (long int k = std::max(0L,min_k); k <= std::min(max_k,nz); k++)
          GRID.SetElement(i,j,k,SOLUTE);
  }

  return Action::OK;
}
