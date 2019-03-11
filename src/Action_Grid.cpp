#include <cmath> // exp
#include <algorithm> // std::max_element
#include "Action_Grid.h"
#include "CpptrajStdio.h"
#include "PDBfile.h"
#include "MaskArray.h"

// CONSTRUCTOR
Action_Grid::Action_Grid() :
  normalize_(NONE),
  // Default particle density (molecules/Ang^3) for water based on 1.0 g/mL
  density_(0.033456),
  max_(0.80),
  madura_(0),
  smooth_(0),
  nframes_(0),
  debug_(0),
  invert_(false),
  pdbfile_(0),
  grid_(0)
{}

void Action_Grid::Help() const {
  mprintf("\t[out <filename>]\n%s\n", GridAction::HelpText);
  mprintf("\t<mask> [normframe | normdensity [density <density>]]\n"
          "\t[pdb <pdbout> [max <fraction>]] \n"
          "\t[[smoothdensity <value>] [invert]] [madura <madura>]\n"
          "  Bin atoms in <mask> into a 3D grid written to <filename>.\n");
}

// Action_Grid::Init()
Action::RetType Action_Grid::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  nframes_ = 0;
  // Get output filename
  std::string filename = actionArgs.GetStringKey("out");
  // Get grid options
  grid_ = GridInit( "GRID", actionArgs, init.DSL() );
  if (grid_ == 0) return Action::ERR;
# ifdef MPI
  if (ParallelGridInit(init.TrajComm(), grid_)) return Action::ERR;
# endif
  // Get extra options
  max_ = actionArgs.getKeyDouble("max", 0.80);
  madura_ = actionArgs.getKeyDouble("madura", 0);
  smooth_ = actionArgs.getKeyDouble("smoothdensity", 0);
  invert_ = actionArgs.hasKey("invert");
  pdbfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("pdb"),"Grid PDB",DataFileList::PDB,true);
  density_ = actionArgs.getKeyDouble("density",0.033456);
  if (actionArgs.hasKey("normframe")) normalize_ = TO_FRAME;
  else if (actionArgs.hasKey("normdensity")) normalize_ = TO_DENSITY;
  else normalize_ = NONE;
  if (normalize_ != NONE && (smooth_ > 0.0 || madura_ > 0.0)) {
    mprinterr("Error: Normalize options are not compatible with smoothdensity/madura options.\n");
    init.DSL().RemoveSet( grid_ );
    return Action::ERR;
  }
  // Get mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: GRID: No mask specified.\n");
    init.DSL().RemoveSet( grid_ );
    return Action::ERR;
  }
  mask_.SetMaskString(maskexpr);

  // Setup output file
  // For backwards compat., if no 'out' assume next string is filename
  if (filename.empty() && actionArgs.Nargs() > 1 && !actionArgs.Marked(1))
    filename = actionArgs.GetStringNext();
  DataFile* outfile = init.DFL().AddDataFile(filename, actionArgs);
  if (outfile != 0) outfile->AddDataSet((DataSet*)grid_);

  // Info
  mprintf("    GRID:\n");
  GridInfo( *grid_ );
  if (outfile != 0) mprintf("\tGrid will be printed to file %s\n", outfile->DataFilename().full());
  mprintf("\tGrid data set: '%s'\n", grid_->legend());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());
  if (pdbfile_ != 0)
      mprintf("\tPseudo-PDB will be printed to %s\n", pdbfile_->Filename().full());
  if (normalize_ == TO_FRAME)
    mprintf("\tGrid will be normalized by number of frames.\n");
  else if (normalize_ == TO_DENSITY)
    mprintf("\tGrid will be normalized to a density of %g molecules/Ang^3.\n", density_);
  // TODO: print extra options

  return Action::OK;
}

// Action_Grid::Setup()
Action::RetType Action_Grid::Setup(ActionSetup& setup) {
  // Setup grid, checks box info.
  if (GridSetup( setup.Top(), setup.CoordInfo() )) return Action::ERR;

  // Setup mask
  if (setup.Top().SetupIntegerMask( mask_ ))
    return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected for topology %s\n", setup.Top().c_str());
    return Action::SKIP;
  }

  // TEST
  Cpptraj::MaskArray mArray;
  mArray.SetType( Cpptraj::MaskArray::BY_RESIDUE );
  mArray.SetupMasks( mask_, setup.Top() );

  return Action::OK;
}

// Action_Grid::DoAction()
Action::RetType Action_Grid::DoAction(int frameNum, ActionFrame& frm) {
  GridFrame( frm.Frm(), mask_, *grid_ );
  ++nframes_;
  return Action::OK;
}

// Action_Grid::print()
void Action_Grid::Print() {
  if (nframes_ < 1) return;
  // Perform normalization and find max.
  double gridMax = 0.0;
  if (normalize_ == NONE) {
    mprintf("    GRID: No normalization");
    if (smooth_ > 0.0) mprintf(", smoothing factor (%g)", smooth_);
    mprintf(".\n");
    for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
      double gridval = (double)(*gval);
      // ----- SMOOTHING -----
      if (smooth_ > 0.0) {
        double yy = gridval - smooth_;
        double xx = yy*yy / (0.2 * smooth_ * smooth_);
        xx = exp( -xx );
        if (invert_) {
          if (gridval > smooth_) // NOTE: Comparison OK? Needs cast?
            gridval = -5.0;
          else
            gridval -= gridval * xx;
          /* COMMENTED OUT IN ORIGINAL PTRAJ CODE
          if (gridInfo->grid[index] < action->darg3) {
            gridInfo->grid[index] = 0.0;
          }
          */
          if (gridval >= 0)
            gridval = smooth_ - gridval;
        } else {
          if (gridval < smooth_)
            gridval = 0;
          else
            gridval -= gridval * xx;
          if (gridval < smooth_)
            gridval = 0;
        }
      }
      // do the madura negative option to expose low density
      if ( madura_ > 0.0 && gridval > 0.0 && gridval < madura_ )
        *gval = (float)-gridval;
      else
        *gval = (float) gridval;
      if ( gridval > gridMax )
        gridMax = gridval;
    }
  } else {
    // Normalize to frames / density.
    mprintf("    GRID: Normalization");
    double dens = 1.0;
    if (normalize_ == TO_DENSITY) {
      dens = grid_->Bin().VoxelVolume() * density_;
      mprintf(" to density %g molecules/Ang^3, voxel volume= %g Ang^3, %g mols/voxel,",
              density_, grid_->Bin().VoxelVolume(), dens);
    } else
      mprintf(" to");
    mprintf(" number of frames %u", nframes_);
    double norm = 1.0 / ((double)nframes_ * dens);
    mprintf(", normalization factor= %g\n",norm);
    for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
      double gridval = (double)(*gval) * norm;
      gridMax = std::max(gridval, gridMax);
      *gval = (float)gridval;
    }
  }
  mprintf("\tGrid max is %g\n", gridMax);
  // PDBfile output
  PrintPDB( gridMax );
}

// Action_Grid::PrintPDB()
void Action_Grid::PrintPDB(double gridMax)
{
  if (gridMax == 0.0) {
    mprinterr("Error: Grid max is 0. No density for PDB write.\n");
    return;
  }
  double norm = 1.0 / gridMax;
  // Calculate normalization if necessary
//  if (norm < 0.0) {
//    norm = (double)*std::max_element(grid_->begin(), grid_->end());
//  }
  // Write PDB
  PDBfile& pdbout = static_cast<PDBfile&>( *pdbfile_ );
  mprintf("\tWriting PDB of grid points > %.2f%% of grid max.\n", max_*100.0);
  int res = 1;
  for (size_t k = 0; k < grid_->NZ(); ++k) {
    for (size_t j = 0; j < grid_->NY(); ++j) {
      for (size_t i = 0; i < grid_->NX(); ++i) {
        double gridval = grid_->GetElement(i, j, k) * norm;
        if (gridval > max_) {
          Vec3 cxyz = grid_->Bin().Center(i,j,k);
          pdbout.WriteATOM(res++, cxyz[0], cxyz[1], cxyz[2], "GRID", gridval);
        }
      }
    }
  }
  // Write grid boundaries
  for (size_t k = 0; k <= grid_->NZ(); k += grid_->NZ())
    for (size_t j = 0; j <= grid_->NY(); j += grid_->NY())
      for (size_t i = 0; i <= grid_->NX(); i += grid_->NX()) {
        Vec3 cxyz = grid_->Bin().Corner(i,j,k);
        pdbout.WriteHET(res, cxyz[0], cxyz[1], cxyz[2]);
      }
  // DEBUG: Write all grid bin corners
  if (debug_ > 1) {
    ++res;
    for (size_t k = 0; k <= grid_->NZ(); k++)
      for (size_t j = 0; j <= grid_->NY(); j++)
        for (size_t i = 0; i <= grid_->NX(); i++) {
          Vec3 cxyz = grid_->Bin().Corner(i,j,k);
          pdbout.WriteATOM(res, cxyz[0], cxyz[1], cxyz[2], "BIN", 0.0);
        }
  }
}
