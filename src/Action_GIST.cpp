#include <cmath>
#include <cfloat> // DBL_MAX
#include "Action_GIST.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_MatrixFlt.h"
#include "DataSet_GridFlt.h"
#include "DataSet_GridDbl.h"
#include "ProgressBar.h"
#include "StringRoutines.h"
#include "DistRoutines.h"
#ifdef _OPENMP
# include <omp.h>
#endif

const double Action_GIST::maxD_ = DBL_MAX;

Action_GIST::Action_GIST() :
  debug_(0),
  numthreads_(1),
#ifdef CUDA
  numberAtoms_(0),
  numberAtomTypes_(0),
  headAtomType_(0),
  solvent_(NULL),
  NBindex_c_(NULL),
  molecule_c_(NULL),
  paramsLJ_c_(NULL),
  max_c_(NULL),
  min_c_(NULL),
  result_w_c_(NULL),
  result_s_c_(NULL),
  result_O_c_(NULL),
  result_N_c_(NULL),
#endif
  gridspacing_(0),
  gridcntr_(0.0),
  griddim_(0.0),
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
  PME_(0),
  U_PME_(0),
  ww_Eij_(0),
  G_max_(0.0),
  CurrentParm_(0),
  datafile_(0),
  eijfile_(0),
  infofile_(0),
  fltFmt_(TextFormat::GDOUBLE),
  intFmt_(TextFormat::INTEGER),
  BULK_DENS_(0.0),
  temperature_(0.0),
  NeighborCut2_(12.25), // 3.5^2
//  system_potential_energy_(0),
//  solute_potential_energy_(0),
  MAX_GRID_PT_(0),
  NSOLVENT_(0),
  N_ON_GRID_(0),
  nMolAtoms_(0),
  NFRAME_(0),
  max_nwat_(0),
  doOrder_(false),
  doEij_(false),
  skipE_(false),
  includeIons_(true)
{}

/** GIST help */
void Action_GIST::Help() const {
  mprintf("\t[doorder] [doeij] [skipE] [skipS] [refdens <rdval>] [temp <tval>]\n"
          "\t[noimage] [gridcntr <xval> <yval> <zval>] [excludeions]\n"
          "\t[griddim <nx> <ny> <nz>] [gridspacn <spaceval>] [neighborcut <ncut>]\n"
          "\t[prefix <filename prefix>] [ext <grid extension>] [out <output suffix>]\n"
          "\t[floatfmt {double|scientific|general}] [floatwidth <fw>] [floatprec <fp>]\n"
          "\t[intwidth <iw>]\n"
          "\t[info <info suffix>]\n");
#         ifdef LIBPME
          mprintf("\t[nopme|pme %s\n\t %s\n\t %s]\n", EwaldOptions::KeywordsCommon1(), EwaldOptions::KeywordsCommon2(), EwaldOptions::KeywordsPME());
#         endif
          mprintf("Perform Grid Inhomogenous Solvation Theory calculation.\n"
#ifdef CUDA
          "The option doeij is not available, when using the CUDA accelerated version,\n"
          "as this would need way too much memory."
#endif
          );
}

/** Init GIST action. */
Action::RetType Action_GIST::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'gist' action does not work with > 1 process (%i processes currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  gist_init_.Start();
  prefix_ = actionArgs.GetStringKey("prefix");
  if (prefix_.empty()) prefix_.assign("gist");
  std::string ext = actionArgs.GetStringKey("ext");
  if (ext.empty()) ext.assign(".dx");
  std::string gistout = actionArgs.GetStringKey("out");
  if (gistout.empty()) gistout.assign(prefix_ + "-output.dat");
  datafile_ = init.DFL().AddCpptrajFile( gistout, "GIST output" );
  if (datafile_ == 0) return Action::ERR;
  // Info file: if not specified use STDOUT
  gistout = actionArgs.GetStringKey("info");
  if (!gistout.empty()) gistout = prefix_ + "-" + gistout;
  infofile_ = init.DFL().AddCpptrajFile( gistout, "GIST info", DataFileList::TEXT, true );
  if (infofile_ == 0) return Action::ERR;
  // Grid files
  DataFile* file_gO = init.DFL().AddDataFile( prefix_ + "-gO" + ext );
  DataFile* file_gH = init.DFL().AddDataFile( prefix_ + "-gH" + ext );
  DataFile* file_Esw = init.DFL().AddDataFile(prefix_ + "-Esw-dens" + ext);
  DataFile* file_Eww = init.DFL().AddDataFile(prefix_ + "-Eww-dens" + ext);
  DataFile* file_dTStrans = init.DFL().AddDataFile(prefix_ + "-dTStrans-dens" + ext);
  DataFile* file_dTSorient = init.DFL().AddDataFile(prefix_ + "-dTSorient-dens" + ext);
  DataFile* file_dTSsix = init.DFL().AddDataFile(prefix_ + "-dTSsix-dens" + ext);
  DataFile* file_neighbor_norm = init.DFL().AddDataFile(prefix_ + "-neighbor-norm" + ext);
  DataFile* file_dipole = init.DFL().AddDataFile(prefix_ + "-dipole-dens" + ext);
  DataFile* file_order_norm = init.DFL().AddDataFile(prefix_ + "-order-norm" + ext);
  DataFile* file_dipolex = init.DFL().AddDataFile(prefix_ + "-dipolex-dens" + ext);
  DataFile* file_dipoley = init.DFL().AddDataFile(prefix_ + "-dipoley-dens" + ext);
  DataFile* file_dipolez = init.DFL().AddDataFile(prefix_ + "-dipolez-dens" + ext);
  // Output format keywords
  std::string floatfmt = actionArgs.GetStringKey("floatfmt");
  if (!floatfmt.empty()) {
    if (floatfmt == "double")
      fltFmt_.SetFormatType(TextFormat::DOUBLE);
    else if (floatfmt == "scientific")
      fltFmt_.SetFormatType(TextFormat::SCIENTIFIC);
    else if (floatfmt == "general")
      fltFmt_.SetFormatType(TextFormat::GDOUBLE);
    else {
      mprinterr("Error: Unrecognized format type for 'floatfmt': %s\n", floatfmt.c_str());
      return Action::ERR;
    }
  }
  fltFmt_.SetFormatWidthPrecision( actionArgs.getKeyInt("floatwidth", 0),
                                   actionArgs.getKeyInt("floatprec", -1) );
  intFmt_.SetFormatWidth( actionArgs.getKeyInt("intwidth", 0) );
  // Other keywords
  double neighborCut = actionArgs.getKeyDouble("neighborcut", 3.5);
  NeighborCut2_ = neighborCut * neighborCut;
  includeIons_ = !actionArgs.hasKey("excludeions");
  imageOpt_.InitImaging( !(actionArgs.hasKey("noimage")), actionArgs.hasKey("nonortho") );
  doOrder_ = actionArgs.hasKey("doorder");
  doEij_ = actionArgs.hasKey("doeij");
#ifdef CUDA
  if (this->doEij_) {
    mprinterr("Error: 'doeij' cannot be specified when using CUDA.\n");
    return Action::ERR;
  }
#endif
  skipE_ = actionArgs.hasKey("skipE");
  if (skipE_) {
    if (doEij_) {
      mprinterr("Error: 'doeij' cannot be specified if 'skipE' is specified.\n");
      return Action::ERR;
    }
  }
  // Parse PME options
  // TODO once PME output is stable, make pme true the default when LIBPME present.
//# ifdef LIBPME
//  usePme_ = true;
//# else
  usePme_ = false;
//# endif
# ifdef CUDA
  // Disable PME for CUDA
  usePme_ = false;
# endif
  if (actionArgs.hasKey("pme"))
    usePme_ = true;
  else if (actionArgs.hasKey("nopme"))
    usePme_ = false;
  // PME and doeij are not compatible
  if (usePme_ && doEij_) {
    mprinterr("Error: 'doeij' cannot be used with PME. Specify 'nopme' to use 'doeij'\n");
    return Action::ERR;
  }
  if (usePme_) {
#   ifdef LIBPME
    pmeOpts_.AllowLjPme(false);
    if (pmeOpts_.GetOptions(EwaldOptions::PME, actionArgs, "GIST")) {
      mprinterr("Error: Getting PME options for GIST failed.\n");
      return Action::ERR;
    }
#   else
    mprinterr("Error: 'pme' with GIST requires compilation with LIBPME.\n");
    return Action::ERR;
#   endif
  }
  DataFile* file_energy_pme = 0;
  DataFile* file_U_energy_pme = 0;
  if (usePme_) {
    file_energy_pme = init.DFL().AddDataFile(prefix_ + "-Water-Etot-pme-dens" + ext);
    file_U_energy_pme = init.DFL().AddDataFile(prefix_ + "-Solute-Etot-pme-dens"+ ext);
  }

  this->skipS_ = actionArgs.hasKey("skipS");

  if (doEij_) {
    eijfile_ = init.DFL().AddCpptrajFile(prefix_ + "-Eww_ij.dat", "GIST Eij matrix file");
    if (eijfile_ == 0) return Action::ERR;
  }
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
  gridspacing_ = actionArgs.getKeyDouble("gridspacn", 0.50);
  // Grid center
  gridcntr_ = Vec3(0.0);
  if ( actionArgs.hasKey("gridcntr") ) {
    gridcntr_[0] = actionArgs.getNextDouble(-1);
    gridcntr_[1] = actionArgs.getNextDouble(-1);
    gridcntr_[2] = actionArgs.getNextDouble(-1);
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
  griddim_ = Vec3((double)nx, (double)ny, (double)nz);
  // Data set name
  std::string dsname = actionArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = init.DSL().GenerateDefaultName("GIST");

  // Set up DataSets.
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

  if (gO_==0 || gH_==0 || Esw_==0 || Eww_==0 || dTStrans_==0 || dTSorient_==0 ||
      dTSsix_==0 || neighbor_norm_==0 || dipole_==0 || order_norm_==0 ||
      dipolex_==0 || dipoley_==0 || dipolez_==0)
    return Action::ERR;

  if (usePme_) {
    PME_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT, MetaData(dsname,"PME"));
    U_PME_ = (DataSet_3D*)init.DSL().AddSet(DataSet::GRID_FLT,MetaData(dsname,"U_PME"));
    if (PME_ == 0 || U_PME_ == 0) return Action::ERR;
  }

  if (doEij_) {
    ww_Eij_ = (DataSet_MatrixFlt*)init.DSL().AddSet(DataSet::MATRIX_FLT, MetaData(dsname, "Eij"));
    if (ww_Eij_ == 0) return Action::ERR;
  }

  // Allocate DataSets. TODO non-orthogonal grids as well
  Vec3 v_spacing( gridspacing_ );
  gO_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  MAX_GRID_PT_ = gO_->Size();
  gH_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  Esw_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  Eww_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dTStrans_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dTSorient_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dTSsix_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  neighbor_norm_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dipole_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);

  order_norm_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dipolex_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dipoley_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);
  dipolez_->Allocate_N_C_D(nx, ny, nz, gridcntr_, v_spacing);

  if (usePme_) {
    PME_->Allocate_N_C_D(nx,ny,nz,gridcntr_,v_spacing);
    U_PME_->Allocate_N_C_D(nx,ny,nz,gridcntr_,v_spacing);
  }

  if (ww_Eij_ != 0) {
    if (ww_Eij_->AllocateTriangle( MAX_GRID_PT_ )) {
      mprinterr("Error: Could not allocate memory for water-water Eij matrix.\n");
      return Action::ERR;
    }
  }

  // Add sets to files
  file_gO->AddDataSet( gO_ );
  file_gH->AddDataSet( gH_ );
  file_Esw->AddDataSet( Esw_ );
  file_Eww->AddDataSet( Eww_ );
  file_dTStrans->AddDataSet( dTStrans_ );
  file_dTSorient->AddDataSet( dTSorient_ );
  file_dTSsix->AddDataSet( dTSsix_ );
  file_neighbor_norm->AddDataSet( neighbor_norm_ );
  file_dipole->AddDataSet( dipole_ );
  file_order_norm->AddDataSet( order_norm_ );
  file_dipolex->AddDataSet( dipolex_ );
  file_dipoley->AddDataSet( dipoley_ );
  file_dipolez->AddDataSet( dipolez_ );
  if (usePme_) {
    file_energy_pme->AddDataSet(PME_);
    file_U_energy_pme->AddDataSet(U_PME_);
  }
  // Set up grid params TODO non-orthogonal as well
  G_max_ = Vec3( (double)nx * gridspacing_ + 1.5,
                 (double)ny * gridspacing_ + 1.5,
                 (double)nz * gridspacing_ + 1.5 );
  N_waters_.assign( MAX_GRID_PT_, 0 );
  N_solute_atoms_.assign( MAX_GRID_PT_, 0);
  N_hydrogens_.assign( MAX_GRID_PT_, 0 );
  voxel_xyz_.resize( MAX_GRID_PT_ ); // [] = X Y Z
  voxel_Q_.resize( MAX_GRID_PT_ ); // [] = W4 X4 Y4 Z4

  numthreads_ = 1;
# ifdef _OPENMP
# pragma omp parallel
  {
  if (omp_get_thread_num() == 0)
    numthreads_ = omp_get_num_threads();
  }
# endif

  if (!skipE_) {
    E_UV_VDW_.resize( numthreads_ );
    E_UV_Elec_.resize( numthreads_ );
    E_VV_VDW_.resize( numthreads_ );
    E_VV_Elec_.resize( numthreads_ );
    neighbor_.resize( numthreads_ );
    for (int thread = 0; thread != numthreads_; thread++) {
      E_UV_VDW_[thread].assign( MAX_GRID_PT_, 0 );
      E_UV_Elec_[thread].assign( MAX_GRID_PT_, 0 );
      E_VV_VDW_[thread].assign( MAX_GRID_PT_, 0 );
      E_VV_Elec_[thread].assign( MAX_GRID_PT_, 0 );
      neighbor_[thread].assign( MAX_GRID_PT_, 0 );
    }
    if (usePme_) {
      E_pme_.assign( MAX_GRID_PT_, 0 );
      U_E_pme_.assign( MAX_GRID_PT_, 0 );
      //E_pme_.resize( numthreads_);
      //U_E_pme_.resize(numthreads_);
      //for (int thread = 0; thread != numthreads_; thread++) {
      //  E_pme_[thread].assign( MAX_GRID_PT_,0);
      //  U_E_pme_[thread].assign( MAX_GRID_PT_,0);
      //}
    }
#   ifdef _OPENMP
    if (doEij_) {
      // Since allocating a separate matrix for every thread will consume a lot
      // of memory and since the Eij matrices tend to be sparse since solute is
      // often present, each thread will record any interaction energies they
      // calculate separately and add to the Eij matrix afterwards to avoid
      // memory clashes. Probably not ideal if the bulk of the grid is water however.
      EIJ_V1_.resize( numthreads_ );
      EIJ_V2_.resize( numthreads_ );
      EIJ_EN_.resize( numthreads_ );
    }
#   endif

    #ifdef CUDA
    if (this->skipE_ && this->doOrder_) {
      mprintf("When the keyword \"skipE\" is supplied, \"doorder\" cannot be"
              " chosen, as both calculations are done on the GPU at the same"
              " time.\nIgnoring \"doorder!\"\n");
    }
    #endif
  }

  //Box gbox;
  //gbox.SetBetaLengths( 90.0, (double)nx * gridspacing_,
  //                           (double)ny * gridspacing_,
  //                           (double)nz * gridspacing_ );
  //grid_.Setup_O_Box( nx, ny, nz, gO_->GridOrigin(), gbox );
  //grid_.Setup_O_D( nx, ny, nz, gO_->GridOrigin(), v_spacing );

  mprintf("    GIST:\n");
  mprintf("\tOutput prefix= '%s', grid output extension= '%s'\n", prefix_.c_str(), ext.c_str());
  mprintf("\tOutput float format string= '%s', output integer format string= '%s'\n", fltFmt_.fmt(), intFmt_.fmt());
  mprintf("\tGIST info written to '%s'\n", infofile_->Filename().full());
  mprintf("\tName for data sets: %s\n", dsname.c_str());
  if (doOrder_)
    mprintf("\tDoing order calculation.\n");
  else
    mprintf("\tSkipping order calculation.\n");
  if (skipE_)
    mprintf("\tSkipping energy calculation.\n");
  else {
    mprintf("\tPerforming energy calculation.\n");
    if (numthreads_ > 1)
      mprintf("\tParallelizing energy calculation with %i threads.\n", numthreads_);
    if (usePme_) {
      mprintf("\tUsing PME.\n");
      pmeOpts_.PrintOptions();
    }
  }
  mprintf("\tCut off for determining solvent O-O neighbors is %f Ang\n", sqrt(NeighborCut2_));
  if (includeIons_)
    mprintf("\tIons will be included in the solute region.\n");
  else
    mprintf("\tIons will be excluded from the calculation.\n");
  if (doEij_) {
    mprintf("\tComputing and printing water-water Eij matrix, output to '%s'\n",
            eijfile_->Filename().full());
    mprintf("\tWater-water Eij matrix size is %s\n",
            ByteString(ww_Eij_->MemUsageInBytes(), BYTE_DECIMAL).c_str());
  } else
    mprintf("\tSkipping water-water Eij matrix.\n");
  mprintf("\tWater reference density: %6.4f molecules/Ang^3\n", BULK_DENS_);
  mprintf("\tSimulation temperature: %6.4f K\n", temperature_);
  if (imageOpt_.UseImage())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tDistances will not be imaged.\n");
  if (imageOpt_.ForceNonOrtho())
    mprintf("\tWill use non-orthogonal imaging routines for all cell types.\n");
  gO_->GridInfo();
  mprintf("\tNumber of voxels: %u, voxel volume: %f Ang^3\n",
          MAX_GRID_PT_, gO_->Bin().VoxelVolume());
  mprintf("#Please cite these papers if you use GIST results in a publication:\n"
          "#    Steven Ramsey, Crystal Nguyen, Romelia Salomon-Ferrer, Ross C. Walker, Michael K. Gilson, and Tom Kurtzman. J. Comp. Chem. 37 (21) 2016\n"
          "#    Crystal Nguyen, Michael K. Gilson, and Tom Young, arXiv:1108.4876v1 (2011)\n"
          "#    Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson,\n"
          "#      J. Chem. Phys. 137, 044101 (2012)\n"
          "#    Lazaridis, J. Phys. Chem. B 102, 3531–3541 (1998)\n"
#ifdef LIBPME
          "#When using the PME-enhanced version of GIST, please cite:\n"
          "#    Lieyang Chen, Anthony Cruz, Daniel R. Roe, Andy C. Simmonett, Lauren Wickstrom, Nanjie Deng, Tom Kurtzman. JCTC (2021) DOI: 10.1021/acs.jctc.0c01185\n"
#endif
#ifdef CUDA
          "#When using the GPU parallelized version of GIST, please cite:\n"
          "#    Johannes Kraml, Anna S. Kamenik, Franz Waibl, Michael Schauperl, Klaus R. Liedl, JCTC (2019)\n"
#endif
          );
# ifdef GIST_USE_NONORTHO_DIST2
  mprintf("DEBUG: Using regular non-orthogonal distance routine.\n");
# endif
  gist_init_.Stop();
  return Action::OK;
}

/// \return True if given floating point values are not equal within a tolerance
static inline bool NotEqual(double v1, double v2) { return ( fabs(v1 - v2) > Constants::SMALL ); }

/** Set up GIST action. */
Action::RetType Action_GIST::Setup(ActionSetup& setup) {
  gist_setup_.Start();
  CurrentParm_ = setup.TopAddress();
  // We need box info
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprinterr("Error: Must have explicit solvent with periodic boundaries!");
    return Action::ERR;
  }
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );
  #ifdef CUDA
  this->numberAtoms_ = setup.Top().Natom();
  this->solvent_ = new bool[this->numberAtoms_];
  #endif

  // Initialize PME
  if (usePme_) {
#   ifdef LIBPME
    if (gistPme_.Init( setup.CoordInfo().TrajBox(), pmeOpts_, debug_ )) {
      mprinterr("Error: GIST PME init failed.\n");
      return Action::ERR;
    }
    // By default all atoms are selected for GIST PME to match up with atom_voxel_ array.
    if (gistPme_.Setup_PME_GIST( setup.Top(), numthreads_, NeighborCut2_ )) {
      mprinterr("Error: GIST PME setup/array allocation failed.\n");
      return Action::ERR;
    }
#   else
    mprinterr("Error: GIST PME requires compilation with LIBPME.\n");
    return Action::ERR;
#   endif
  }

  // Get molecule number for each solvent molecule
  //mol_nums_.clear();
  O_idxs_.clear();
  A_idxs_.clear();
  atom_voxel_.clear();
  atomIsSolute_.clear();
  atomIsSolventO_.clear();
  U_idxs_.clear();
  // NOTE: these are just guesses
  O_idxs_.reserve( setup.Top().Nsolvent() );
  A_idxs_.reserve( setup.Top().Natom() );
  // atom_voxel_ and atomIsSolute will be indexed by atom #
  atom_voxel_.assign( setup.Top().Natom(), OFF_GRID_ );
  atomIsSolute_.assign(setup.Top().Natom(), false);
  atomIsSolventO_.assign(setup.Top().Natom(), false);
  U_idxs_.reserve(setup.Top().Natom()-setup.Top().Nsolvent()*nMolAtoms_);
  unsigned int midx = 0;
  unsigned int NsolventAtoms = 0;
  unsigned int NsoluteAtoms = 0;
  bool isFirstSolvent = true;
  for (Topology::mol_iterator mol = setup.Top().MolStart();
                              mol != setup.Top().MolEnd(); ++mol, ++midx)
  {
    if (mol->IsSolvent()) {
      // NOTE: We assume the oxygen is the first atom!
      int o_idx = mol->MolUnit().Front();
      #ifdef CUDA
      this->headAtomType_ = setup.Top()[o_idx].TypeIndex();
      #endif
      // Check that molecule has correct # of atoms
      unsigned int molNumAtoms = (unsigned int)mol->NumAtoms();
      if (nMolAtoms_ == 0) {
        nMolAtoms_ = molNumAtoms;
        mprintf("\tEach solvent molecule has %u atoms\n", nMolAtoms_);
      } else if (molNumAtoms != nMolAtoms_) {
        mprinterr("Error: All solvent molecules must have same # atoms.\n"
                  "Error: Molecule '%s' has %u atoms, expected %u.\n",
                  setup.Top().TruncResNameNum( setup.Top()[o_idx].ResNum() ).c_str(),
                  molNumAtoms, nMolAtoms_);
        return Action::ERR;
      }
      //mol_nums_.push_back( midx ); // TODO needed?
      // Check that first atom is actually Oxygen
      if (setup.Top()[o_idx].Element() != Atom::OXYGEN) {
        mprinterr("Error: Molecule '%s' is not water or does not have oxygen atom.\n",
                  setup.Top().TruncResNameNum( setup.Top()[o_idx].ResNum() ).c_str());
        return Action::ERR;
      }
      O_idxs_.push_back( o_idx );
      atomIsSolventO_[o_idx] = true;
      // Check that the next two atoms are Hydrogens
      if (setup.Top()[o_idx+1].Element() != Atom::HYDROGEN ||
          setup.Top()[o_idx+2].Element() != Atom::HYDROGEN)
      {
        mprinterr("Error: Molecule '%s' does not have hydrogen atoms.\n",
                  setup.Top().TruncResNameNum( setup.Top()[o_idx].ResNum() ).c_str());
        return Action::ERR;
      }
      // Save all atom indices for energy calc, including extra points
      for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
        A_idxs_.push_back( o_idx + IDX );
        atomIsSolute_[A_idxs_.back()] = false; // The identity of the atom is water
        atom_voxel_[A_idxs_.back()] = OFF_GRID_;
        #ifdef CUDA
        this->molecule_.push_back( setup.Top()[o_idx + IDX ].MolNum() );
        this->charges_.push_back( setup.Top()[o_idx + IDX ].Charge() );
        this->atomTypes_.push_back( setup.Top()[o_idx + IDX ].TypeIndex() );
        this->solvent_[ o_idx + IDX ] = true;
        #endif
      }
      NsolventAtoms += nMolAtoms_;
      // If first solvent molecule, save charges. If not, check that charges match.
      if (isFirstSolvent) {
        double q_sum = 0.0;
        Q_.reserve( nMolAtoms_ );
        for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
          Q_.push_back( setup.Top()[o_idx+IDX].Charge() );
          q_sum += Q_.back();
          //mprintf("DEBUG: Q= %20.10E  q_sum= %20.10E\n", setup.Top()[o_idx+IDX].Charge(), q_sum);
        }
        // Sanity checks.
        // NOTE: We know indices 1 and 2 are hydrogens (with 0 being oxygen); this is checked above.
        if (NotEqual(Q_[1], Q_[2]))
          mprintf("Warning: Charges on water hydrogens do not match (%g, %g).\n", Q_[1], Q_[2]);
        if (fabs( q_sum ) > 0.0)
          mprintf("Warning: Charges on water do not sum to 0 (%g)\n", q_sum);
        //mprintf("DEBUG: Water charges: O=%g  H1=%g  H2=%g\n", q_O_, q_H1_, q_H2_);
      } else {
        for (unsigned int IDX = 0; IDX < nMolAtoms_; IDX++) {
          double q_atom = setup.Top()[o_idx+IDX].Charge();
          if (NotEqual(Q_[IDX], q_atom)) {
            mprintf("Warning: Charge on water '%s' (%g) does not match first water (%g).\n",
                  setup.Top().TruncResAtomName( o_idx+IDX ).c_str(), q_atom, Q_[IDX]);
          }
        }
      }
      isFirstSolvent = false;
    } else {
      // This is a non-solvent molecule. Save atom indices. May want to exclude
      // if only 1 atom (probably ion).
      if (mol->NumAtoms() > 1 || includeIons_) {
        for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                  seg != mol->MolUnit().segEnd(); ++seg)
        {
          for (int u_idx = seg->Begin(); u_idx != seg->End(); ++u_idx) {
            A_idxs_.push_back( u_idx );
            atomIsSolute_[A_idxs_.back()] = true; // the identity of the atom is solute
            NsoluteAtoms++;
            U_idxs_.push_back( u_idx ); // store the solute atom index for locating voxel index
            #ifdef CUDA
            this->molecule_.push_back( setup.Top()[ u_idx ].MolNum() );
            this->charges_.push_back( setup.Top()[ u_idx ].Charge() );
            this->atomTypes_.push_back( setup.Top()[ u_idx ].TypeIndex() );
            this->solvent_[ u_idx ] = false;
            #endif
          }
        }
      }
    }
  }
  NSOLVENT_ = O_idxs_.size();
  mprintf("\t%zu solvent molecules, %u solvent atoms, %u solute atoms (%zu total).\n",
          O_idxs_.size(), NsolventAtoms, NsoluteAtoms, A_idxs_.size());
  if (doOrder_ && NSOLVENT_ < 5) {
    mprintf("Warning: Less than 5 solvent molecules. Cannot perform order calculation.\n");
    doOrder_ = false;
  }
  // Allocate space for saving indices of water atoms that are on the grid
  // Estimate how many solvent molecules can possibly fit onto the grid.
  // Add some extra voxels as a buffer.
  double max_voxels = (double)MAX_GRID_PT_ + (1.10 * (double)MAX_GRID_PT_);
  double totalVolume = max_voxels * gO_->Bin().VoxelVolume();
  double max_mols = totalVolume * BULK_DENS_;
  //mprintf("\tEstimating grid can fit a max of %.0f solvent molecules (w/ 10%% buffer).\n",
  //        max_mols);
  OnGrid_idxs_.reserve( (size_t)max_mols * (size_t)nMolAtoms_ );
  N_ON_GRID_ = 0;

  if (!skipE_) {
    if (imageOpt_.ImagingEnabled())
      mprintf("\tImaging enabled for energy distance calculations.\n");
    else
      mprintf("\tNo imaging will be performed for energy distance calculations.\n");
  }

#ifdef CUDA
  NonbondParmType nb = setup.Top().Nonbond();
  this->NBIndex_ = nb.NBindex();
  this->numberAtomTypes_ = nb.Ntypes();
  for (unsigned int i = 0; i < nb.NBarray().size(); ++i) {
    this->lJParamsA_.push_back( (float) nb.NBarray().at(i).A() );
    this->lJParamsB_.push_back( (float) nb.NBarray().at(i).B() );
  }

  try {
    allocateCuda(((void**)&this->NBindex_c_), this->NBIndex_.size() * sizeof(int));
    allocateCuda((void**)&this->max_c_, 3 * sizeof(float));
    allocateCuda((void**)&this->min_c_, 3 * sizeof(float));
    allocateCuda((void**)&this->result_w_c_, this->numberAtoms_ * sizeof(float));
    allocateCuda((void**)&this->result_s_c_, this->numberAtoms_ * sizeof(float));
    allocateCuda((void**)&this->result_O_c_, this->numberAtoms_ * 4 * sizeof(int));
    allocateCuda((void**)&this->result_N_c_, this->numberAtoms_ * sizeof(int));
  } catch (CudaException &e) {
    mprinterr("Error: Could not allocate memory on GPU!\n");
    this->freeGPUMemory();
    return Action::ERR;
  }
  try {
    this->copyToGPU();
  } catch (CudaException &e) {
    return Action::ERR;
  }
#endif

  gist_setup_.Stop();
  return Action::OK;
}

const Vec3 Action_GIST::x_lab_ = Vec3(1.0, 0.0, 0.0);
const Vec3 Action_GIST::y_lab_ = Vec3(0.0, 1.0, 0.0);
const Vec3 Action_GIST::z_lab_ = Vec3(0.0, 0.0, 1.0);
const double Action_GIST::QFAC_ = Constants::ELECTOAMBER * Constants::ELECTOAMBER;
const int Action_GIST::OFF_GRID_ = -1;

/* Calculate the charge-charge, vdw interaction using pme, frame by frame
 * 
 */
void Action_GIST::NonbondEnergy_pme(Frame const& frameIn)
{
# ifdef LIBPME
  // Two energy terms for the whole system
  //double ene_pme_all = 0.0;
  //double ene_vdw_all = 0.0;
  // pointer to the E_pme_, where has the voxel-wise pme energy for water
  double* E_pme_grid = &E_pme_[0];
  // pointer to U_E_pme_, where has the voxel-wise pme energy for solute
  double* U_E_pme_grid = &U_E_pme_[0]; 

  gistPme_.CalcNonbondEnergy_GIST(frameIn, atom_voxel_, atomIsSolute_, atomIsSolventO_,
                                  E_UV_VDW_, E_UV_Elec_, E_VV_VDW_, E_VV_Elec_,
                                  neighbor_);

//  system_potential_energy_ += ene_pme_all + ene_vdw_all;

  // Water energy on the GIST grid
  double pme_sum = 0.0;

  for (unsigned int gidx=0; gidx < N_ON_GRID_; gidx++ ) 
  {
    int a = OnGrid_idxs_[gidx]; // index of the atom of on-grid solvent;
    int a_voxel = atom_voxel_[a]; // index of the voxel
    double nonbond_energy = gistPme_.E_of_atom(a);
    pme_sum += nonbond_energy;
    E_pme_grid[a_voxel] += nonbond_energy;  
  }

  // Solute energy on the GIST grid
  double solute_on_grid_sum = 0.0; // To sum up the potential energy on solute atoms that on the grid

  for (unsigned int uidx=0; uidx < U_onGrid_idxs_.size(); uidx++ )
  {
    int u = U_onGrid_idxs_[uidx]; // index of the solute atom on the grid
    int u_voxel = atom_voxel_[u];
    double u_nonbond_energy = gistPme_.E_of_atom(u);
    solute_on_grid_sum += u_nonbond_energy; 
    U_E_pme_grid[u_voxel] += u_nonbond_energy;
  }

/*
  // Total solute energy
  double solute_sum = 0.0;

  for (unsigned int uidx=0; uidx < U_idxs_.size(); uidx++)
  {
    int u = U_idxs_[uidx];
    double u_nonbond_energy = gistPme_.E_of_atom(u);
    solute_sum += u_nonbond_energy;
    solute_potential_energy_ += u_nonbond_energy; // used to calculated the ensemble energy for all solute, will print out in terminal
  }
*/
  //mprintf("The total potential energy on water atoms: %f \n", pme_sum);
# else /*LIBPME */
  mprinterr("Error: Compiled without LIBPME\n");
  return;
# endif /*LIBPME */
}

/** Non-bonded energy calc. */
void Action_GIST::Ecalc(double rij2, double q1, double q2, NonbondType const& LJ,
                        double& Evdw, double& Eelec)
{
  double rij = sqrt(rij2);
  // VDW
  double r2    = 1.0 / rij2;
  double r6    = r2 * r2 * r2;
  double r12   = r6 * r6;
  double f12   = LJ.A() * r12;  // A/r^12
  double f6    = LJ.B() * r6;   // B/r^6
         Evdw  = f12 - f6;      // (A/r^12)-(B/r^6)
  // Coulomb
  double qiqj  = QFAC_ * q1 * q2;
         Eelec = qiqj / rij;
}

/** Calculate the energy between all solute/solvent atoms and solvent atoms
  * on the grid. This is done after the intial GIST calculations
  * so that all waters have voxels assigned in atom_voxel_.
  * NOTE: This routine modifies the coordinates in OnGrid_XYZ_ when the cell
  *       has nonorthogonal shape in order to properly satsify the minimum
  *       image convention, so any calculations that rely on the on grid
  *       coordinates (like Order()) must be done *BEFORE* this routine.
  */
void Action_GIST::NonbondEnergy(Frame const& frameIn, Topology const& topIn)
{
  // Set up imaging info.
  if (imageOpt_.ImagingType() == ImageOption::NONORTHO) {
    // Wrap on-grid water coords back to primary cell TODO openmp
    double* ongrid_xyz = &OnGrid_XYZ_[0];
    int maxXYZ = (int)OnGrid_XYZ_.size();
    int idx;
#   ifdef _OPENMP
#   pragma omp parallel private(idx)
    {
#   pragma omp for
#   endif
    for (idx = 0; idx < maxXYZ; idx += 3)
    {
      double* XYZ = ongrid_xyz + idx;
      // Convert to frac coords
      frameIn.BoxCrd().FracCell().TimesVec( XYZ, XYZ );
      // Wrap to primary cell
      XYZ[0] = XYZ[0] - floor(XYZ[0]);
      XYZ[1] = XYZ[1] - floor(XYZ[1]);
      XYZ[2] = XYZ[2] - floor(XYZ[2]);
      // Convert back to Cartesian
      frameIn.BoxCrd().UnitCell().TransposeMult( XYZ, XYZ );
    }
#   ifdef _OPENMP
    }
#   endif
  }

//  mprintf("DEBUG: NSolventAtoms= %zu  NwatAtomsOnGrid= %u\n", O_idxs_.size()*nMolAtoms_, N_ON_GRID_);

  double* E_UV_VDW  = &(E_UV_VDW_[0][0]);
  double* E_UV_Elec = &(E_UV_Elec_[0][0]);
  double* E_VV_VDW  = &(E_VV_VDW_[0][0]);
  double* E_VV_Elec = &(E_VV_Elec_[0][0]);
  float* Neighbor = &(neighbor_[0][0]);
  double Evdw, Eelec;
  int aidx;
  int maxAidx = (int)A_idxs_.size();
  // Loop over all solute + solvent atoms
# ifdef _OPENMP
  int mythread;
  Iarray* eij_v1 = 0;
  Iarray* eij_v2 = 0;
  Farray* eij_en = 0;
# pragma omp parallel private(aidx, mythread, E_UV_VDW, E_UV_Elec, E_VV_VDW, E_VV_Elec, Neighbor, Evdw, Eelec, eij_v1, eij_v2, eij_en)
  {
  mythread = omp_get_thread_num();
  E_UV_VDW = &(E_UV_VDW_[mythread][0]);
  E_UV_Elec = &(E_UV_Elec_[mythread][0]);
  E_VV_VDW = &(E_VV_VDW_[mythread][0]);
  E_VV_Elec = &(E_VV_Elec_[mythread][0]);
  Neighbor = (&neighbor_[mythread][0]);
  if (doEij_) {
    eij_v1 = &(EIJ_V1_[mythread]);
    eij_v2 = &(EIJ_V2_[mythread]);
    eij_en = &(EIJ_EN_[mythread]);
    eij_v1->clear();
    eij_v2->clear();
    eij_en->clear();
  }
# pragma omp for
# endif
  for (aidx = 0; aidx < maxAidx; aidx++)
  {
    int a1 = A_idxs_[aidx];            // Index of atom1
    int a1_voxel = atom_voxel_[a1];    // Voxel of atom1
    int a1_mol = topIn[ a1 ].MolNum(); // Molecule # of atom 1
    Vec3 A1_XYZ( frameIn.XYZ( a1 ) );  // Coord of atom1
    double qA1 = topIn[ a1 ].Charge(); // Charge of atom1
    bool a1IsO = atomIsSolventO_[a1];
    std::vector<Vec3> vImages;
    if (imageOpt_.ImagingType() == ImageOption::NONORTHO) {
      // Convert to frac coords
      Vec3 vFrac = frameIn.BoxCrd().FracCell() * A1_XYZ;
      // Wrap to primary unit cell
      vFrac[0] = vFrac[0] - floor(vFrac[0]);
      vFrac[1] = vFrac[1] - floor(vFrac[1]);
      vFrac[2] = vFrac[2] - floor(vFrac[2]);
      // Calculate all images of this atom
      vImages.reserve(27);
      for (int ix = -1; ix != 2; ix++)
        for (int iy = -1; iy != 2; iy++)
          for (int iz = -1; iz != 2; iz++)
            // Convert image back to Cartesian
            vImages.push_back( frameIn.BoxCrd().UnitCell().TransposeMult( vFrac + Vec3(ix,iy,iz) ) );
    }
    // Loop over all solvent atoms on the grid
    for (unsigned int gidx = 0; gidx < N_ON_GRID_; gidx++)
    {
      int a2 = OnGrid_idxs_[gidx];              // Index of on-grid solvent
      int a2_mol = topIn[ a2 ].MolNum();        // Molecule # of on-grid solvent
      if (a1_mol != a2_mol)
      {
        int a2_voxel = atom_voxel_[a2];                  // Voxel of on-grid solvent
        const double* A2_XYZ = (&OnGrid_XYZ_[0])+gidx*3; // Coord of on-grid solvent
        if (atomIsSolute_[a1]) {
          // Solute to on-grid solvent energy
          // Calculate distance
          //gist_nonbond_dist_.Start();
          double rij2;
          if (imageOpt_.ImagingType() == ImageOption::NONORTHO) {
#           ifdef GIST_USE_NONORTHO_DIST2
            rij2 = DIST2_ImageNonOrtho(A1_XYZ, A2_XYZ, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell());
#           else
            rij2 = maxD_;
            for (std::vector<Vec3>::const_iterator vCart = vImages.begin();
                                                   vCart != vImages.end(); ++vCart)
            {
              double x = (*vCart)[0] - A2_XYZ[0];
              double y = (*vCart)[1] - A2_XYZ[1];
              double z = (*vCart)[2] - A2_XYZ[2];
              rij2 = std::min(rij2, x*x + y*y + z*z);
            }
#           endif
          } else if (imageOpt_.ImagingType() == ImageOption::ORTHO)
            rij2 = DIST2_ImageOrtho( A1_XYZ, A2_XYZ, frameIn.BoxCrd() );
          else
            rij2 = DIST2_NoImage( A1_XYZ, A2_XYZ );
          //gist_nonbond_dist_.Stop();
          //gist_nonbond_UV_.Start();
          // Calculate energy
          Ecalc( rij2, qA1, topIn[ a2 ].Charge(), topIn.GetLJparam(a1, a2), Evdw, Eelec );
          E_UV_VDW[a2_voxel]  += Evdw;
          E_UV_Elec[a2_voxel] += Eelec;
          //gist_nonbond_UV_.Stop();
        } else {
          // Off-grid/on-grid solvent to on-grid solvent energy
          // Only do the energy calculation if not previously done or atom1 not on grid
          if (a2 != a1 && (a2 > a1 || a1_voxel == OFF_GRID_))
          {
            // Calculate distance
            //gist_nonbond_dist_.Start();
            double rij2;
            if (imageOpt_.ImagingType() == ImageOption::NONORTHO) {
#            ifdef GIST_USE_NONORTHO_DIST2
             rij2 = DIST2_ImageNonOrtho(A1_XYZ, A2_XYZ, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell());
#            else
             rij2 = maxD_;
              for (std::vector<Vec3>::const_iterator vCart = vImages.begin();
                                                     vCart != vImages.end(); ++vCart)
              {
                double x = (*vCart)[0] - A2_XYZ[0];
                double y = (*vCart)[1] - A2_XYZ[1];
                double z = (*vCart)[2] - A2_XYZ[2];
                rij2 = std::min(rij2, x*x + y*y + z*z);
              }
#             endif
            } else if (imageOpt_.ImagingType() == ImageOption::ORTHO)
              rij2 = DIST2_ImageOrtho( A1_XYZ, A2_XYZ, frameIn.BoxCrd() );
            else
              rij2 = DIST2_NoImage( A1_XYZ, A2_XYZ );
            //gist_nonbond_dist_.Stop();
            //gist_nonbond_VV_.Start();
            // Calculate energy
            Ecalc( rij2, qA1, topIn[ a2 ].Charge(), topIn.GetLJparam(a1, a2), Evdw, Eelec );
            //mprintf("DEBUG1: v1= %i v2= %i EVV %i %i Vdw= %f Elec= %f\n", a2_voxel, a1_voxel, a2, a1, Evdw, Eelec);
            E_VV_VDW[a2_voxel] += Evdw;
            E_VV_Elec[a2_voxel] += Eelec;
            // Store water neighbor using only O-O distance
            bool is_O_O = (a1IsO && atomIsSolventO_[a2]);
            if (is_O_O && rij2 < NeighborCut2_)
              Neighbor[a2_voxel] += 1.0;
            // If water atom1 was also on the grid update its energy as well.
            if ( a1_voxel != OFF_GRID_ ) {
              E_VV_VDW[a1_voxel] += Evdw;
              E_VV_Elec[a1_voxel] += Eelec;
              if (is_O_O && rij2 < NeighborCut2_)
                Neighbor[a1_voxel] += 1.0;
              if (doEij_) {
                if (a1_voxel != a2_voxel) {
#                 ifdef _OPENMP
                  eij_v1->push_back( a1_voxel );
                  eij_v2->push_back( a2_voxel );
                  eij_en->push_back( Evdw + Eelec );
#                 else
                  ww_Eij_->UpdateElement(a1_voxel, a2_voxel, Evdw + Eelec);
#                 endif
                }
              }
            }
            //gist_nonbond_VV_.Stop();
          }
        }
      } // END a1 and a2 not in same molecule
    } // End loop over all solvent atoms on grid
  } // End loop over all solvent + solute atoms
# ifdef _OPENMP
  } // END pragma omp parallel
  if (doEij_) {
    // Add any Eijs to matrix
    for (unsigned int thread = 0; thread != EIJ_V1_.size(); thread++)
      for (unsigned int idx = 0; idx != EIJ_V1_[thread].size(); idx++)
        ww_Eij_->UpdateElement(EIJ_V1_[thread][idx], EIJ_V2_[thread][idx], EIJ_EN_[thread][idx]);
  }
# endif
}

/** GIST order calculation. */
void Action_GIST::Order(Frame const& frameIn) {
  // Loop over all solvent molecules that are on the grid
  for (unsigned int gidx = 0; gidx < N_ON_GRID_; gidx += nMolAtoms_)
  {
    int oidx1 = OnGrid_idxs_[gidx];
    int voxel1 = atom_voxel_[oidx1];
    Vec3 XYZ1( (&OnGrid_XYZ_[0])+gidx*3 );
    // Find coordinates for 4 closest neighbors to this water (on or off grid).
    // TODO set up overall grid in DoAction.
    // TODO initialize WAT?
    Vec3 WAT[4];
    double d1 = maxD_;
    double d2 = maxD_;
    double d3 = maxD_;
    double d4 = maxD_;
    for (unsigned int sidx2 = 0; sidx2 < NSOLVENT_; sidx2++)
    {
      int oidx2 = O_idxs_[sidx2];
      if (oidx2 != oidx1)
      {
        const double* XYZ2 = frameIn.XYZ( oidx2 );
        double dist2 = DIST2_NoImage( XYZ1.Dptr(), XYZ2 );
        if        (dist2 < d1) {
          d4 = d3; d3 = d2; d2 = d1; d1 = dist2;
          WAT[3] = WAT[2]; WAT[2] = WAT[1]; WAT[1] = WAT[0]; WAT[0] = XYZ2;
        } else if (dist2 < d2) {
          d4 = d3; d3 = d2; d2 = dist2;
          WAT[3] = WAT[2]; WAT[2] = WAT[1]; WAT[1] = XYZ2;
        } else if (dist2 < d3) {
          d4 = d3; d3 = dist2;
          WAT[3] = WAT[2]; WAT[2] = XYZ2;
        } else if (dist2 < d4) {
          d4 = dist2;
          WAT[3] = XYZ2;
        }
      }
    }
    // Compute the tetrahedral order parameter
    double sum = 0.0;
    for (int mol1 = 0; mol1 < 3; mol1++) {
      for (int mol2 = mol1 + 1; mol2 < 4; mol2++) {
        Vec3 v1 = WAT[mol1] - XYZ1;
        Vec3 v2 = WAT[mol2] - XYZ1;
        double r1 = v1.Magnitude2();
        double r2 = v2.Magnitude2();
        double cos = (v1* v2) / sqrt(r1 * r2);
        sum += (cos + 1.0/3)*(cos + 1.0/3);
      }
    }
    order_norm_->UpdateVoxel(voxel1, (1.0 - (3.0/8)*sum));
    //mprintf("DBG: gidx= %u  oidx1=%i  voxel1= %i  XYZ1={%g, %g, %g}  sum= %g\n", gidx, oidx1, voxel1, XYZ1[0], XYZ1[1], XYZ1[2], sum);
  } // END loop over all solvent molecules
}

/** GIST action */
Action::RetType Action_GIST::DoAction(int frameNum, ActionFrame& frm) {
  gist_action_.Start();
  NFRAME_++;
  // TODO only !skipE?
  N_ON_GRID_ = 0;
  OnGrid_idxs_.clear();
  OnGrid_XYZ_.clear();

  // Determine imaging type
# ifdef DEBUG_GIST
  //mprintf("DEBUG: Is_X_Aligned_Ortho() = %i  Is_X_Aligned() = %i\n", (int)frm.Frm().BoxCrd().Is_X_Aligned_Ortho(), (int)frm.Frm().BoxCrd().Is_X_Aligned());
  frm.Frm().BoxCrd().UnitCell().Print("Ucell");
  frm.Frm().BoxCrd().FracCell().Print("Frac");
# endif
  if (imageOpt_.ImagingEnabled())
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );
# ifdef DEBUG_GIST
  switch (imageOpt_.ImagingType()) {
    case ImageOption::NO_IMAGE : mprintf("DEBUG: No Image.\n"); break;
    case ImageOption::ORTHO    : mprintf("DEBUG: Orthogonal image.\n"); break;
    case ImageOption::NONORTHO : mprintf("DEBUG: Nonorthogonal image.\n"); break;
  }
# endif
  // CUDA necessary information

  size_t bin_i, bin_j, bin_k;
  Vec3 const& Origin = gO_->Bin().GridOrigin();
  // Loop over each solvent molecule
  for (unsigned int sidx = 0; sidx < NSOLVENT_; sidx++)
  {
    gist_grid_.Start();
    int oidx = O_idxs_[sidx];
    for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++)
      atom_voxel_[oidx+IDX] = OFF_GRID_;
    const double* O_XYZ  = frm.Frm().XYZ( oidx );
    // Get vector of water oxygen to grid origin.
    Vec3 W_G( O_XYZ[0] - Origin[0],
              O_XYZ[1] - Origin[1],
              O_XYZ[2] - Origin[2] );
    gist_grid_.Stop();
    // Check if water oxygen is no more then 1.5 Ang from grid
    // NOTE: using <= to be consistent with original code
    if ( W_G[0] <= G_max_[0] && W_G[0] >= -1.5 &&
         W_G[1] <= G_max_[1] && W_G[1] >= -1.5 &&
         W_G[2] <= G_max_[2] && W_G[2] >= -1.5 )
    {
      const double* H1_XYZ = frm.Frm().XYZ( oidx + 1 );
      const double* H2_XYZ = frm.Frm().XYZ( oidx + 2 );
      // Try to bin the oxygen
      if ( gO_->Bin().Calc( O_XYZ[0], O_XYZ[1], O_XYZ[2], bin_i, bin_j, bin_k ) )
      {
        // Oxygen is inside the grid. Record the voxel.
        // NOTE hydrogens/EP always assigned to same voxel for energy purposes.
        int voxel = (int)gO_->CalcIndex(bin_i, bin_j, bin_k);
        const double* wXYZ = O_XYZ;
        for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
          atom_voxel_[oidx+IDX] = voxel;
          //OnGrid_idxs_[N_ON_GRID_+IDX] = oidx + IDX;
          OnGrid_idxs_.push_back( oidx+IDX );
          OnGrid_XYZ_.push_back( wXYZ[0] );
          OnGrid_XYZ_.push_back( wXYZ[1] );
          OnGrid_XYZ_.push_back( wXYZ[2] );
          wXYZ+=3;
        }
        N_ON_GRID_ += nMolAtoms_;
        //mprintf("DEBUG1: Water atom %i voxel %i\n", oidx, voxel);
        N_waters_[voxel]++;
        max_nwat_ = std::max( N_waters_[voxel], max_nwat_ );
        // ----- EULER ---------------------------
        gist_euler_.Start();
        // Record XYZ coords of water atoms (nonEP) in voxel TODO need EP?
        voxel_xyz_[voxel].push_back( O_XYZ[0] );
        voxel_xyz_[voxel].push_back( O_XYZ[1] );
        voxel_xyz_[voxel].push_back( O_XYZ[2] );
        // Get O-HX vectors
        Vec3 H1_wat( H1_XYZ[0]-O_XYZ[0], H1_XYZ[1]-O_XYZ[1], H1_XYZ[2]-O_XYZ[2] );
        Vec3 H2_wat( H2_XYZ[0]-O_XYZ[0], H2_XYZ[1]-O_XYZ[1], H2_XYZ[2]-O_XYZ[2] );
        H1_wat.Normalize();
        H2_wat.Normalize();

        Vec3 ar1 = H1_wat.Cross( x_lab_ ); // ar1 = V cross U 
        Vec3 sar = ar1;                    // sar = V cross U
        ar1.Normalize();
        //mprintf("------------------------------------------\n");
        //H1_wat.Print("DEBUG: H1_wat");
        //x_lab_.Print("DEBUG: x_lab_");
        //ar1.Print("DEBUG: ar1");
        //sar.Print("DEBUG: sar");
        double dp1 = x_lab_ * H1_wat; // V dot U
        double theta = acos(dp1);
        double sign = sar * H1_wat;
        //mprintf("DEBUG0: dp1= %f  theta= %f  sign= %f\n", dp1, theta, sign);
        // NOTE: Use SMALL instead of 0 to avoid issues with denormalization
        if (sign > Constants::SMALL)
          theta /= 2.0;
        else
          theta /= -2.0;
        double w1 = cos(theta);
        double sin_theta = sin(theta);
        //mprintf("DEBUG0: theta= %f  w1= %f  sin_theta= %f\n", theta, w1, sin_theta);
        double x1 = ar1[0] * sin_theta;
        double y1 = ar1[1] * sin_theta;
        double z1 = ar1[2] * sin_theta;
        double w2 = w1;
        double x2 = x1;
        double y2 = y1;
        double z2 = z1;

        Vec3 H_temp;
        H_temp[0] = ((w2*w2+x2*x2)-(y2*y2+z2*z2))*H1_wat[0];
        H_temp[0] = (2*(x2*y2 + w2*z2)*H1_wat[1]) + H_temp[0];
        H_temp[0] = (2*(x2*z2-w2*y2)*H1_wat[2]) + H_temp[0];

        H_temp[1] = 2*(x2*y2 - w2*z2)* H1_wat[0];
        H_temp[1] = ((w2*w2-x2*x2+y2*y2-z2*z2)*H1_wat[1]) + H_temp[1];
        H_temp[1] = (2*(y2*z2+w2*x2)*H1_wat[2]) +H_temp[1];

        H_temp[2] = 2*(x2*z2+w2*y2) *H1_wat[0];
        H_temp[2] = (2*(y2*z2-w2*x2)*H1_wat[1]) + H_temp[2];
        H_temp[2] = ((w2*w2-x2*x2-y2*y2+z2*z2)*H1_wat[2]) + H_temp[2];

        H1_wat = H_temp;

        Vec3 H_temp2;
        H_temp2[0] = ((w2*w2+x2*x2)-(y2*y2+z2*z2))*H2_wat[0];
        H_temp2[0] = (2*(x2*y2 + w2*z2)*H2_wat[1]) + H_temp2[0];
        H_temp2[0] = (2*(x2*z2-w2*y2)*H2_wat[2]) +H_temp2[0];

        H_temp2[1] = 2*(x2*y2 - w2*z2) *H2_wat[0];
        H_temp2[1] = ((w2*w2-x2*x2+y2*y2-z2*z2)*H2_wat[1]) +H_temp2[1];
        H_temp2[1] = (2*(y2*z2+w2*x2)*H2_wat[2]) +H_temp2[1];

        H_temp2[2] = 2*(x2*z2+w2*y2)*H2_wat[0];
        H_temp2[2] = (2*(y2*z2-w2*x2)*H2_wat[1]) +H_temp2[2];
        H_temp2[2] = ((w2*w2-x2*x2-y2*y2+z2*z2)*H2_wat[2]) + H_temp2[2];

        H2_wat = H_temp2;

        Vec3 ar2 = H_temp.Cross(H_temp2);
        ar2.Normalize();
        double dp2 = ar2 * z_lab_;
        theta = acos(dp2);

        sar = ar2.Cross( z_lab_ );
        sign = sar * H_temp;

        if (sign < 0)
          theta /= 2.0;
        else
          theta /= -2.0;

        double w3 = cos(theta);
        sin_theta = sin(theta);
        double x3 = x_lab_[0] * sin_theta;
        double y3 = x_lab_[1] * sin_theta;
        double z3 = x_lab_[2] * sin_theta;

        double w4 = w1*w3 - x1*x3 - y1*y3 - z1*z3;
        double x4 = w1*x3 + x1*w3 + y1*z3 - z1*y3;
        double y4 = w1*y3 - x1*z3 + y1*w3 + z1*x3;
        double z4 = w1*z3 + x1*y3 - y1*x3 + z1*w3;

        voxel_Q_[voxel].push_back( w4 );
        voxel_Q_[voxel].push_back( x4 );
        voxel_Q_[voxel].push_back( y4 );
        voxel_Q_[voxel].push_back( z4 );
        //mprintf("DEBUG1: sidx= %u  voxel= %i  wxyz4= %g %g %g %g\n", sidx, voxel, w4, x4, y4, z4);
        //mprintf("DEBUG2: wxyz3= %g %g %g %g  wxyz2= %g %g %g %g  wxyz1= %g %g %g\n",
        //        w3, x3, y3, z3,
        //        w2, x2, y2, z2,
        //        w1, x1, y1, z1);
        // NOTE: No need for nw_angle_ here, it is same as N_waters_
        gist_euler_.Stop();
        // ----- DIPOLE --------------------------
        gist_dipole_.Start();
        //mprintf("DEBUG1: voxel %i dipole %f %f %f\n", voxel,
        //        O_XYZ[0]*q_O_ + H1_XYZ[0]*q_H1_ + H2_XYZ[0]*q_H2_,
        //        O_XYZ[1]*q_O_ + H1_XYZ[1]*q_H1_ + H2_XYZ[1]*q_H2_,
        //        O_XYZ[2]*q_O_ + H1_XYZ[2]*q_H1_ + H2_XYZ[2]*q_H2_);
        double DPX = 0.0;
        double DPY = 0.0;
        double DPZ = 0.0;
        for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
          const double* XYZ = frm.Frm().XYZ( oidx+IDX );
          DPX += XYZ[0] * Q_[IDX];
          DPY += XYZ[1] * Q_[IDX];
          DPZ += XYZ[2] * Q_[IDX];
        }
        dipolex_->UpdateVoxel(voxel, DPX);
        dipoley_->UpdateVoxel(voxel, DPY);
        dipolez_->UpdateVoxel(voxel, DPZ);
        gist_dipole_.Stop();
        // ---------------------------------------
      }

      // Water is at most 1.5A away from grid, so we need to check for H
      // even if O is outside grid.
      if (gO_->Bin().Calc( H1_XYZ[0], H1_XYZ[1], H1_XYZ[2], bin_i, bin_j, bin_k ) )
        N_hydrogens_[ (int)gO_->CalcIndex(bin_i, bin_j, bin_k) ]++;
      if (gO_->Bin().Calc( H2_XYZ[0], H2_XYZ[1], H2_XYZ[2], bin_i, bin_j, bin_k ) )
        N_hydrogens_[ (int)gO_->CalcIndex(bin_i, bin_j, bin_k) ]++;
    } // END water is within 1.5 Ang of grid
  } // END loop over each solvent molecule

  // Do solute grid assignment for PME
  if (usePme_) {
    U_onGrid_idxs_.clear();
    gist_grid_.Start();
    for (unsigned int s = 0; s != U_idxs_.size(); s++)
    {
      int uidx = U_idxs_[s]; // the solute atom index
      atom_voxel_[uidx] = OFF_GRID_;
      const double* u_XYZ = frm.Frm().XYZ( uidx );
      // get the vector of this solute atom to the grid origin
      Vec3 U_G( u_XYZ[0] - Origin[0],
                u_XYZ[1] - Origin[1],
                u_XYZ[2] - Origin[2]);
      //size_t bin_i, bin_j, bin_k;

      if ( U_G[0] <= G_max_[0] && U_G[0] >= -1.5 &&
           U_G[1] <= G_max_[1] && U_G[1] >= -1.5 &&
           U_G[2] <= G_max_[2] && U_G[2] >- -1.5)
      {
        if ( gO_->Bin().Calc(u_XYZ[0],u_XYZ[1],u_XYZ[2],bin_i,bin_j,bin_k))  // used the gO class function to calcaute voxel index
        {
          int voxel = (int)gO_->CalcIndex(bin_i,bin_j,bin_k);
          atom_voxel_[uidx] = voxel;   // asign the voxel index to the solute atom
          //U_ON_GRID_ +=1;           // add +1 to the number of atom on the GIST Grid
          N_solute_atoms_[voxel] +=1;  // add +1 to the solute atom num in this voxel
          U_onGrid_idxs_.push_back(uidx); // The index of the solute atom on GIST Grid
        }
      }
    }
    gist_grid_.Stop();
  }

# ifndef CUDA
  // Do order calculation if requested.
  // Do not do this for CUDA since CUDA nonbond routine handles the order calc.
  // NOTE: This has to be done before the nonbond energy calc since
  //       the nonbond calc can modify the on-grid coordinates (for minimum
  //       image convention when cell is non-orthogonal).
  gist_order_.Start();
  if (doOrder_) Order(frm.Frm());
  gist_order_.Stop();
# endif
  // Do nonbond energy calc if not skipping energy
  gist_nonbond_.Start();
  if (!skipE_) {
    if (usePme_) {
      // PME
      NonbondEnergy_pme( frm.Frm() );
    } else {
      // Non-PME
#     ifdef CUDA
      NonbondCuda(frm);
#     else
      NonbondEnergy(frm.Frm(), *CurrentParm_);
#     endif
    }
  }
  gist_nonbond_.Stop();

  gist_action_.Stop();
  return Action::OK;
}

/** Translational entropy calc between given water and all waters in voxel 2.
  * \param VX voxel 1 water X
  * \param VY voxel 1 water Y
  * \param VZ voxel 1 water Z
  * \param W4 voxel 1 water W4
  * \param X4 voxel 1 water X4
  * \param Y4 voxel 1 water Y4
  * \param Z4 voxel 1 water Z4
  * \param voxel2 Index of second voxel
  */
void Action_GIST::TransEntropy(float VX, float VY, float VZ,
                               float W4, float X4, float Y4, float Z4,
                               int voxel2, double& NNd, double& NNs) const
{
  int nw_tot = N_waters_[voxel2];
  Farray const& V_XYZ = voxel_xyz_[voxel2];
  Farray const& V_Q   = voxel_Q_[voxel2];
  for (int n1 = 0; n1 != nw_tot; n1++)
  {
    int i1 = n1 * 3; // index into V_XYZ for n1
    double dx = (double)(VX - V_XYZ[i1  ]);
    double dy = (double)(VY - V_XYZ[i1+1]);
    double dz = (double)(VZ - V_XYZ[i1+2]);
    double dd = dx*dx+dy*dy+dz*dz;
    if (dd < NNd && dd > 0) { NNd = dd; }

    int q1 = n1 * 4; // index into V_Q for n1
    double rR = 2.0 * acos( fabs(W4 * V_Q[q1  ] +
                            X4 * V_Q[q1+1] +
                            Y4 * V_Q[q1+2] +
                            Z4 * V_Q[q1+3] )); //add fabs for quaternions distance calculation
    double ds = rR*rR + dd;
    if (ds < NNs && ds > 0) { NNs = ds; }
  }
}

// Action_GIST::SumEVV()
void Action_GIST::SumEVV() {
  if (E_VV_VDW_.size() > 1) {
    for (unsigned int gr_pt = 0; gr_pt != MAX_GRID_PT_; gr_pt++) {
      for (unsigned int thread = 1; thread < E_VV_VDW_.size(); thread++) {
        E_UV_VDW_[0][gr_pt]  += E_UV_VDW_[thread][gr_pt];
        E_UV_Elec_[0][gr_pt] += E_UV_Elec_[thread][gr_pt];
        E_VV_VDW_[0][gr_pt]  += E_VV_VDW_[thread][gr_pt];
        E_VV_Elec_[0][gr_pt] += E_VV_Elec_[thread][gr_pt];
        neighbor_[0][gr_pt]  += neighbor_[thread][gr_pt];
      }
    }
  }
}

/** Calculate average voxel energy for PME grids. */
void Action_GIST::CalcAvgVoxelEnergy_PME(double Vvox, DataSet_GridFlt& PME_dens, DataSet_GridFlt& U_PME_dens, Farray& PME_norm)
const
{
  double PME_tot =0.0;
  double U_PME_tot = 0.0;
  mprintf("\t Calculating average voxel energies: \n");
  ProgressBar E_progress(MAX_GRID_PT_);
  for ( unsigned int gr_pt =0; gr_pt < MAX_GRID_PT_; gr_pt++)
  {
    E_progress.Update(gr_pt);
    int nw_total = N_waters_[gr_pt];
    if (nw_total >=1)
    {
      PME_dens[gr_pt] = E_pme_[gr_pt] / (NFRAME_ * Vvox);
      PME_norm[gr_pt] = E_pme_[gr_pt] / nw_total;
      PME_tot += PME_dens[gr_pt];

    }else{
      PME_dens[gr_pt]=0;
      PME_norm[gr_pt]=0; 
    }
    int ns_total = N_solute_atoms_[gr_pt];  
    if (ns_total >=1)
    {
      U_PME_dens[gr_pt] = U_E_pme_[gr_pt] / (NFRAME_ * Vvox);
      U_PME_tot += U_PME_dens[gr_pt];

    }else{
      U_PME_dens[gr_pt]=0;
    }
  }
  PME_tot *=Vvox;
  U_PME_tot *=Vvox;

  infofile_->Printf("Ensemble total water energy on the grid: %9.5f Kcal/mol \n", PME_tot);
  infofile_->Printf("Ensemble total solute energy on the grid: %9.5f Kcal/mol \n",U_PME_tot);

//  infofile_->Printf("Ensemble solute's total potential energy : %9.5f Kcal/mol \n", solute_potential_energy_ / NFRAME_);
//  infofile_->Printf("Ensemble system's total potential energy: %9.5f Kcal/mol \n", system_potential_energy_/NFRAME_);
}

/** Calculate average voxel energy for GIST grids. */
void Action_GIST::CalcAvgVoxelEnergy(double Vvox, DataSet_GridFlt& Eww_dens, DataSet_GridFlt& Esw_dens,
                                     Farray& Eww_norm, Farray& Esw_norm,
                                     DataSet_GridDbl& qtet,
                                     DataSet_GridFlt& neighbor_norm, Farray& neighbor_dens)
{
    #ifndef CUDA
    Darray const& E_UV_VDW = E_UV_VDW_[0];
    Darray const& E_UV_Elec = E_UV_Elec_[0];
    Darray const& E_VV_VDW = E_VV_VDW_[0];
    Darray const& E_VV_Elec = E_VV_Elec_[0];
    #endif
    Farray const& Neighbor = neighbor_[0];
    #ifndef CUDA
    // Sum values from other threads if necessary
    SumEVV();
    #endif
    double Eswtot = 0.0;
    double Ewwtot = 0.0;
    mprintf("\tCalculating average voxel energies:\n");
    ProgressBar E_progress( MAX_GRID_PT_ );
    for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++)
    {
      E_progress.Update( gr_pt );

      //mprintf("DEBUG1: VV vdw=%f elec=%f\n", E_VV_VDW_[gr_pt], E_VV_Elec_[gr_pt]);
      int nw_total = N_waters_[gr_pt]; // Total number of waters that have been in this voxel.
      if (nw_total > 0) {
        #ifndef CUDA
        Esw_dens[gr_pt] = (E_UV_VDW[gr_pt]  + E_UV_Elec[gr_pt]) / (NFRAME_ * Vvox);
        Esw_norm[gr_pt] = (E_UV_VDW[gr_pt]  + E_UV_Elec[gr_pt]) / nw_total;
        Eww_dens[gr_pt] = (E_VV_VDW[gr_pt]  + E_VV_Elec[gr_pt]) / (2 * NFRAME_ * Vvox);
        Eww_norm[gr_pt] = (E_VV_VDW[gr_pt]  + E_VV_Elec[gr_pt]) / (2 * nw_total);
        #else
        double esw = this->Esw_->operator[](gr_pt);
        double eww = this->Eww_->operator[](gr_pt);
        Esw_dens[gr_pt] = esw / (this->NFRAME_ * Vvox);
        Esw_norm[gr_pt] = esw / nw_total;
        Eww_dens[gr_pt] = eww / (this->NFRAME_ * Vvox);
        Eww_norm[gr_pt] = eww / nw_total;
        #endif
        Eswtot += Esw_dens[gr_pt];
        Ewwtot += Eww_dens[gr_pt];

      } else {
        Esw_dens[gr_pt]=0;
        Esw_norm[gr_pt]=0;
        Eww_norm[gr_pt]=0;
        Eww_dens[gr_pt]=0;
      }
      // Compute the average number of water neighbor and average order parameter.
      if (nw_total > 0) {
        qtet[gr_pt] /= nw_total;
        //mprintf("DEBUG1: neighbor= %8.1f  nw_total= %8i\n", neighbor[gr_pt], nw_total);
        neighbor_norm[gr_pt] = (double)Neighbor[gr_pt] / nw_total;
      }
      neighbor_dens[gr_pt] = (double)Neighbor[gr_pt] / (NFRAME_ * Vvox);
    } // END loop over all grid points (voxels)
    Eswtot *= Vvox;
    Ewwtot *= Vvox;
    infofile_->Printf("Total water-solute energy of the grid: Esw = %9.5f kcal/mol\n", Eswtot);
    infofile_->Printf("Total unreferenced water-water energy of the grid: Eww = %9.5f kcal/mol\n",
                      Ewwtot);
}

/** Handle averaging for grids and output from GIST. */
void Action_GIST::Print() {
  gist_print_.Start();
  double Vvox = gO_->Bin().VoxelVolume();

  mprintf("    GIST OUTPUT:\n");

  // The variables are kept outside, so that they are declared for later use.
  // Calculate orientational entropy
  DataSet_GridFlt& dTSorient_dens = static_cast<DataSet_GridFlt&>( *dTSorient_ );
  Farray dTSorient_norm( MAX_GRID_PT_, 0.0 );
  double dTSorienttot = 0;
  int nwtt = 0;
  double dTSo = 0;
  if (! this->skipS_) {
    // LOOP over all voxels
    mprintf("\tCalculating orientational entropy:\n");
    ProgressBar oe_progress( MAX_GRID_PT_ );
    for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
      oe_progress.Update( gr_pt );
      dTSorient_dens[gr_pt] = 0;
      dTSorient_norm[gr_pt] = 0;
      int nw_total = N_waters_[gr_pt]; // Total number of waters that have been in this voxel.
      nwtt += nw_total;
      //mprintf("DEBUG1: %u nw_total %i\n", gr_pt, nw_total);
      if (nw_total > 1) {
        for (int n0 = 0; n0 < nw_total; n0++)
        {
          double NNr = 10000;
          int q0 = n0 * 4; // Index into voxel_Q_ for n0
          for (int n1 = 0; n1 < nw_total; n1++)
          {
            if (n0 != n1) {
              int q1 = n1 * 4; // Index into voxel_Q_ for n1
              //mprintf("DEBUG1:\t\t q1= %8i {%12.4f %12.4f %12.4f %12.4f} q0= %8i {%12.4f %12.4f %12.4f %12.4f}\n",
              //        q1, voxel_Q_[gr_pt][q1  ], voxel_Q_[gr_pt][q1+1], voxel_Q_[gr_pt][q1+2], voxel_Q_[gr_pt][q1+3],
              //        q0, voxel_Q_[gr_pt][q0  ], voxel_Q_[gr_pt][q0+1], voxel_Q_[gr_pt][q0+2], voxel_Q_[gr_pt][q0+3]);
              double rR = 2.0 * acos(  fabs(voxel_Q_[gr_pt][q1  ] * voxel_Q_[gr_pt][q0  ]
                                   + voxel_Q_[gr_pt][q1+1] * voxel_Q_[gr_pt][q0+1]
                                   + voxel_Q_[gr_pt][q1+2] * voxel_Q_[gr_pt][q0+2]
                                   + voxel_Q_[gr_pt][q1+3] * voxel_Q_[gr_pt][q0+3] )); // add fabs for quaternion distance calculation
              //mprintf("DEBUG1:\t\t %8i %8i %g\n", n0, n1, rR);
              if (rR > 0 && rR < NNr) NNr = rR;
            }
          } // END inner loop over all waters for this voxel

          if (NNr < 9999 && NNr > 0) {
            double dbl = log(NNr*NNr*NNr*nw_total / (3.0*Constants::TWOPI));
            //mprintf("DEBUG1: %u  nw_total= %i  NNr= %f  dbl= %f\n", gr_pt, nw_total, NNr, dbl);
            dTSorient_norm[gr_pt] += dbl;
            dTSo += dbl;
          }
        } // END outer loop over all waters for this voxel
        //mprintf("DEBUG1: dTSorient_norm %f\n", dTSorient_norm[gr_pt]);
        dTSorient_norm[gr_pt] = Constants::GASK_KCAL * temperature_ *
                                ((dTSorient_norm[gr_pt]/nw_total) + Constants::EULER_MASC);
        double dtso_norm_nw = (double)dTSorient_norm[gr_pt] * (double)nw_total;
        dTSorient_dens[gr_pt] = (dtso_norm_nw / (NFRAME_ * Vvox));
        dTSorienttot += dTSorient_dens[gr_pt];
        //mprintf("DEBUG1: %f\n", dTSorienttot);
      }
    } // END loop over all grid points (voxels)
    dTSorienttot *= Vvox;
    infofile_->Printf("Maximum number of waters found in one voxel for %d frames = %d\n",
                      NFRAME_, max_nwat_);
    infofile_->Printf("Total referenced orientational entropy of the grid:"
                      " dTSorient = %9.5f kcal/mol, Nf=%d\n", dTSorienttot, NFRAME_);
  }
  // Compute translational entropy for each voxel
  double dTStranstot = 0.0;
  double dTSt = 0.0;
  double dTSs = 0.0;
  int nwts = 0;
  unsigned int nx = gO_->NX();
  unsigned int ny = gO_->NY();
  unsigned int nz = gO_->NZ();
  unsigned int addx = ny * nz;
  unsigned int addy = nz;
  unsigned int addz = 1;
  DataSet_GridFlt& gO = static_cast<DataSet_GridFlt&>( *gO_ );
  DataSet_GridFlt& gH = static_cast<DataSet_GridFlt&>( *gH_ );
  DataSet_GridFlt& dTStrans = static_cast<DataSet_GridFlt&>( *dTStrans_ );
  DataSet_GridFlt& dTSsix = static_cast<DataSet_GridFlt&>( *dTSsix_ );
  Farray dTStrans_norm( MAX_GRID_PT_, 0.0 );
  Farray dTSsix_norm( MAX_GRID_PT_, 0.0 );

  // Loop over all grid points
  if (! this->skipS_)
    mprintf("\tCalculating translational entropy:\n");
  else
    mprintf("Calculating Densities:\n");
  ProgressBar te_progress( MAX_GRID_PT_ );
  for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
    te_progress.Update( gr_pt );
    int numplane = gr_pt / addx;
    double W_dens = 1.0 * N_waters_[gr_pt] / (NFRAME_*Vvox);
    gO[gr_pt] = W_dens / BULK_DENS_;
    gH[gr_pt] = 1.0 * N_hydrogens_[gr_pt] / (NFRAME_*Vvox*2*BULK_DENS_);
    if (! this->skipS_) {
      int nw_total = N_waters_[gr_pt]; // Total number of waters that have been in this voxel.
      for (int n0 = 0; n0 < nw_total; n0++)
      {
        double NNd = 10000;
        double NNs = 10000;
        int i0 = n0 * 3; // index into voxel_xyz_ for n0
        float VX = voxel_xyz_[gr_pt][i0  ];
        float VY = voxel_xyz_[gr_pt][i0+1];
        float VZ = voxel_xyz_[gr_pt][i0+2];
        int q0 = n0 * 4;  // index into voxel_Q_ for n0
        float W4 = voxel_Q_[gr_pt][q0  ];
        float X4 = voxel_Q_[gr_pt][q0+1];
        float Y4 = voxel_Q_[gr_pt][q0+2];
        float Z4 = voxel_Q_[gr_pt][q0+3];
        // First do own voxel
        for (int n1 = 0; n1 < nw_total; n1++) {
          if ( n1 != n0) {
            int i1 = n1 * 3; // index into voxel_xyz_ for n1
            double dx = (double)(VX - voxel_xyz_[gr_pt][i1  ]);
            double dy = (double)(VY - voxel_xyz_[gr_pt][i1+1]);
            double dz = (double)(VZ - voxel_xyz_[gr_pt][i1+2]);
            double dd = dx*dx+dy*dy+dz*dz;
            if (dd < NNd && dd > 0) { NNd = dd; }
            int q1 = n1 * 4; // index into voxel_Q_ for n1
            double rR = 2 * acos( fabs(W4*voxel_Q_[gr_pt][q1  ] +
                                  X4*voxel_Q_[gr_pt][q1+1] +
                                  Y4*voxel_Q_[gr_pt][q1+2] +
                                  Z4*voxel_Q_[gr_pt][q1+3] )); //add fabs for quaternion distance calculation
            double ds = rR*rR + dd;
            if (ds < NNs && ds > 0) { NNs = ds; }
          }
        } // END self loop over all waters for this voxel
        //mprintf("DEBUG1: self NNd=%f NNs=%f\n", NNd, NNs);
        // Determine which directions are possible.
        bool cannotAddZ = (nz == 0 || ( gr_pt%nz == nz-1 ));
        bool cannotAddY = ((nz == 0 || ny-1 == 0) || ( gr_pt%(nz*(ny-1)+(numplane*addx)) < nz));
        bool cannotAddX = (gr_pt >= addx * (nx-1) && gr_pt < addx * nx );
        bool cannotSubZ = (nz == 0 || gr_pt%nz == 0);
        bool cannotSubY = ((nz == 0 || ny == 0) || (gr_pt%addx < nz));
        bool cannotSubX = ((nz == 0 || ny == 0) || (gr_pt < addx));
        bool boundary = ( cannotAddZ || cannotAddY || cannotAddX ||
                          cannotSubZ || cannotSubY || cannotSubX );
        if (!boundary) {
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addy, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addy, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addz + addy, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addz - addy, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addz + addy, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addz - addy, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addz + addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addz - addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addz + addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addz - addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addy + addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addy - addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addy + addx, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addy - addx, NNd, NNs);

          // add the 8 more voxels for NNr searching

          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addx + addy + addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addx + addy - addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addx - addy + addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt + addx - addy - addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addx + addy + addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addx + addy - addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addx - addy + addz, NNd, NNs);
          TransEntropy(VX, VY, VZ, W4, X4, Y4, Z4, gr_pt - addx - addy - addz, NNd, NNs);


          NNd = sqrt(NNd);
          NNs = sqrt(NNs);

          if (NNd < 3 && NNd > 0/*NNd < 9999 && NNd > 0*/) {
            double dbl = log((NNd*NNd*NNd*NFRAME_*4*Constants::PI*BULK_DENS_)/3);
            dTStrans_norm[gr_pt] += dbl;
            dTSt += dbl;
            dbl = log((NNs*NNs*NNs*NNs*NNs*NNs*NFRAME_*Constants::PI*BULK_DENS_)/48);
            dTSsix_norm[gr_pt] += dbl;
            dTSs += dbl;
            //mprintf("DEBUG1: dbl=%f NNs=%f\n", dbl, NNs);
          }
        }
      } // END loop over all waters for this voxel
      if (dTStrans_norm[gr_pt] != 0) {
        nwts += nw_total;
        dTStrans_norm[gr_pt] = Constants::GASK_KCAL*temperature_*( (dTStrans_norm[gr_pt]/nw_total) +
                                                                  Constants::EULER_MASC );
        dTSsix_norm[gr_pt] = Constants::GASK_KCAL*temperature_*( (dTSsix_norm[gr_pt]/nw_total) +
                                                                Constants::EULER_MASC );
      }
      double dtst_norm_nw = (double)dTStrans_norm[gr_pt] * (double)nw_total;
      dTStrans[gr_pt] = (dtst_norm_nw / (NFRAME_*Vvox));
      double dtss_norm_nw = (double)dTSsix_norm[gr_pt] * (double)nw_total;
      dTSsix[gr_pt] = (dtss_norm_nw / (NFRAME_*Vvox));
      dTStranstot += dTStrans[gr_pt];
    } // END loop over all grid points (voxels)
  }
  if (!this->skipS_) {
    dTStranstot *= Vvox;
    double dTSst = 0.0;
    double dTStt = 0.0;
    if (nwts > 0) {
      dTSst = Constants::GASK_KCAL*temperature_*((dTSs/nwts) + Constants::EULER_MASC);
      dTStt = Constants::GASK_KCAL*temperature_*((dTSt/nwts) + Constants::EULER_MASC);
    }
    double dTSot = Constants::GASK_KCAL*temperature_*((dTSo/nwtt) + Constants::EULER_MASC);
    infofile_->Printf("watcount in vol = %d\n", nwtt);
    infofile_->Printf("watcount in subvol = %d\n", nwts);
    infofile_->Printf("Total referenced translational entropy of the grid:"
                      " dTStrans = %9.5f kcal/mol, Nf=%d\n", dTStranstot, NFRAME_);
    infofile_->Printf("Total 6d if all one vox: %9.5f kcal/mol\n", dTSst);
    infofile_->Printf("Total t if all one vox: %9.5f kcal/mol\n", dTStt);
    infofile_->Printf("Total o if all one vox: %9.5f kcal/mol\n", dTSot);
  }
  // Compute average voxel energy. Allocate these sets even if skipping energy
  // to be consistent with previous output.
  DataSet_GridFlt& PME_dens = static_cast<DataSet_GridFlt&>( *PME_);
  DataSet_GridFlt& U_PME_dens = static_cast<DataSet_GridFlt&>( *U_PME_);
  DataSet_GridFlt& Esw_dens = static_cast<DataSet_GridFlt&>( *Esw_ );
  DataSet_GridFlt& Eww_dens = static_cast<DataSet_GridFlt&>( *Eww_ );
  DataSet_GridFlt& neighbor_norm = static_cast<DataSet_GridFlt&>( *neighbor_norm_ );
  DataSet_GridDbl& qtet = static_cast<DataSet_GridDbl&>( *order_norm_ );
  Farray Esw_norm( MAX_GRID_PT_, 0.0 );
  Farray Eww_norm( MAX_GRID_PT_, 0.0 );
  Farray PME_norm( MAX_GRID_PT_,0.0);
  Farray neighbor_dens( MAX_GRID_PT_, 0.0 );
  if (!skipE_) {
    if (usePme_) {
      CalcAvgVoxelEnergy_PME(Vvox, PME_dens, U_PME_dens, PME_norm);
    }// else {
      CalcAvgVoxelEnergy(Vvox, Eww_dens, Esw_dens, Eww_norm, Esw_norm, qtet,
                         neighbor_norm, neighbor_dens);
    //}
  }
  // Compute average dipole density.
  DataSet_GridFlt& pol = static_cast<DataSet_GridFlt&>( *dipole_ );
  DataSet_GridDbl& dipolex = static_cast<DataSet_GridDbl&>( *dipolex_ );
  DataSet_GridDbl& dipoley = static_cast<DataSet_GridDbl&>( *dipoley_ );
  DataSet_GridDbl& dipolez = static_cast<DataSet_GridDbl&>( *dipolez_ );
  for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++)
  {
    dipolex[gr_pt] /= (Constants::DEBYE_EA * NFRAME_ * Vvox);
    dipoley[gr_pt] /= (Constants::DEBYE_EA * NFRAME_ * Vvox);
    dipolez[gr_pt] /= (Constants::DEBYE_EA * NFRAME_ * Vvox);
    pol[gr_pt] = sqrt( dipolex[gr_pt]*dipolex[gr_pt] +
                       dipoley[gr_pt]*dipoley[gr_pt] +
                       dipolez[gr_pt]*dipolez[gr_pt] );
  }

  // Write the GIST output file.
  // TODO: Make a data file format?
  if (datafile_ != 0) {
    mprintf("\tWriting GIST results for each voxel:\n");
    // Create the format strings.
    std::string fmtstr =
                intFmt_.Fmt() + // grid point
          " " + fltFmt_.Fmt() + // grid X
          " " + fltFmt_.Fmt() + // grid Y
          " " + fltFmt_.Fmt() + // grid Z
          " " + intFmt_.Fmt() + // # waters
          " " + fltFmt_.Fmt() + // gO
          " " + fltFmt_.Fmt() + // gH
          " " + fltFmt_.Fmt() + // dTStrans
          " " + fltFmt_.Fmt() + // dTStrans_norm
          " " + fltFmt_.Fmt() + // dTSorient_dens
          " " + fltFmt_.Fmt() + // dTSorient_norm
          " " + fltFmt_.Fmt() + // dTSsix
          " " + fltFmt_.Fmt() + // dTSsix_norm
          " " + fltFmt_.Fmt() + // Esw_dens
          " " + fltFmt_.Fmt() + // Esw_norm
          " " + fltFmt_.Fmt() + // Eww_dens
          " " + fltFmt_.Fmt();  // EWW_norm
    if (usePme_) {
      fmtstr +=
          " " + fltFmt_.Fmt() + // PME_dens
        + " " + fltFmt_.Fmt();  // PME_norm
    }
    fmtstr +=
          " " + fltFmt_.Fmt() + // dipolex
          " " + fltFmt_.Fmt() + // dipoley
          " " + fltFmt_.Fmt() + // dipolez
          " " + fltFmt_.Fmt() + // pol
          " " + fltFmt_.Fmt() + // neighbor_dens
          " " + fltFmt_.Fmt() + // neighbor_norm
          " " + fltFmt_.Fmt() + // qtet
          " \n";                // NEWLINE
    if (debug_ > 0) mprintf("DEBUG: Fmt='%s'\n", fmtstr.c_str());
    const char* gistOutputVersion;
    if (usePme_)
      gistOutputVersion = "v3";
    else
      gistOutputVersion = "v2";
    // Do the header
    datafile_->Printf("GIST Output %s "
                      "spacing=%.4f center=%.6f,%.6f,%.6f dims=%i,%i,%i \n"
                      "voxel xcoord ycoord zcoord population g_O g_H"
                      " dTStrans-dens(kcal/mol/A^3) dTStrans-norm(kcal/mol)"
                      " dTSorient-dens(kcal/mol/A^3) dTSorient-norm(kcal/mol)"
                      " dTSsix-dens(kcal/mol/A^3) dTSsix-norm(kcal/mol)"
                      " Esw-dens(kcal/mol/A^3) Esw-norm(kcal/mol)"
                      " Eww-dens(kcal/mol/A^3) Eww-norm-unref(kcal/mol)",
                      gistOutputVersion, gridspacing_,
                      gridcntr_[0], gridcntr_[1], gridcntr_[2],
                      (int)griddim_[0], (int)griddim_[1], (int)griddim_[2]);
    if (usePme_)
      datafile_->Printf(" PME-dens(kcal/mol/A^3) PME-norm(kcal/mol)");
    datafile_->Printf(" Dipole_x-dens(D/A^3) Dipole_y-dens(D/A^3) Dipole_z-dens(D/A^3)"
                      " Dipole-dens(D/A^3) neighbor-dens(1/A^3) neighbor-norm order-norm\n");
    // Loop over voxels
    ProgressBar O_progress( MAX_GRID_PT_ );
    for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
      O_progress.Update( gr_pt );
      size_t i, j, k;
      gO_->ReverseIndex( gr_pt, i, j, k );
      Vec3 XYZ = gO_->Bin().Center( i, j, k );
      
      
      if (usePme_) {
        datafile_->Printf(fmtstr.c_str(),
                          gr_pt, XYZ[0], XYZ[1], XYZ[2], N_waters_[gr_pt], gO[gr_pt], gH[gr_pt],
                          dTStrans[gr_pt], dTStrans_norm[gr_pt],
                          dTSorient_dens[gr_pt], dTSorient_norm[gr_pt],
                          dTSsix[gr_pt], dTSsix_norm[gr_pt],
                          Esw_dens[gr_pt], Esw_norm[gr_pt],
                          Eww_dens[gr_pt], Eww_norm[gr_pt],
                          PME_dens[gr_pt], PME_norm[gr_pt],
                          dipolex[gr_pt], dipoley[gr_pt], dipolez[gr_pt],
                          pol[gr_pt], neighbor_dens[gr_pt], neighbor_norm[gr_pt], qtet[gr_pt]);
      } else {
        datafile_->Printf(fmtstr.c_str(),
                          gr_pt, XYZ[0], XYZ[1], XYZ[2], N_waters_[gr_pt], gO[gr_pt], gH[gr_pt],
                          dTStrans[gr_pt], dTStrans_norm[gr_pt],
                          dTSorient_dens[gr_pt], dTSorient_norm[gr_pt],
                          dTSsix[gr_pt], dTSsix_norm[gr_pt],
                          Esw_dens[gr_pt], Esw_norm[gr_pt],
                          Eww_dens[gr_pt], Eww_norm[gr_pt],
                          dipolex[gr_pt], dipoley[gr_pt], dipolez[gr_pt],
                          pol[gr_pt], neighbor_dens[gr_pt], neighbor_norm[gr_pt], qtet[gr_pt]);
      }
    } // END loop over voxels
  } // END datafile_ not null

  // Write water-water interaction energy matrix
  if (ww_Eij_ != 0) {
    DataSet_MatrixFlt& ww_Eij = static_cast<DataSet_MatrixFlt&>( *ww_Eij_ );
    double fac = 1.0 / (double)(NFRAME_ * 2);
    for (unsigned int idx = 0; idx != ww_Eij.Size(); idx++) {
      if (fabs(ww_Eij[idx]) < Constants::SMALL)
        ww_Eij[idx] = 0.0;
      else {
        double val = (double)ww_Eij[idx];
        ww_Eij[idx] = (float)(val * fac);
      }
    }
    // Eij matrix output, skip any zeros.
    for (unsigned int a = 1; a < MAX_GRID_PT_; a++) {
      for (unsigned int l = 0; l < a; l++) {
        double dbl = ww_Eij_->GetElement(a, l);
        if (dbl != 0)
          eijfile_->Printf("%10d %10d %12.5E\n", a, l, dbl);
      }
    }
  }
  gist_print_.Stop();
  double total = gist_init_.Total() + gist_setup_.Total() +
                 gist_action_.Total() + gist_print_.Total();
  mprintf("\tGIST timings:\n");
  gist_init_.WriteTiming(1,    "Init:  ", total);
  gist_setup_.WriteTiming(1,   "Setup: ", total);
  gist_action_.WriteTiming(1,  "Action:", total);
  gist_grid_.WriteTiming(2,    "Grid:   ", gist_action_.Total());
  gist_nonbond_.WriteTiming(2, "Nonbond:", gist_action_.Total());
# ifdef LIBPME
  if (usePme_)
    gistPme_.Timing( gist_nonbond_.Total() );
# endif
  //gist_nonbond_dist_.WriteTiming(3, "Dist2:", gist_nonbond_.Total());
  //gist_nonbond_UV_.WriteTiming(3, "UV:", gist_nonbond_.Total());
  //gist_nonbond_VV_.WriteTiming(3, "VV:", gist_nonbond_.Total());
  //gist_nonbond_OV_.WriteTiming(3, "OV:", gist_nonbond_.Total());
  gist_euler_.WriteTiming(2,   "Euler:  ", gist_action_.Total());
  gist_dipole_.WriteTiming(2,  "Dipole: ", gist_action_.Total());
  gist_order_.WriteTiming(2,   "Order: ", gist_action_.Total());
  gist_print_.WriteTiming(1,   "Print:", total);
  mprintf("TIME:\tTotal: %.4f s\n", total);
  #ifdef CUDA
  this->freeGPUMemory();
  #endif
}

#ifdef CUDA
void Action_GIST::NonbondCuda(ActionFrame frm) {
  // Simply to get the information for the energetic calculations
  std::vector<float> eww_result(this->numberAtoms_);
  std::vector<float> esw_result(this->numberAtoms_);
  std::vector<std::vector<int> > order_indices;
  this->gist_nonbond_.Start();

  float *recip = NULL;
  float *ucell = NULL;
  int boxinfo;

  // Check Boxinfo and write the necessary data into recip, ucell and boxinfo.
  switch(imageOpt_.ImagingType()) {
    case ImageOption::NONORTHO:
      recip = new float[9];
      ucell = new float[9];
      for (int i = 0; i < 9; ++i) {
        ucell[i] = (float) frm.Frm().BoxCrd().UnitCell()[i];
        recip[i] = (float) frm.Frm().BoxCrd().FracCell()[i];
      }
      boxinfo = 2;
      break;
    case ImageOption::ORTHO:
      recip = new float[9];
      recip[0] = frm.Frm().BoxCrd().Param(Box::X);
      recip[1] = frm.Frm().BoxCrd().Param(Box::Y);
      recip[2] = frm.Frm().BoxCrd().Param(Box::Z);
      ucell = NULL;
      boxinfo = 1;
      break;
    case ImageOption::NO_IMAGE:
      recip = NULL;
      ucell = NULL;
      boxinfo = 0;
      break;
    default:
      mprinterr("Error: Unexpected box information found.");
      return;
  }

  std::vector<int> result_o = std::vector<int>(4 * this->numberAtoms_);
  std::vector<int> result_n = std::vector<int>(this->numberAtoms_);
  // Call the GPU Wrapper, which subsequently calls the kernel, after setup operations.
  // Must create arrays from the vectors, does that by getting the address of the first element of the vector.
  std::vector<std::vector<float> > e_result = doActionCudaEnergy(frm.Frm().xAddress(), this->NBindex_c_, this->numberAtomTypes_, this->paramsLJ_c_, this->molecule_c_, boxinfo, recip, ucell, this->numberAtoms_, this->min_c_,
                                                    this->max_c_, this->headAtomType_,this->NeighborCut2_, &(result_o[0]), &(result_n[0]), this->result_w_c_,
                                                    this->result_s_c_, this->result_O_c_, this->result_N_c_, this->doOrder_);
  eww_result = e_result.at(0);
  esw_result = e_result.at(1);

  if (this->doOrder_) {
    int counter = 0;
    for (unsigned int i = 0; i < (4 * this->numberAtoms_); i += 4) {
      ++counter;
      std::vector<int> temp;
      for (unsigned int j = 0; j < 4; ++j) {
        temp.push_back(result_o.at(i + j));
      }
      order_indices.push_back(temp);
    }
  }

  delete[] recip; // Free memory
  delete[] ucell; // Free memory

  for (unsigned int sidx = 0; sidx < NSOLVENT_; sidx++) {
    int headAtomIndex = O_idxs_[sidx];
    size_t bin_i, bin_j, bin_k;
    const double *vec = frm.Frm().XYZ(headAtomIndex);
    int voxel = -1;
    if (this->gO_->Bin().Calc(vec[0], vec[1], vec[2], bin_i, bin_j, bin_k)) {
      voxel = this->gO_->CalcIndex(bin_i, bin_j, bin_k);
      this->neighbor_.at(0).at(voxel) += result_n.at(headAtomIndex);
      // This is not nice, as it assumes that O is set before the two Hydrogens
      // might be the case, but is still not nice (in my opinion)
      for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
        this->Esw_->UpdateVoxel(voxel, esw_result.at(headAtomIndex + IDX));
        this->Eww_->UpdateVoxel(voxel, eww_result.at(headAtomIndex + IDX));
      }
      // Order calculation
      if (this->doOrder_) {
        double sum = 0;
        Vec3 cent( frm.Frm().xAddress() + (headAtomIndex) * 3 );
        std::vector<Vec3> vectors;
        switch(imageOpt_.ImagingType()) {
          case ImageOption::NONORTHO:
          case ImageOption::ORTHO:
            {
              Vec3 vec(frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(0) * 3));
              vectors.push_back( MinImagedVec(vec, cent, frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell()));
              vec = Vec3(frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(1) * 3));
              vectors.push_back( MinImagedVec(vec, cent, frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell()));
              vec = Vec3(frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(2) * 3));
              vectors.push_back( MinImagedVec(vec, cent, frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell()));
              vec = Vec3(frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(3) * 3));
              vectors.push_back( MinImagedVec(vec, cent, frm.Frm().BoxCrd().UnitCell(), frm.Frm().BoxCrd().FracCell()));
            }
            break;
          default:
            vectors.push_back( Vec3( frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(0) * 3) ) - cent );
            vectors.push_back( Vec3( frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(1) * 3) ) - cent );
            vectors.push_back( Vec3( frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(2) * 3) ) - cent );
            vectors.push_back( Vec3( frm.Frm().xAddress() + (order_indices.at(headAtomIndex).at(3) * 3) ) - cent );
        }

        for (int i = 0; i < 3; ++i) {
          for (int j = i + 1; j < 4; ++j) {
            double cosThet = (vectors.at(i) * vectors.at(j)) / sqrt(vectors.at(i).Magnitude2() * vectors.at(j).Magnitude2());
            sum += (cosThet + 1.0/3) * (cosThet + 1.0/3);
          }
        }
        this->order_norm_->UpdateVoxel(voxel, 1.0 - (3.0/8.0) * sum);
      }
    }

  }
  this->gist_nonbond_.Stop();
}

/**
 * Frees all the Memory on the GPU.
 */
void Action_GIST::freeGPUMemory(void) {
  freeCuda(this->NBindex_c_);
  freeCuda(this->molecule_c_);
  freeCuda(this->paramsLJ_c_);
  freeCuda(this->max_c_);
  freeCuda(this->min_c_);
  freeCuda(this->result_w_c_);
  freeCuda(this->result_s_c_);
  freeCuda(this->result_O_c_);
  freeCuda(this->result_N_c_);
  this->NBindex_c_   = NULL;
  this->molecule_c_  = NULL;
  this->paramsLJ_c_  = NULL;
  this->max_c_     = NULL;
  this->min_c_     = NULL;
  this->result_w_c_= NULL;
  this->result_s_c_= NULL;
  this->result_O_c_  = NULL;
  this->result_N_c_  = NULL;
}

/**
 * Copies data from the CPU to the GPU.
 * @throws: CudaException
 */
void Action_GIST::copyToGPU(void) {
  try {
    copyMemoryToDevice(&(this->NBIndex_[0]), this->NBindex_c_, this->NBIndex_.size() * sizeof(int));
    copyMemoryToDeviceStruct(&(this->charges_[0]), &(this->atomTypes_[0]), this->solvent_, &(this->molecule_[0]), this->numberAtoms_, &(this->molecule_c_),
                              &(this->lJParamsA_[0]), &(this->lJParamsB_[0]), this->lJParamsA_.size(), &(this->paramsLJ_c_));
  } catch (CudaException &ce) {
    this->freeGPUMemory();
    mprinterr("Error: Could not copy data to the device.\n");
    throw ce;
  } catch (std::exception &e) {
    this->freeGPUMemory();
    throw e;
  }
}
#endif
