#include <cmath>
#include <cfloat> // DBL_MAX
#include <algorithm>
#include "Action_GIST.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_MatrixFlt.h"
#include "DataSet_GridFlt.h"
#include "DataSet_GridDbl.h"
#include "ProgressBar.h"
#include "StringRoutines.h"
#include "DistRoutines.h"
#include "GistEntropyUtils.h"
#ifdef _OPENMP
# include <omp.h>
#endif

// Note: The Order calculation is not updated for solvents other than water.
// E.g., it does not use rigidAtomIndices[0].
// It will not crash, but also not produce useful results.
//
// The neighbor calculation is also not updated. The distance cutoff is hard-coded to 3.5 
// Angstrom, which is reasonable only for water.
//
// In the CUDA kernel, the neighbor calc uses the atom type to choose which atom to use.
// This is not reasonable with non-water solvents.

// Crashes with zero solvent molecules

// TO-DO: check for linear solvent molecules

const double Action_GIST::maxD_ = DBL_MAX;
#define GIST_TINY 1e-10

std::vector<std::string> split_string(std::string s, const std::string& delimiter)
{
  std::vector<std::string> ret;
  if (s.size() == 0) {
    return ret;
  }
  size_t pos = 0;
  while ((pos = s.find(delimiter)) != std::string::npos) {
      ret.emplace_back(s.substr(0, pos));
      s.erase(0, pos + delimiter.length());
  }
  ret.emplace_back(s);
  return ret;
}

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
  griddim_(),
  rigidAtomNames_(3),
  Esw_(0),
  Eww_(0),
  dTStrans_(0),
  dTSorient_(0),
  dTSsix_(0),
  neighbor_(0),
  dipole_(0),
  order_(0),
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
  exactNnVolume_(false),
  useCom_(true),
  setupSuccessful_(false)
{}

/** GIST help */
void Action_GIST::Help() const {
  mprintf("\t[doorder] [doeij] [skipE] [skipS] [refdens <rdval>] [temp <tval>]\n"
          "\t[noimage] [gridcntr <xval> <yval> <zval>]\n"
          "\t[griddim <nx> <ny> <nz>] [gridspacn <spaceval>] [neighborcut <ncut>]\n"
          "\t[prefix <filename prefix>] [ext <grid extension>] [out <output suffix>]\n"
          "\t[floatfmt {double|scientific|general}] [floatwidth <fw>] [floatprec <fp>]\n"
          "\t[intwidth <iw>] [oldnnvolume] [nnsearchlayers <nlayers>] [solute <mask>] [solventmols <str>]\n"
          "\t[rigidatoms <i1> <i2> <i3> [nocom]\n"
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
  DFL_ = &init.DFL();
  DSL_ = init.DslPtr();
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'gist' action does not work with > 1 process (%i processes currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  gist_init_.Start();
  prefix_ = actionArgs.GetStringKey("prefix", "gist");
  ext_ = actionArgs.GetStringKey("ext", ".dx");
  std::string gistout = actionArgs.GetStringKey("out", prefix_ + "-output.dat");
  datafile_ = init.DFL().AddCpptrajFile( gistout, "GIST output" );
  if (datafile_ == 0) return Action::ERR;
  // Info file: if not specified use STDOUT
  std::string info = actionArgs.GetStringKey("info");
  if (!info.empty()) info = prefix_ + "-" + info;
  infofile_ = init.DFL().AddCpptrajFile( info, "GIST info", DataFileList::TEXT, true );
  if (infofile_ == 0) return Action::ERR;

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
  exactNnVolume_ = !actionArgs.hasKey("oldnnvolume");
  nNnSearchLayers_ = actionArgs.getKeyInt("nnsearchlayers", 1);
  imageOpt_.InitImaging( !(actionArgs.hasKey("noimage")), actionArgs.hasKey("nonortho") );
  doOrder_ = actionArgs.hasKey("doorder");
  doEij_ = actionArgs.hasKey("doeij");
  useCom_ = !actionArgs.hasKey("nocom");
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
  ArgList centerArgs = actionArgs.GetNstringKey("gridcntr", 3);
  if ( centerArgs.empty() ) {
    mprintf("Warning: No grid center values specified, using default (origin)\n");
    gridcntr_ = Vec3(0.0);
  } else {
    if ( !validDouble(centerArgs[0]) || !validDouble(centerArgs[1]) || !validDouble(centerArgs[2]) ) {
      mprinterr("Invalid grid center: %s %s %s\n", centerArgs[0].c_str(), centerArgs[1].c_str(), centerArgs[2].c_str());
      return Action::ERR;
    }
    gridcntr_[0] = centerArgs.getNextDouble(-1);
    gridcntr_[1] = centerArgs.getNextDouble(-1);
    gridcntr_[2] = centerArgs.getNextDouble(-1);
  }
  // Grid dimensions
  ArgList dimArgs = actionArgs.GetNstringKey("griddim", 3);
  if ( dimArgs.empty() ) {
    griddim_[0] = 40;
    griddim_[1] = 40;
    griddim_[2] = 40;
    mprintf("Warning: No grid dimension values specified, using default (40,40,40)\n");
  } else {
    if ( !validInteger(dimArgs[0]) || !validInteger(dimArgs[1]) || !validInteger(dimArgs[2]) ) {
      mprinterr("Invalid grid dimensions: %s %s %s\n", dimArgs[0].c_str(), dimArgs[1].c_str(), dimArgs[2].c_str());
      return Action::ERR;
    }
    griddim_[0] = dimArgs.getNextInteger(-1);
    griddim_[1] = dimArgs.getNextInteger(-1);
    griddim_[2] = dimArgs.getNextInteger(-1);
  }
  if ( griddim_[0] < 1 || griddim_[1] < 1 || griddim_[2] < 1 ) {
    mprinterr("Error: grid dimensions must be >0, but are %d %d %d.\n", griddim_[0], griddim_[1], griddim_[2]);
    return Action::ERR;
  }
  ArgList indArgs = actionArgs.GetNstringKey("rigidatoms", 3);
  if ( !indArgs.empty() ) {
    rigidAtomNames_[0] = indArgs.GetStringNext();
    rigidAtomNames_[1] = indArgs.GetStringNext();
    rigidAtomNames_[2] = indArgs.GetStringNext();
  }
  soluteMask_ = actionArgs.GetStringKey("solute", "");
  solventNames_ = split_string(actionArgs.GetStringKey("solventmols"), ",");
  // Data set name
  dsname_ = actionArgs.GetStringKey("name");
  if (dsname_.empty())
    dsname_ = init.DSL().GenerateDefaultName("GIST");

  // Set up DataSets.

  Esw_ = AddDatasetAndFile("Esw", prefix_ + "-Esw-dens" + ext_, DataSet::GRID_FLT);
  Eww_ = AddDatasetAndFile("Eww", prefix_ + "-Eww-dens" + ext_, DataSet::GRID_FLT);
  dTStrans_ = AddDatasetAndFile("dTStrans", prefix_ + "-dTStrans-dens" + ext_, DataSet::GRID_FLT);
  dTSorient_ = AddDatasetAndFile("dTSorient", prefix_ + "-dTSorient-dens" + ext_, DataSet::GRID_FLT);
  dTSsix_ = AddDatasetAndFile("dTSsix", prefix_ + "-dTSsix-dens" + ext_, DataSet::GRID_FLT);
  neighbor_ = AddDatasetAndFile("neighbor", prefix_ + "-neighbor-norm" + ext_, DataSet::GRID_FLT);
  dipole_ = AddDatasetAndFile("dipole", prefix_ + "-dipole-dens" + ext_, DataSet::GRID_FLT);
  order_ = AddDatasetAndFile("order", prefix_ + "-order-norm" + ext_, DataSet::GRID_DBL);
  dipolex_ = AddDatasetAndFile("dipolex", prefix_ + "-dipolex-dens" + ext_, DataSet::GRID_DBL);
  dipoley_ = AddDatasetAndFile("dipoley", prefix_ + "-dipoley-dens" + ext_, DataSet::GRID_DBL);
  dipolez_ = AddDatasetAndFile("dipolez", prefix_ + "-dipolez-dens" + ext_, DataSet::GRID_DBL);

  if (!Esw_ || !Eww_ || !dTStrans_ || !dTSorient_ || !dTSsix_ || !neighbor_ || !dipole_ || !order_
      || !dipolex_ || !dipoley_ || !dipolez_) {
    return Action::ERR;
  }

  if (usePme_) {
    PME_ = AddDatasetAndFile("PME", prefix_ + "-Water-Etot-pme-dens" + ext_, DataSet::GRID_DBL);
    U_PME_ = AddDatasetAndFile("U_PME", prefix_ + "-Solute-Etot-pme-dens"+ ext_, DataSet::GRID_DBL);

    if (!PME_ || !U_PME_) {
      return Action::ERR;
    }
  }

  if (!createMoleculeDatasets()) {
    mprinterr("Failed to create datasets for molecular densities;\n");
    return Action::ERR;
  }

  gridBin_ = &Eww_->Bin();

  if (doEij_) {
    ww_Eij_ = (DataSet_MatrixFlt*)init.DSL().AddSet(DataSet::MATRIX_FLT, MetaData(dsname_, "Eij"));
    if (ww_Eij_ == 0) return Action::ERR;
  }

  MAX_GRID_PT_ = griddim_[0] * griddim_[1] * griddim_[2];

  if (ww_Eij_ != 0) {
    if (ww_Eij_->AllocateTriangle( MAX_GRID_PT_ )) {
      mprinterr("Error: Could not allocate memory for water-water Eij matrix.\n");
      return Action::ERR;
    }
  }

  // Set up grid params TODO non-orthogonal as well
  G_max_ = Vec3( (double)griddim_[0] * gridspacing_ + 1.5,
                 (double)griddim_[1] * gridspacing_ + 1.5,
                 (double)griddim_[2] * gridspacing_ + 1.5 );
  N_solvent_.assign( MAX_GRID_PT_, 0 );
  N_main_solvent_.assign( MAX_GRID_PT_, 0 );
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
    neighborPerThread_.resize( numthreads_ );
    for (int thread = 0; thread != numthreads_; thread++) {
      E_UV_VDW_[thread].assign( MAX_GRID_PT_, 0 );
      E_UV_Elec_[thread].assign( MAX_GRID_PT_, 0 );
      E_VV_VDW_[thread].assign( MAX_GRID_PT_, 0 );
      E_VV_Elec_[thread].assign( MAX_GRID_PT_, 0 );
      neighborPerThread_[thread].assign( MAX_GRID_PT_, 0 );
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
  mprintf("\tOutput prefix= '%s', grid output extension= '%s'\n", prefix_.c_str(), ext_.c_str());
  mprintf("\tOutput float format string= '%s', output integer format string= '%s'\n", fltFmt_.fmt(), intFmt_.fmt());
  mprintf("\tGIST info written to '%s'\n", infofile_->Filename().full());
  mprintf("\tName for data sets: %s\n", dsname_.c_str());
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
  Esw_->GridInfo();
  mprintf("\tNumber of voxels: %u, voxel volume: %f Ang^3\n",
          MAX_GRID_PT_, gridBin_->VoxelVolume());
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

/**
 * Adds a dataset to the global DataSetList, and a datafile to the global DataFileList. Stores them in dataSets3D_.
 *
 * The DataType MUST be a subclass of DataSet_3D (GRID_FLT, GRID_DBL, GRID_...)!!!
 * \param name key in dataSets3D_
 * \param filename name of the output file (usuallz ...dx)
 * \param dtype DataSet::GRID_DBL or another GRID_ datatype.
 */
DataSet_3D* Action_GIST::AddDatasetAndFile(const std::string& name, const std::string& filename, DataSet::DataType dtype)
{
  DataFile* file = DFL_->AddDataFile( filename );
  DataSet_3D* dataset = static_cast<DataSet_3D*>(DSL_->AddSet(dtype, MetaData(dsname_, name)));
  if (dataset == NULL || file == NULL) {
    return 0;
  }
  Vec3 v_spacing( gridspacing_ );
  dataset->Allocate_N_C_D(griddim_[0], griddim_[1], griddim_[2], gridcntr_, v_spacing);
  file->AddDataSet( dataset );
  return dataset;
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

  // NOTE: these are just guesses
  O_idxs_.reserve( setup.Top().Nsolvent() );
  atomIsSolute_.assign(setup.Top().Natom(), false);
  atomIsSolventO_.assign(setup.Top().Natom(), false);
  U_idxs_.reserve(setup.Top().Natom()-setup.Top().Nsolvent()*nMolAtoms_);

  setSolventType(setup.Top());

  setSoluteSolvent(setup.Top());

  bool isFirstSolvent = true;
  for (Topology::mol_iterator mol = setup.Top().MolStart();
                              mol != setup.Top().MolEnd(); ++mol)
  {
    int mol_start = mol->MolUnit().Front();
    if (!atomIsSolute_[mol_start]) {
      O_idxs_.push_back( mol_start );
    }
    if (isMainSolvent(mol_start)) {
      int error;
      if (isFirstSolvent) {
        error = setSolventProperties(*mol, setup.Top());
        analyzeSolventElements(*mol, setup.Top());
        if (!setRigidAtomIndices(*mol, setup.Top())) { mprinterr("Failed to set indices of rigid atoms.\n"); error = 1; }
        if (rigidAtomIndices_[0] < 0 || rigidAtomIndices_[0] >= (int)nMolAtoms_
            || rigidAtomIndices_[1] < 0 || rigidAtomIndices_[1] >= (int)nMolAtoms_
            || rigidAtomIndices_[2] < 0 || rigidAtomIndices_[2] >= (int)nMolAtoms_) {
          mprinterr("All rigidatomindices must be between 1 and the number of atoms per solvent.\n");
          error = 1;
        }
        if (!createAtomDensityDatasets()) { mprinterr("Failed to create datasets for atomic densities.\n"); error = 1; }
        isFirstSolvent = false;
      } else {
        error = checkSolventProperties(*mol, setup.Top());
      }
      if (error != 0) {
        mprinterr("Error: In molecule %s.\n",
                  setup.Top().TruncResNameNum( setup.Top()[mol_start].ResNum() ).c_str());
        return Action::ERR;
      }
      atomIsSolventO_[mol_start+rigidAtomIndices_[0]] = true;
    }
  }
  #ifdef CUDA
  for (int i = 0; i != setup.Top().Natom(); ++i) {
    this->molecule_.push_back( setup.Top()[i].MolNum() );
    this->charges_.push_back( setup.Top()[i].Charge() );
    this->atomTypes_.push_back( setup.Top()[i].TypeIndex() );
  }
  #endif
  NSOLVENT_ = O_idxs_.size();
  int NsolventAtoms = NSOLVENT_ * nMolAtoms_;
  mprintf("\t%u solvent molecules, %u solvent atoms, %zu solute atoms (%d total).\n",
          NSOLVENT_, NsolventAtoms, U_idxs_.size(), setup.Top().Natom());
  if (doOrder_ && NSOLVENT_ < 5) {
    mprintf("Warning: Less than 5 solvent molecules. Cannot perform order calculation.\n");
    doOrder_ = false;
  }
  // Allocate space for saving indices of water atoms that are on the grid
  // Estimate how many solvent molecules can possibly fit onto the grid.
  // Add some extra voxels as a buffer.
  double max_voxels = 1.10 * (double)MAX_GRID_PT_;
  double totalVolume = max_voxels * gridBin_->VoxelVolume();
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
    mprinterr("Error: Could not copy memory to GPU!\n");
    return Action::ERR;
  }
#endif

  /* for (int i = 0; i < setup.Top().Natom(); ++i) { */
  /*   mprintf("%d\n", solventType_[i]); */
  /* } */
  setupSuccessful_ = true;
  gist_setup_.Stop();
  return Action::OK;
}

int Action_GIST::setSolventProperties(const Molecule& mol, const Topology& top)
{
  int o_idx = mol.MolUnit().Front();
  nMolAtoms_ = mol.NumAtoms();
  mprintf("\tEach solvent molecule has %u atoms\n", nMolAtoms_);
  if (top[o_idx].Element() != Atom::OXYGEN ||
      top[o_idx+1].Element() != Atom::HYDROGEN ||
      top[o_idx+2].Element() != Atom::HYDROGEN)
  {
    mprintf("\tFirst solvent molecule '%s' is not water.\n",
            top.TruncResNameNum( top[o_idx].ResNum() ).c_str());
    mprintf("Warning: The neighbor and order calculations are not meant to be used with"
            " solvents other than water!\n");
  }
  #ifdef CUDA
  this->headAtomType_ = top[o_idx].TypeIndex();
  #endif
  double q_sum = 0.0;
  Q_.reserve( nMolAtoms_ );
  for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
    Q_.push_back( top[o_idx+IDX].Charge() );
    q_sum += Q_.back();
  }
  // Sanity checks.
  if (fabs( q_sum ) > 0.0) {
    mprintf("Warning: Charges on solvent do not sum to 0 (%g)\n", q_sum);
  }
  return 0;
}

int Action_GIST::checkSolventProperties(const Molecule& mol, const Topology& top) const
{
  int o_idx = mol.MolUnit().Front();
  if (mol.NumAtoms() != nMolAtoms_) {
    mprinterr("Error: All solvent molecules must have same # atoms.\n"
              "Error: A Molecule has %u atoms, expected %u.\n",
              mol.NumAtoms(), nMolAtoms_);
    return 1;
  }
  for (unsigned int IDX = 0; IDX < nMolAtoms_; IDX++) {
    double q_atom = top[o_idx+IDX].Charge();
    if (NotEqual(Q_[IDX], q_atom)) {
      mprintf("Warning: Charge on water '%s' (%g) does not match first water (%g).\n",
            top.TruncResAtomName( o_idx+IDX ).c_str(), q_atom, Q_[IDX]);
    }
  }
  return 0;
}

void Action_GIST::analyzeSolventElements(const Molecule& mol, const Topology& top)
{
  int o_idx = mol.MolUnit().Front();
  for (unsigned int i_mol=0; i_mol<nMolAtoms_; ++i_mol) {
    std::string elem = top[o_idx+i_mol].ElementName();
    bool found_element = false;
    for (unsigned int i_elem = 0; i_elem < solventInfo_.unique_elements.size(); ++i_elem) {
      if (solventInfo_.unique_elements[i_elem] == elem) {
        solventInfo_.i_element.push_back(i_elem);
        ++solventInfo_.element_count[i_elem];
        found_element = true;
        break;
      }
    }
    if (!found_element) {
      solventInfo_.unique_elements.push_back(elem);
      solventInfo_.i_element.push_back(solventInfo_.unique_elements.size() - 1);
      solventInfo_.element_count.push_back(1);
    }
  }
}

bool Action_GIST::setRigidAtomIndices(const Molecule& mol, const Topology& top)
{
  unsigned int mol_begin = mol.MolUnit().Front();
  unsigned int mol_end = mol.MolUnit().Back();

  if (rigidAtomNames_[0].length() == 0 || rigidAtomNames_[1].length() == 0 || rigidAtomNames_[2].length() == 0) {
    // Set rigidAtomIndices_ automatically
    std::vector<int> n_bonds_heavy(mol_end - mol_begin);
    std::vector<int> n_bonds_all(mol_end - mol_begin);

    for (unsigned int atom = mol_begin; atom != mol_end; ++atom) {
      for (Atom::bond_iterator bonded = top[atom].bondbegin(); bonded != top[atom].bondend(); ++bonded) {
        if (top[*bonded].Element() != Atom::HYDROGEN) {
          ++n_bonds_heavy[atom - mol_begin];
        }
        ++n_bonds_all[atom - mol_begin];
      }
    }
    int max_n_heavy = *std::max_element(n_bonds_heavy.begin(), n_bonds_heavy.end());
    bool use_all = max_n_heavy < 2;
    const std::vector<int>& n_bonds = use_all ? n_bonds_all : n_bonds_heavy;
    int i_most_bonded = std::max_element(n_bonds.begin(), n_bonds.end()) - n_bonds.begin();
    if (n_bonds[i_most_bonded] < 2) {
      mprinterr("Error: No atom with more than 1 bond found. This implementation needs at least 3 atoms in the main solvent.\n");
      return false;
    }
    rigidAtomIndices_[0] = i_most_bonded;
    const Atom& most_bonded = top[mol_begin + i_most_bonded];
    int i_output = 1;
    for (Atom::bond_iterator it = most_bonded.bondbegin(); i_output != 3; ++it) {
      if (use_all || top[*it].Element() != Atom::HYDROGEN) {
        rigidAtomIndices_[i_output] = *it - mol_begin;
        ++i_output;
      }
    }
  } else {
    // Set rigidAtomIndices_ from the specified names
    for (unsigned int i = 0; i < 3; ++i) {
      bool found_atom = false;
      for (unsigned int atom = mol_begin; atom != mol_end; ++atom) {
        if (top[atom].Name() == rigidAtomNames_[i]) {
          rigidAtomIndices_[i] = atom - mol_begin;
          found_atom = true;
          break;
        }
      }
      if (!found_atom) {
        mprinterr("Error: Solvent atom %s was not found.\n", rigidAtomNames_[i].c_str());
        return false;
      }
    }
  }
  mprintf("\tUsing atoms %s-%s-%s as rigid substructure of the solvent.\n",
          *top[mol_begin+rigidAtomIndices_[1]].Name(), *top[mol_begin+rigidAtomIndices_[0]].Name(), *top[mol_begin+rigidAtomIndices_[2]].Name());
  return true;
}

bool Action_GIST::createMoleculeDatasets()
{
  bool all_successful = true;
  for (unsigned int i = 0; i != solventNames_.size(); ++i) {
    std::string mol = solventNames_[i];
    molDensitySets_.push_back(AddDatasetAndFile("g_mol_" + mol, prefix_ + "-g-mol-" + mol + ext_, DataSet::GRID_FLT));
    if (!molDensitySets_.back()) {
      all_successful = false;
    }
  }
  return all_successful;
}

bool Action_GIST::createAtomDensityDatasets()
{
  int n_unique_elements = solventInfo_.unique_elements.size();
  bool all_successful = true;
  atomDensitySets_.assign(n_unique_elements, NULL);
  for (int i = 0; i < n_unique_elements; ++i) {
    std::string elem_name = solventInfo_.unique_elements[i];
    atomDensitySets_[i] = AddDatasetAndFile("g" + elem_name, prefix_ + "-g" + elem_name + ext_, DataSet::GRID_FLT);
    if (!atomDensitySets_[i]) {
      all_successful = false;
    }
  }
  return all_successful;
}

/**
 * \brief Set the solventType_ for each atom based on solventNames_;
 * 
 * For each atom i, solventType_[i] will be the index in solventNames_ 
 * that matches the atom's molecule name. Atoms that don't match any 
 * solvent will be UNKNOWN_MOLECULE.
 * 
 * \param top : current topology
 */
void Action_GIST::setSolventType(const Topology& top)
{
  size_t n_solvents = solventNames_.size();
  std::vector<CharMask> masks;
  masks.resize(n_solvents);
  for (size_t i = 0; i < n_solvents; ++i) {
    masks[i].SetMaskString(std::string(":") + solventNames_[i]);
    top.SetupCharMask(masks[i]);
  }
  solventType_.assign(top.Natom(), UNKNOWN_MOLECULE_);
  for (Topology::mol_iterator mol = top.MolStart(); mol != top.MolEnd(); ++mol) {
    int first = mol->MolUnit().Front();
    int last = mol->MolUnit().Back();
    for (size_t i = 0; i < n_solvents; ++i) {
      if (masks[i].AtomInCharMask(first)) {
        for (int atom = first; atom != last; ++atom) {
          solventType_[atom] = i;
        }
      }
    }
  }
}

void Action_GIST::setSoluteSolvent(const Topology& top)
{
  bool useMask = soluteMask_.length() > 0;
  CharMask isSolute_(soluteMask_);
  top.SetupCharMask(isSolute_);
  // If the mask is not used, we set the values via mol->IsSolvent.
  // Then, we copy the values to atomIsSolute and solvent_;
  if (!useMask) {
    for (Topology::mol_iterator mol=top.MolStart(); mol!=top.MolEnd(); ++mol) {
      bool is_solute = !mol->IsSolvent();
      for (int atom=mol->MolUnit().Front(); atom!=mol->MolUnit().Back(); ++atom) {
        atomIsSolute_[atom] = is_solute;
      }
    }
  } else {
    for (int i = 0; i != top.Natom(); ++i) {
      atomIsSolute_[i] = isSolute_.AtomInCharMask(i);
    }
  }
  for (int i = 0; i != top.Natom(); ++i) {
    if (atomIsSolute_[i]) {
      U_idxs_.push_back(i);
    }
    #ifdef CUDA
    this->solvent_[i] = !atomIsSolute_[i];
    #endif
  }

}

const Vec3 Action_GIST::x_lab_ = Vec3(1.0, 0.0, 0.0);
const Vec3 Action_GIST::y_lab_ = Vec3(0.0, 1.0, 0.0);
const Vec3 Action_GIST::z_lab_ = Vec3(0.0, 0.0, 1.0);
const double Action_GIST::QFAC_ = Constants::ELECTOAMBER * Constants::ELECTOAMBER;
const int Action_GIST::OFF_GRID_ = -1;
const int Action_GIST::UNKNOWN_MOLECULE_ = -1;

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
                                  neighborPerThread_);

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
  float* Neighbor = &(neighborPerThread_[0][0]);
  double Evdw, Eelec;
  int aidx;
  int maxAidx = frameIn.Natom();
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
  Neighbor = (&neighborPerThread_[mythread][0]);
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
    int a1 = aidx;            // Index of atom1
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
    if (!isMainSolvent(oidx1)) { continue; }
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
      if (isMainSolvent(oidx2) && oidx2 != oidx1)
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
    order_->UpdateVoxel(voxel1, (1.0 - (3.0/8)*sum));
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
  atom_voxel_.assign( frm.Frm().Natom(), OFF_GRID_ );

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

  Vec3 const& Origin = gridBin_->GridOrigin();
  // Loop over each solvent molecule
  for (Topology::mol_iterator mol = CurrentParm_->MolStart(); mol != CurrentParm_->MolEnd(); ++mol)
  {
    gist_grid_.Start();
    int mol_first = mol->MolUnit().Front();
    int mol_end = mol->MolUnit().Back();
    if (atomIsSolute_[mol_first]) { continue; }
    Vec3 mol_center = calcMolCenter(frm, mol_first, mol_end);
    // frm.Frm().VCenterOfMass(oidx, oidx+nMolAtoms_);
    Vec3 W_G = mol_center - Origin;
    gist_grid_.Stop();
    // Check if water oxygen is no more then 1.5 Ang from grid
    // NOTE: using <= to be consistent with original code
    if ( W_G[0] <= G_max_[0] && W_G[0] >= -1.5 &&
         W_G[1] <= G_max_[1] && W_G[1] >= -1.5 &&
         W_G[2] <= G_max_[2] && W_G[2] >= -1.5 )
    {
      // Try to bin the oxygen
      int voxel = calcVoxelIndex(mol_center[0], mol_center[1], mol_center[2]);
      if ( voxel != OFF_GRID_ )
      {
        // Oxygen is inside the grid. Record the voxel.
        const double* wXYZ = frm.Frm().XYZ( mol_first );
        for (int atom = mol_first; atom != mol_end; ++atom) {
          atom_voxel_[atom] = voxel;
          OnGrid_idxs_.push_back( atom );
          OnGrid_XYZ_.push_back( wXYZ[0] );
          OnGrid_XYZ_.push_back( wXYZ[1] );
          OnGrid_XYZ_.push_back( wXYZ[2] );
          wXYZ+=3;
          ++N_ON_GRID_;
        }
        ++N_solvent_[voxel];
        max_nwat_ = std::max( N_solvent_[voxel], max_nwat_ );

        if (solventType_[mol_first] != UNKNOWN_MOLECULE_) {
            molDensitySets_[solventType_[mol_first]]->UpdateVoxel(voxel, 1.0);
        }

        if (isMainSolvent(mol_first)) {
          // ----- EULER ---------------------------
          gist_euler_.Start();
          ++N_main_solvent_[voxel];
          // Record XYZ coords of water atoms (nonEP) in voxel TODO need EP?
          if (!skipS_) {
            voxel_xyz_[voxel].push_back( mol_center[0] );
            voxel_xyz_[voxel].push_back( mol_center[1] );
            voxel_xyz_[voxel].push_back( mol_center[2] );
            // Get O-HX vectors
            const double* O_XYZ  = frm.Frm().XYZ( mol_first + rigidAtomIndices_[0] );
            const double* H1_XYZ = frm.Frm().XYZ( mol_first + rigidAtomIndices_[1] );
            const double* H2_XYZ = frm.Frm().XYZ( mol_first + rigidAtomIndices_[2] );
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
          }
          //mprintf("DEBUG1: sidx= %u  voxel= %i  wxyz4= %g %g %g %g\n", sidx, voxel, w4, x4, y4, z4);
          //mprintf("DEBUG2: wxyz3= %g %g %g %g  wxyz2= %g %g %g %g  wxyz1= %g %g %g\n",
          //        w3, x3, y3, z3,
          //        w2, x2, y2, z2,
          //        w1, x1, y1, z1);
          // NOTE: No need for nw_angle_ here, it is same as N_solvent_
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
          for (int IDX = 0; IDX != mol_end-mol_first; IDX++) {
            const double* XYZ = frm.Frm().XYZ( mol_first+IDX );
            DPX += XYZ[0] * Q_.at(IDX);
            DPY += XYZ[1] * Q_.at(IDX);
            DPZ += XYZ[2] * Q_.at(IDX);
          }
          dipolex_->UpdateVoxel(voxel, DPX);
          dipoley_->UpdateVoxel(voxel, DPY);
          dipolez_->UpdateVoxel(voxel, DPZ);
          gist_dipole_.Stop();
        }
        // ---------------------------------------
      }

      // Water is at most 1.5A away from grid, so we need to check for H
      // even if O is outside grid.
      if (isMainSolvent(mol_first)) {
        for (int i = 0; i != mol_end-mol_first; ++i) {
          const double* xyz = frm.Frm().XYZ(mol_first + i);
          int vox = calcVoxelIndex(xyz[0], xyz[1], xyz[2]);
          if ( vox != OFF_GRID_ ) {
            int i_elem = solventInfo_.i_element[i];
            atomDensitySets_[i_elem]->UpdateVoxel(vox, 1.0);
          }
        }
      }
    } // END water is within 1.5 Ang of grid
  } // END loop over each solvent molecule

  // Do solute grid assignment for PME
  if (usePme_) {
    U_onGrid_idxs_.clear();
    gist_grid_.Start();
    for (unsigned int s = 0; s != U_idxs_.size(); s++)
    {
      int uidx = U_idxs_[s]; // the solute atom index
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
        int voxel = calcVoxelIndex(u_XYZ[0], u_XYZ[1], u_XYZ[2]);
        if ( voxel != OFF_GRID_ )
        {
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

bool Action_GIST::isMainSolvent(int atom) const
{
  if (solventNames_.size() == 0) {
    return !atomIsSolute_[atom];
  } else {
    return solventType_[atom] == 0;
  }
}

Vec3 Action_GIST::calcMolCenter(const ActionFrame& frm, int begin, int end) const
{
  if (useCom_) {
    return frm.Frm().VCenterOfMass(begin, end);
  } else {
    return Vec3(frm.Frm().XYZ(begin + rigidAtomIndices_[0]));
  }
}

int Action_GIST::calcVoxelIndex(double x, double y, double z) {
  size_t i, j, k;
  if ( gridBin_->Calc(x, y, z, i, j, k) ) {
    return int(i)*griddim_[1]*griddim_[2] + int(j)*griddim_[2] + int(k);
  }
  return OFF_GRID_;
}

// Action_GIST::SumEVV()
void Action_GIST::SumEVV() {
  for (unsigned int thread = 0; thread < E_VV_VDW_.size(); thread++) {
    for (unsigned int gr_pt = 0; gr_pt != MAX_GRID_PT_; gr_pt++) {
      Esw_->UpdateVoxel(gr_pt, E_UV_VDW_[thread][gr_pt] + E_UV_Elec_[thread][gr_pt]);
      Eww_->UpdateVoxel(gr_pt, (E_VV_VDW_[thread][gr_pt] + E_VV_Elec_[thread][gr_pt]) * 0.5);
      neighbor_->UpdateVoxel(gr_pt, neighborPerThread_[thread][gr_pt]);
    }
  }
}

Action_GIST::Farray Action_GIST::DataSetAsArray(const DataSet_3D& ds) const
{
  Farray ret(ds.Size());
  for (size_t i = 0; i < ds.Size(); ++i) {
    ret[i] = ds[i];
  }
  return ret;
}

void Action_GIST::ScaleDataSet(DataSet_3D& ds, double factor) const
{
  for (size_t i = 0; i < ds.Size(); ++i) {
    ds.SetGrid(i, ds[i]*factor);
  }
}

template<typename T>
std::vector<T> Action_GIST::NormalizeDataSet(const DataSet_3D& ds, const Action_GIST::Iarray& norm) const
{
  std::vector<T> ret(ds.Size());
  for (size_t i = 0; i < ds.Size(); ++i) {
    if (norm[i] == 0) {
      ret[i] = 0.0;
    } else {
      ret[i] = ds[i] / norm[i];
    }
  }
  return ret;
}

template<typename T>
std::vector<T> Action_GIST::WeightDataSet(const DataSet_3D& ds, double factor) const
{
  std::vector<T> ret(ds.Size());
  for (size_t i = 0; i < ds.Size(); ++i) {
    ret[i] = ds[i] * factor;
  }
  return ret;
}

template<typename T>
std::vector<T> Action_GIST::DensityWeightDataSet(const DataSet_3D& ds) const
{
  return WeightDataSet<T>(ds, 1.0 / (NFRAME_ * gridspacing_ * gridspacing_ * gridspacing_));
}

template<typename ARRAY_TYPE>
void Action_GIST::CopyArrayToDataSet(const ARRAY_TYPE& arr, DataSet_3D& ds) const
{
  for (size_t i = 0; i < arr.size(); ++i) {
    ds.SetGrid(i, arr[i]);
  }
}

double Action_GIST::SumDataSet(const DataSet_3D& ds) const
{
  double total = 0.0;
  for (size_t i = 0; i < ds.Size(); ++i) {
    total += ds[i];
  }
  return total;
}

/** Handle averaging for grids and output from GIST. */
void Action_GIST::Print() {
  if (!setupSuccessful_) {
    mprintf("Skipping GIST output because setup failed.\n");
    return;
  }
  gist_print_.Start();
  double Vvox = gridBin_->VoxelVolume();

  mprintf("    GIST OUTPUT:\n");

  // The variables are kept outside, so that they are declared for later use.
  // Calculate orientational entropy
  int nwtt = 0;
  if (! this->skipS_) {
    // LOOP over all voxels
    mprintf("\tCalculating orientational entropy:\n");
    ParallelProgress oe_progress( MAX_GRID_PT_ );
    int n_finished = 0;
#   ifdef _OPENMP
#   pragma omp parallel shared(n_finished) firstprivate(oe_progress)
    {
    oe_progress.SetThread( omp_get_thread_num() );
#   pragma omp for
#   endif
    for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
      oe_progress.Update( n_finished );
      int nw_total = N_main_solvent_[gr_pt]; // Total number of waters that have been in this voxel.
      //mprintf("DEBUG1: %u nw_total %i\n", gr_pt, nw_total);
      if (nw_total > 1) {
        double sorient_norm = 0.0;
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
            if (exactNnVolume_) {
              sorient_norm += log((NNr - sin(NNr)) * nw_total / Constants::PI);
            } else {
              sorient_norm += log(NNr*NNr*NNr*nw_total / (3.0*Constants::TWOPI));
            }
            //mprintf("DEBUG1: %u  nw_total= %i  NNr= %f  dbl= %f\n", gr_pt, nw_total, NNr, dbl);
          }
        } // END outer loop over all waters for this voxel
#       ifdef _OPENMP
#       pragma omp critical
#       endif
        {
          //mprintf("DEBUG1: dTSorient_norm %f\n", dTSorient_norm[gr_pt]);
          nwtt += nw_total;
          ++n_finished;
          dTSorient_->SetGrid(gr_pt, Constants::GASK_KCAL * temperature_ * nw_total
                                  * (sorient_norm/nw_total + Constants::EULER_MASC));
          //mprintf("DEBUG1: %f\n", dTSorienttot);
        }
      } else { // nw_total <= 1
#       ifdef _OPENMP
#       pragma omp atomic
#       endif
        ++n_finished;
      }
    } // END loop over all grid points (voxels)
#   ifdef _OPENMP
    }
#   endif
    oe_progress.Finish();
    infofile_->Printf("Maximum number of waters found in one voxel for %d frames = %d\n",
                      NFRAME_, max_nwat_);
    infofile_->Printf("Total referenced orientational entropy of the grid:"
                      " dTSorient = %9.5f kcal/mol, Nf=%d\n", SumDataSet(*dTSorient_) / NFRAME_, NFRAME_);
  }
  // Compute translational entropy for each voxel
  int nwts = 0;
  int nx = griddim_[0];
  int ny = griddim_[1];
  int nz = griddim_[2];

  Vec3 grid_origin(gridBin_->Center(0, 0, 0));

  // Loop over all grid points
  if (! this->skipS_)
    mprintf("\tCalculating translational entropy:\n");
  else
    mprintf("Calculating Densities:\n");
  for (unsigned int i_ds = 0; i_ds < solventInfo_.unique_elements.size(); ++i_ds) {
    ScaleDataSet(*atomDensitySets_[i_ds], 1.0 / (double(NFRAME_)*Vvox*BULK_DENS_*double(solventInfo_.element_count[i_ds])));
  }
  ParallelProgress te_progress( MAX_GRID_PT_ );
  int n_finished = 0;
# ifdef _OPENMP
# pragma omp parallel shared(n_finished) firstprivate(te_progress)
  {
  te_progress.SetThread( omp_get_thread_num() );
# pragma omp for
# endif
  for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
    te_progress.Update( n_finished );
    if (! this->skipS_) {
      int nw_total = N_main_solvent_[gr_pt]; // Total number of waters that have been in this voxel.
      int ix = gr_pt / (ny * nz);
      int iy = (gr_pt / nz) % ny;
      int iz = gr_pt % nz;
      bool boundary = ( ix == 0 || iy == 0 || iz == 0 || ix == (nx-1) || iy == (ny-1) || iz == (nz-1) );

      if ( !boundary ) {
        double strans_norm = 0.0;
        double ssix_norm = 0.0;
        int vox_nwts = 0;
        for (int n0 = 0; n0 < nw_total; ++n0)
        {
          Vec3 center(voxel_xyz_[gr_pt][3*n0], voxel_xyz_[gr_pt][3*n0+1], voxel_xyz_[gr_pt][3*n0+2]);
          int q0 = n0 * 4;  // index into voxel_Q_ for n0
          float W4 = voxel_Q_[gr_pt][q0  ];
          float X4 = voxel_Q_[gr_pt][q0+1];
          float Y4 = voxel_Q_[gr_pt][q0+2];
          float Z4 = voxel_Q_[gr_pt][q0+3];
          std::pair<double, double> NN = GistEntropyUtils::searchGridNearestNeighbors6D(
            center, W4, X4, Y4, Z4,
            voxel_xyz_, voxel_Q_,
            nx, ny, nz, grid_origin, gridspacing_,
            nNnSearchLayers_, n0);
          // It sometimes happens that we get numerically 0 values. 
          // Using a minimum distance changes the result only by a tiny amount 
          // (since those cases are rare), and avoids -inf values in the output.
          double NNd = std::max(sqrt(NN.first), GIST_TINY);
          double NNs = std::max(sqrt(NN.second), GIST_TINY);

          bool has_neighbor = NN.first < GistEntropyUtils::GIST_HUGE;
          if (has_neighbor) {
            ++vox_nwts;
            strans_norm += log((NNd*NNd*NNd*NFRAME_*4*Constants::PI*BULK_DENS_)/3);
            double sixDens = (NNs*NNs*NNs*NNs*NNs*NNs*NFRAME_*Constants::PI*BULK_DENS_) / 48;
            if (exactNnVolume_) {
              sixDens /= GistEntropyUtils::sixVolumeCorrFactor(NNs);
            }
            ssix_norm += log(sixDens);
            //mprintf("DEBUG1: dbl=%f NNs=%f\n", dbl, NNs);
          }
        } // END loop over all waters for this voxel
#       ifdef _OPENMP
#       pragma omp critical
#       endif
        {
          nwts += vox_nwts;
          ++n_finished;
          if (strans_norm != 0) {
            dTStrans_->SetGrid(gr_pt, Constants::GASK_KCAL * temperature_ * nw_total
                                    * (strans_norm/nw_total + Constants::EULER_MASC));
            dTSsix_->SetGrid(gr_pt, Constants::GASK_KCAL * temperature_ * nw_total
                                  * (ssix_norm/nw_total + Constants::EULER_MASC));
          }
        }
      } else { // boundary
#       ifdef _OPENMP
#       pragma omp atomic
#       endif
        ++n_finished;
      }
    }
  } // END loop over all grid points (voxels)
  #ifdef _OPENMP
  }
  #endif
  te_progress.Finish();
  if (!this->skipS_) {
    infofile_->Printf("watcount in vol = %d\n", nwtt);
    infofile_->Printf("watcount in subvol = %d\n", nwts);
    infofile_->Printf("Total referenced translational entropy of the grid:"
                      " dTStrans = %9.5f kcal/mol, Nf=%d\n", SumDataSet(*dTStrans_) / NFRAME_, NFRAME_);
    infofile_->Printf("Total 6d if all one vox: %9.5f kcal/mol\n", SumDataSet(*dTSsix_) / nwts);
    infofile_->Printf("Total t if all one vox: %9.5f kcal/mol\n", SumDataSet(*dTStrans_) / nwts);
    infofile_->Printf("Total o if all one vox: %9.5f kcal/mol\n", SumDataSet(*dTSorient_) / nwtt);
  }
  // free some memory before allocating all those Farrays for the -dens and -norm data.
  voxel_xyz_.clear();
  voxel_xyz_.shrink_to_fit();
  voxel_Q_.clear();
  voxel_Q_.shrink_to_fit();
  // Sum values from other threads if necessary
  #ifdef CUDA
  if (usePme_) { SumEVV(); }
  #else
  SumEVV();
  #endif

  // Remove solute-solvent energy in voxels without solvent (i.e., at the solute)
  for (size_t i = 0; i < Esw_->Size(); ++i) {
    if (N_solvent_[i] == 0) {
      Esw_->SetGrid(i, 0.0);
    }
  }

  for (unsigned int i = 0; i != solventNames_.size(); ++i) {
    ScaleDataSet(*molDensitySets_[i], 1.0 / (Vvox * BULK_DENS_ * NFRAME_));
  }
  // Compute average voxel energy. Allocate these sets even if skipping energy
  // to be consistent with previous output.

  Darray PME_norm, PME_dens, U_PME_dens;
  if (usePme_ && !skipE_) {
    CopyArrayToDataSet(E_pme_, *PME_);
    CopyArrayToDataSet(U_E_pme_, *U_PME_);
    mprintf("\t Calculating average voxel energies: \n");
    PME_norm = NormalizeDataSet<double>(*PME_, N_solvent_);
    PME_dens = DensityWeightDataSet<double>(*PME_);
    U_PME_dens = DensityWeightDataSet<double>(*U_PME_);
    infofile_->Printf("Ensemble total water energy on the grid: %9.5f Kcal/mol \n", SumDataSet(*PME_) / NFRAME_);
    infofile_->Printf("Ensemble total solute energy on the grid: %9.5f Kcal/mol \n", SumDataSet(*U_PME_) / NFRAME_);
  }
  if (!skipE_) {
    mprintf("\tCalculating average voxel energies:\n");
  }
  Farray Eww_norm = NormalizeDataSet<float>(*Eww_, N_solvent_);
  Farray Eww_dens = WeightDataSet<float>(*Eww_, 1.0 / (NFRAME_ * Vvox));
  Farray Esw_norm = NormalizeDataSet<float>(*Esw_, N_solvent_);
  Farray Esw_dens = DensityWeightDataSet<float>(*Esw_);
  if (!skipE_) {
    infofile_->Printf("Total water-solute energy of the grid: Esw = %9.5f kcal/mol\n", SumDataSet(*Esw_) / NFRAME_);
    infofile_->Printf("Total unreferenced water-water energy of the grid: Eww = %9.5f kcal/mol\n",
                      SumDataSet(*Eww_) / NFRAME_);
  }
  Farray neighbor_norm = NormalizeDataSet<float>(*neighbor_, N_main_solvent_);
  Farray neighbor_dens = DensityWeightDataSet<float>(*neighbor_);
  Farray dTStrans_norm = NormalizeDataSet<float>(*dTStrans_, N_main_solvent_);
  Farray dTSorient_norm = NormalizeDataSet<float>(*dTSorient_, N_main_solvent_);
  Farray dTSsix_norm = NormalizeDataSet<float>(*dTSsix_, N_main_solvent_);
  Farray dTStrans_dens = DensityWeightDataSet<float>(*dTStrans_);
  Farray dTSorient_dens = DensityWeightDataSet<float>(*dTSorient_);
  Farray dTSsix_dens = DensityWeightDataSet<float>(*dTSsix_);
  // Compute average dipole density.
  Farray dipolex_dens = WeightDataSet<float>(*dipolex_, 1.0 / (Constants::DEBYE_EA * NFRAME_ * Vvox));
  Farray dipoley_dens = WeightDataSet<float>(*dipoley_, 1.0 / (Constants::DEBYE_EA * NFRAME_ * Vvox));
  Farray dipolez_dens = WeightDataSet<float>(*dipolez_, 1.0 / (Constants::DEBYE_EA * NFRAME_ * Vvox));
  Darray order_norm = NormalizeDataSet<double>(*order_, N_main_solvent_);

  // Write final values to the datasets
  CopyArrayToDataSet(Esw_dens, *Esw_);
  CopyArrayToDataSet(Eww_dens, *Eww_);
  if (usePme_ && !skipE_) {
    CopyArrayToDataSet(PME_dens, *PME_);
    CopyArrayToDataSet(U_PME_dens, *U_PME_);
  }
  CopyArrayToDataSet(Eww_dens, *Eww_);
  CopyArrayToDataSet(dTStrans_dens, *dTStrans_);
  CopyArrayToDataSet(dTSorient_dens, *dTSorient_);
  CopyArrayToDataSet(dTSsix_dens, *dTSsix_);
  CopyArrayToDataSet(dipolex_dens, *dipolex_);
  CopyArrayToDataSet(dipoley_dens, *dipoley_);
  CopyArrayToDataSet(dipolez_dens, *dipolez_);
  CopyArrayToDataSet(order_norm, *order_);
  CopyArrayToDataSet(neighbor_norm, *neighbor_);

  for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++)
  {
    dipole_->SetGrid(gr_pt, sqrt( dipolex_dens[gr_pt] * dipolex_dens[gr_pt] +
                                  dipoley_dens[gr_pt] * dipoley_dens[gr_pt] +
                                  dipolez_dens[gr_pt] * dipolez_dens[gr_pt] ));
  }

  // Write the GIST output file.
  // TODO: Make a data file format?
  if (datafile_ != 0) {
    mprintf("\tWriting GIST results for each voxel:\n");
    const char* gistOutputVersion = "v4";
    // Do the header
    datafile_->Printf("GIST Output %s "
                      "spacing=%.4f center=%.6f,%.6f,%.6f dims=%i,%i,%i \n"
                      "voxel xcoord ycoord zcoord population",
                      gistOutputVersion, gridspacing_,
                      gridcntr_[0], gridcntr_[1], gridcntr_[2],
                      griddim_[0], griddim_[1], griddim_[2]);
    for (unsigned int i = 0; i < solventInfo_.unique_elements.size(); ++i) {
      datafile_->Printf( (std::string(" g_") + solventInfo_.unique_elements[i]).c_str() );
    }
    for (unsigned int i = 0; i != solventNames_.size(); ++i) {
      datafile_->Printf( (std::string(" g_mol_") + solventNames_[i]).c_str() );
    }
    datafile_->Printf(" dTStrans-dens(kcal/mol/A^3) dTStrans-norm(kcal/mol)"
                      " dTSorient-dens(kcal/mol/A^3) dTSorient-norm(kcal/mol)"
                      " dTSsix-dens(kcal/mol/A^3) dTSsix-norm(kcal/mol)"
                      " Esw-dens(kcal/mol/A^3) Esw-norm(kcal/mol)"
                      " Eww-dens(kcal/mol/A^3) Eww-norm-unref(kcal/mol)");
    if (usePme_)
      datafile_->Printf(" PME-dens(kcal/mol/A^3) PME-norm(kcal/mol)");
    datafile_->Printf(" Dipole_x-dens(D/A^3) Dipole_y-dens(D/A^3) Dipole_z-dens(D/A^3)"
                      " Dipole-dens(D/A^3) neighbor-dens(1/A^3) neighbor-norm order-norm\n");
    // Loop over voxels
    DataFilePrinter printer(*datafile_, fltFmt_, intFmt_);
    ProgressBar O_progress( MAX_GRID_PT_ );
    for (unsigned int gr_pt = 0; gr_pt < MAX_GRID_PT_; gr_pt++) {
      O_progress.Update( gr_pt );
      size_t i, j, k;
      Esw_->ReverseIndex( gr_pt, i, j, k );
      Vec3 XYZ = Esw_->Bin().Center( i, j, k );
      
      printer << gr_pt << XYZ[0] << XYZ[1] << XYZ[2] << N_solvent_[gr_pt];
      for (unsigned int i = 0; i < solventInfo_.unique_elements.size(); ++i) {
        printer << (*atomDensitySets_[i])[gr_pt];
      }
      for (unsigned int i = 0; i != solventNames_.size(); ++i) {
        printer << (*molDensitySets_[i])[gr_pt];
      }
      printer << dTStrans_dens[gr_pt] << dTStrans_norm[gr_pt] << dTSorient_dens[gr_pt] << dTSorient_norm[gr_pt] 
              << dTSsix_dens[gr_pt] << dTSsix_norm[gr_pt] << Esw_dens[gr_pt] << Esw_norm[gr_pt] << Eww_dens[gr_pt] << Eww_norm[gr_pt];
      if (usePme_) {
        printer << PME_dens[gr_pt] << PME_norm[gr_pt];
      }
      printer << dipolex_dens[gr_pt] << dipoley_dens[gr_pt] << dipolez_dens[gr_pt] 
              << (*dipole_)[gr_pt] << neighbor_dens[gr_pt] << neighbor_norm[gr_pt] << order_norm[gr_pt];
      printer.newline();
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
    int voxel = atom_voxel_[headAtomIndex];
    if (voxel != OFF_GRID_) {
      neighbor_->UpdateVoxel(voxel, result_n.at(headAtomIndex));
      for (unsigned int IDX = 0; IDX != nMolAtoms_; IDX++) {
        Esw_->UpdateVoxel(voxel, esw_result.at(headAtomIndex + IDX));
        Eww_->UpdateVoxel(voxel, eww_result.at(headAtomIndex + IDX));
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
        order_->UpdateVoxel(voxel, 1.0 - (3.0/8.0) * sum);
      }
    }

  }
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
