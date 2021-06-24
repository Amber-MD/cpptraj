#include <cmath> // pow, exp, sqrt
#include <algorithm> // std::min, std::max
#include "Action_Volmap.h"
#include "Constants.h" // PI
#include "CpptrajStdio.h"
#ifdef VOLMAP_DOUBLE
# include "DataSet_GridDbl.h"
# define VOLMAP_DS_T DataSet_GridDbl
# define VOLMAP_T double
#else
# include "DataSet_GridFlt.h"
# define VOLMAP_DS_T DataSet_GridFlt
# define VOLMAP_T float
#endif
#ifdef _OPENMP
# include <omp.h>
#endif
//#if defined(VOLMAP_USEFASTEXPS) || defined(VOLMAP_USEFASTEXP64)
//# incl ude "FastExp_Schraudolph.h"
//#elif defined(VOLMAP_USEFASTEXPIEEE)
//# incl ude "FastExp_fastexp.h"
//#endif

const double Action_Volmap::sqrt_8_pi_cubed_ = sqrt(8.0*Constants::PI*Constants::PI*Constants::PI);

// CONSTRUCTOR
Action_Volmap::Action_Volmap() :
  radiiType_(UNSPECIFIED),
  dx_(0.0),
  dy_(0.0),
  dz_(0.0),
  xmin_(0.0),
  ymin_(0.0),
  zmin_(0.0),
  debug_(0),
  Nframes_(0),
  setupGridOnMask_(false),
  spheremode_(false),
  grid_(0),
  total_volume_(0),
  calcpeaks_(false),
  peakfile_(0),
  peakdata_(0),
  peakcut_(0.05),
  buffer_(3.0),
  radscale_(1.0),
  stepfac_(4.1),
  splineDx_(0.01) // Recommendation from Roe & Brooks JMGM 2021
{}

void Action_Volmap::Help() const {
  RawHelp();
  mprintf("    filename  : Output file name\n"
          "    dx, dy, dz: grid spacing in the x-, y-, and z-dimensions, respectively.\n"
          "  The grid size can be determined either by the size (x,y,z in Angstroms)\n"
          "  or by a rectangular prism enclosing a mask with <buffer> clearance\n"
          "  in every dimension. The density is calculated from the atoms in the\n"
          "  required <mask>. If a <buffer> is given, the grid is centered on the\n"
          "  centermask if provided, or the required mask if not.\n");
}

void Action_Volmap::RawHelp() const {
  mprintf("\t[out <filename>] <mask> [radscale <factor>] [stepfac <fac>]\n"
          "\t[sphere] [radii {vdw | element}] [splinedx <spacing>]\n"
          "\t[calcpeaks] [peakcut <cutoff>] [peakfile <xyzfile>]\n"
          "\t{ data <existing set> |\n"
          "\t  name <setname> <dx> [<dy> <dz>]\n"
          "\t    { size <x,y,z> [center <x,y,z>] |\n"
          "\t      centermask <mask> [buffer <buffer>] |\n"
          "\t      boxref <reference> }\n"
          "\t}\n"
         );
}

// Action_Volmap::Init()
Action::RetType Action_Volmap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  // Get specific keywords
  peakcut_ = actionArgs.getKeyDouble("peakcut", 0.05);
  spheremode_ = actionArgs.hasKey("sphere");
  radscale_ = 1.0;
  stepfac_ = 4.1;
  if (spheremode_) {
    radscale_ = 0.5;
    stepfac_ = 1.0;
  }
  //splineDx_ = actionArgs.getKeyDouble("splinedx", 1.0/5000.0);
  splineDx_ = actionArgs.getKeyDouble("splinedx", 0.01);
  radscale_ = 1.0 / actionArgs.getKeyDouble("radscale", radscale_);
  stepfac_ = actionArgs.getKeyDouble("stepfac", stepfac_);
  std::string radarg = actionArgs.GetStringKey("radii");
  radiiType_ = UNSPECIFIED;
  if (!radarg.empty()) {
    if (radarg == "vdw")
      radiiType_ = VDW;
    else if (radarg == "element")
      radiiType_ = ELEMENT;
    else {
      mprinterr("Error: Unrecognized radii type: %s\n", radarg.c_str());
      return Action::ERR;
    }
  }
  // Determine how to set up grid: previous data set, 'size'/'center', or 'centermask'
  setupGridOnMask_ = false;
  enum SetupMode { DATASET=0, SIZE_CENTER, CENTERMASK, BOXREF, NMODE };
  SetupMode mode = NMODE;
  const char* SetupKey[5] = { "data", "size", "centermask", "boxref", 0 };
  std::string setup_arg;
  for (int i = 0; i < (int)NMODE; i++)
  {
    setup_arg = actionArgs.GetStringKey( SetupKey[i] );
    if (!setup_arg.empty()) {
      mode = (SetupMode)i;
      break;
    }
  }
  if (mode == NMODE) {
    mprinterr("Error: To set up grid must specify one of:");
    for (int i = 0; i < (int)NMODE; i++) mprinterr(" '%s'", SetupKey[i]);
    mprinterr("\n");
    return Action::ERR;
  }

# ifdef VOLMAP_DOUBLE
  static const DataSet::DataType GridDataType = DataSet::GRID_DBL;
# else
  static const DataSet::DataType GridDataType = DataSet::GRID_FLT;
# endif

  if (mode != DATASET) {
    Vec3 Sizes(0.0);
    Vec3 Center(0.0);
    std::string setname = actionArgs.GetStringKey("name");
    if (mode == SIZE_CENTER) {
      // 'size' specified. 
      ArgList sizeargs = ArgList(setup_arg, ",");
      Sizes[0] = sizeargs.getNextDouble(0.0);
      Sizes[1] = sizeargs.getNextDouble(0.0);
      Sizes[2] = sizeargs.getNextDouble(0.0);
      if (Sizes[0] <= 0 || Sizes[1] <= 0 || Sizes[2] <= 0) {
        mprinterr("Error: Volmap: Illegal grid sizes [%s]\n", setup_arg.c_str());
        return Action::ERR;
      }
      // Get 'center' arg.
      ArgList centerargs = ArgList(actionArgs.GetStringKey("center"), ",");
      Center[0] = centerargs.getNextDouble(0.0);
      Center[1] = centerargs.getNextDouble(0.0);
      Center[2] = centerargs.getNextDouble(0.0);
    } else if (mode == CENTERMASK) {
      // Center mask specified. Get buffer size.
      setupGridOnMask_ = true;
      if (centermask_.SetMaskString( setup_arg )) return Action::ERR;
      buffer_ = actionArgs.getKeyDouble("buffer", 3.0);
      if (buffer_ < 0) {
        mprinterr("Error: Volmap: The buffer must be non-negative.\n");
        return Action::ERR;
      }
    } else if (mode == BOXREF) {
      // Use reference box coords.
      DataSet_Coords_REF* REF = (DataSet_Coords_REF*)init.DSL().FindSetOfType( setup_arg, DataSet::REF_FRAME );
      if (REF == 0) {
        mprinterr("Error: Reference '%s' not found.\n", setup_arg.c_str());
        return Action::ERR;
      }
      if (REF->CoordsInfo().TrajBox().Is_X_Aligned_Ortho()) {
        mprinterr("Error: Reference '%s' does not have orthogonal box information.\n",
                  setup_arg.c_str());
        return Action::ERR;
      }
      Sizes = REF->CoordsInfo().TrajBox().Lengths();
      Center = REF->CoordsInfo().TrajBox().Center();
    }
    // Get grid resolutions
    dx_ = actionArgs.getNextDouble(0.0);
    dy_ = actionArgs.getNextDouble(dx_);
    dz_ = actionArgs.getNextDouble(dy_);
    if (dx_ <= 0.0 || dy_ <= 0.0 || dz_ <= 0.0) {
      mprinterr("Error: At least one of <dx> <dy> <dz> spacings must be specified\n"
                "Error:   if 'data' not specified.\n");
      return Action::ERR;
    }
    // Allocate grid dataset
    grid_ = (VOLMAP_DS_T*)init.DSL().AddSet(GridDataType, setname, "VOLMAP");
    if (grid_ == 0) {
      mprinterr("Error: Could not create grid dataset '%s'\n", setname.c_str());
      return Action::ERR;
    }
    if (!setupGridOnMask_) {
      // Allocate grid memory
      if ( grid_->Allocate_X_C_D(Sizes, Center, Vec3(dx_,dy_,dz_)) ) {
        mprinterr("Error: Could not allocate grid dataset '%s'\n", grid_->legend());
        Sizes.Print("Sizes");
        Center.Print("Center");
        return Action::ERR;
      }
      Vec3 const& oxyz = grid_->Bin().GridOrigin();
      xmin_ = oxyz[0];
      ymin_ = oxyz[1];
      zmin_ = oxyz[2];
    }
  } else {
    // Get existing grid dataset
    grid_ = (VOLMAP_DS_T*)init.DSL().FindSetOfType( setup_arg, GridDataType );
    if (grid_ == 0) {
      mprinterr("Error: Could not find grid data set with name '%s'\n",
                setup_arg.c_str());
      return Action::ERR;
    }
    if (!grid_->Bin().IsOrthoGrid()) {
      mprinterr("Error: Currently only works with orthogonal grids.\n");
      return Action::ERR;
    }
    dx_ = grid_->Bin().DX();
    dy_ = grid_->Bin().DY();
    dz_ = grid_->Bin().DZ();
    Vec3 const& oxyz = grid_->Bin().GridOrigin();
    xmin_ = oxyz[0];
    ymin_ = oxyz[1];
    zmin_ = oxyz[2];
  }
  // Warn for small grid spacing
  if (dx_ < 0.4 || dy_ < 0.4 || dz_ < 0.4)
    mprintf("Warning: Grid spacings smaller than 0.4 Ang. may be very slow;\n"
            "Warning:  consider using larger grid spacings.\n");
  //std::string density = actionArgs.GetStringKey("density"); // FIXME obsolete?
  // Get the required mask
  std::string reqmask = actionArgs.GetMaskNext();
  if (reqmask.empty()) {
     mprinterr("Error: No atom mask specified.\n");
     return Action::ERR;
  }
  if (densitymask_.SetMaskString(reqmask)) return Action::ERR;
  // See if peaks are being determined
  calcpeaks_ = actionArgs.hasKey("calcpeaks");
  peakfile_ = 0;
  std::string pfilename = actionArgs.GetStringKey("peakfile");
  if (!pfilename.empty()) {
    peakfile_ = init.DFL().AddDataFile(pfilename, actionArgs);
    if (peakfile_ == 0) {
      mprinterr("Error: Unable to allocate peak file.\n");
      return Action::ERR;
    }
    calcpeaks_ = true;
  }
  if (calcpeaks_) {
    peakdata_ = init.DSL().AddSet(DataSet::VECTOR_SCALAR, MetaData(grid_->Meta().Name(), "peaks"));
    if (peakdata_ == 0) {
      mprinterr("Error: Unable to allocate peak data set.\n");
      return Action::ERR;
    }
#   ifdef MPI
    peakdata_->SetNeedsSync( false );
#   endif
    if (peakfile_ != 0) peakfile_->AddDataSet( peakdata_ );
  }
  // Get output filename (newer syntax)
  std::string outfilename = actionArgs.GetStringKey("out");
  if (outfilename.empty())
    outfilename = actionArgs.GetStringNext(); // Backwards compat.
  DataFile* outfile = init.DFL().AddDataFile( outfilename, actionArgs );
  if (outfile != 0) outfile->AddDataSet( grid_ );

  // Create total volume set
  total_volume_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(grid_->Meta().Name(), "totalvol"));
  if (total_volume_ == 0) return Action::ERR;

# ifdef _OPENMP
  // Always need to allocate temp grid space for each thread.
  int numthreads = 0;
# pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
      numthreads = omp_get_num_threads();
  }
  GRID_THREAD_.resize( numthreads );
  // If not setting up grid via mask later, allocate temp space now.
  if (!setupGridOnMask_)
    for (Garray::iterator gt = GRID_THREAD_.begin(); gt != GRID_THREAD_.end(); ++gt)
      gt->resize( grid_->NX(), grid_->NY(), grid_->NZ() );
# endif

  // Info
  mprintf("    VOLMAP: Grid spacing will be %.2fx%.2fx%.2f Angstroms\n", dx_, dy_, dz_);
# ifdef VOLMAP_DOUBLE
  mprintf("\tUsing double precision grid.\n");
# else
  mprintf("\tUsing single precision grid.\n");
# endif
  if (setupGridOnMask_)
    mprintf("\tGrid centered around mask %s with %.2f Ang. clearance\n",
            centermask_.MaskExpression().c_str(), buffer_);
  else
    grid_->GridInfo();
  mprintf("\tGridding atoms in mask '%s'\n", densitymask_.MaskString());
  if (radiiType_ == VDW)
    mprintf("\tUsing van der Waals radii.\n");
  else if (radiiType_ == ELEMENT)
    mprintf("\tUsing elemental radii.\n");
  else if (radiiType_ == UNSPECIFIED)
    mprintf("\tUsing van der Waals radii if present, elemental radii otherwise.\n");
  if (spheremode_)
    mprintf("\tWhen smearing Gaussian, voxels farther than radii/2 will be skipped.\n");
  mprintf("\tDividing radii by %f\n", 1.0/radscale_);
  mprintf("\tFactor for determining number of bins to smear Gaussian is %f\n", stepfac_);
# if defined(VOLMAP_USEEXP)
  mprintf("\tUsing system exp() function for evaluating Gaussians.\n");
//# elif defined(VOLMAP_USEFASTEXPS)
//  mprintf("\tUsing exp() from N. Schraudolph, Neural Computation 11, 853–862 (1999).\n");
//# elif defined(VOLMAP_USEFASTEXP64)
//  mprintf("\tUsing 64 bit version of exp() from N. Schraudolph, Neural Computation 11, 853–862 (1999).\n");
//# elif defined(VOLMAP_USEFASTEXPIEEE)
//  mprintf("\tUsing exp() from Schraudolph & Malossi et al.\n");
# else /* VOLMAP_USEEXP */
  mprintf("\tExponential for Gaussians will be approximated using cubic splines\n"
          "\t  with a spacing of %g Ang.\n", splineDx_);
  mprintf("# Citation: Roe, D. R.; Brooks, B. R.; \"Improving the Speed of Volumetric\n"
          "#           Density Map Generation via Cubic Spline Interpolation\".\n"
          "#           Journal of Molecular Graphics and Modelling (2021).\n");
# if defined(VOLMAP_USEACCURATE)
  mprintf("\tSplines using more accurate but slower table lookup.\n");
# elif defined(VOLMAP_USEXTABLE)
  mprintf("\tSplines using less accurate lookup with tabled X values.\n");
# endif
# endif /* VOLMAP_USEEXP */
  if (outfile != 0)
    mprintf("\tDensity will wrtten to '%s'\n", outfile->DataFilename().full());
  mprintf("\tGrid dataset name is '%s'\n", grid_->legend());
  mprintf("\tTotal grid volume dataset name is '%s'\n", total_volume_->legend());
  if (calcpeaks_) {
    mprintf("\tDensity peaks above %.3f will be saved to %s\n", peakcut_, peakdata_->legend());
    if (peakfile_ != 0)
      mprintf("\tDensity peaks will be printed to %s in XYZ format.\n",
              peakfile_->DataFilename().full());
  }
# ifdef _OPENMP
  if (GRID_THREAD_.size() > 1)
    mprintf("\tParallelizing calculation with %zu threads.\n", GRID_THREAD_.size());
# endif

  return Action::OK;
}

// Action_Volmap::Setup()
Action::RetType Action_Volmap::Setup(ActionSetup& setup) {
  // Set up our density mask and make sure it's not empty
  if (setup.Top().SetupIntegerMask( densitymask_ ))
    return Action::ERR;
  if (densitymask_.None()) {
    mprinterr("Error: Volmap: Density mask selection empty!\n");
    return Action::ERR;
  }
  mprintf("\tVolmap: Grid mask [%s] selects %d atoms.\n", densitymask_.MaskString(),
          densitymask_.Nselected());
  // If we did not specify a size, make sure we have a valid centermask_
  if (setupGridOnMask_) {
    if (setup.Top().SetupIntegerMask( centermask_ ))
      return Action::ERR;
    // The masks must be populated
    if (centermask_.None()) {
      mprinterr("Error: Volmap: mask selection(s) empty!\n");
      return Action::ERR;
    }
    mprintf("\tVolmap: Centered mask [%s] selects %d atoms.\n", centermask_.MaskString(),
            centermask_.Nselected());
  }
  // Set up our radii_
  halfradii_.clear();
  halfradii_.reserve( densitymask_.Nselected() );
  Atoms_.clear();
  Atoms_.reserve( densitymask_.Nselected() );
  RadiiType radiiToUse = radiiType_;
  if (radiiToUse == VDW && !setup.Top().Nonbond().HasNonbond()) {
    mprinterr("Error: VDW radii specified but no VDW radii present in '%s'.\n",
              setup.Top().c_str());
    return Action::ERR;
  }
  if (radiiToUse == UNSPECIFIED) {
    if (setup.Top().Nonbond().HasNonbond())
      radiiToUse = VDW;
    else
      radiiToUse = ELEMENT;
  }
  double maxRad = 0;
  for (AtomMask::const_iterator atom = densitymask_.begin(); atom != densitymask_.end(); ++atom)
  {
    double rad = 0.0;
    if (radiiToUse == VDW)
      rad = setup.Top().GetVDWradius(*atom);
    else if (radiiToUse == ELEMENT)
      rad = setup.Top()[*atom].ElementRadius();
    if (rad > 0.0) {
      halfradii_.push_back( (double)(rad * radscale_ / 2) );
      Atoms_.push_back( *atom );
      // For determining function range.
      if (halfradii_.back() > maxRad) maxRad = halfradii_.back();
    }
  }
  // Try to determine the max value we may encounter for the exponential.
  //mprintf("\tMax observed radius: %g Ang\n", maxRad);
  int nxstep = (int) ceil(stepfac_ * maxRad / dx_)+1;
  int nystep = (int) ceil(stepfac_ * maxRad / dy_)+1;
  int nzstep = (int) ceil(stepfac_ * maxRad / dz_)+1;
  double maxx = (double)nxstep * dx_;
  double maxy = (double)nystep * dy_;
  double maxz = (double)nzstep * dz_;
  double maxDist = maxx*maxx + maxy*maxy + maxz*maxz;
  //mprintf("DEBUG: nx= %i  ny= %i  nz= %i\n", nxstep, nystep, nzstep);
  //mprintf("DEBUG: %g %g %g %g\n", maxx, maxy, maxz, maxDist);
  maxDist *= (-1.0 / (2.0 * maxRad * maxRad));
  //mprintf("DEBUG: max= %g\n", maxDist);
# ifndef VOLMAP_USEEXP
  // Set up the interpolation table
  if (table_.FillTable( exp, splineDx_, maxDist, 1.0 )) return Action::ERR;
  table_.PrintMemUsage("\t");
  if (debug_ > 0)
    table_.PrintTableInfo("DEBUG: ");
# endif

  if ((int)Atoms_.size() < densitymask_.Nselected())
    mprintf("Warning: %i atoms have 0.0 radii and will be skipped.\n",
            densitymask_.Nselected() - (int)Atoms_.size());
  if (Atoms_.empty()) {
    mprinterr("Error: No atoms with radii > 0.0\n");
    return Action::ERR;
  }
  // DEBUG
//for (AtomMask::const_iterator it = densitymask_.begin(); it != densitymask_.end(); it++)
//  mprintf("Radius of atom %d is %f\n", *it, 2 * halfradii_[*it]);
  
  return Action::OK;
}

// Action_Volmap::DoAction()
Action::RetType Action_Volmap::DoAction(int frameNum, ActionFrame& frm) {
  // If this is our first frame, then we may have to set up our grid from our masks
  if (Nframes_ == 0) {
    if (setupGridOnMask_) {
      // Determine min/max coord values for atoms in centermask. Calculate
      // geometric center while doing this.
      double xmin, xmax, ymin, ymax, zmin, zmax;
#     ifdef MPI
      if (trajComm_.Master()) {
#     endif
        AtomMask::const_iterator it = centermask_.begin();
        Vec3 cxyz = Vec3(frm.Frm().XYZ(*it));
        xmin = xmax = cxyz[0];
        ymin = ymax = cxyz[1];
        zmin = zmax = cxyz[2];
        ++it;
        for (; it != centermask_.end(); it++) {
          Vec3 pt = Vec3(frm.Frm().XYZ(*it));
          cxyz += pt;
          xmin = std::min(xmin, pt[0]);
          xmax = std::max(xmax, pt[0]);
          ymin = std::min(ymin, pt[1]);
          ymax = std::max(ymax, pt[1]);
          zmin = std::min(zmin, pt[2]);
          zmax = std::max(zmax, pt[2]);
        }
        cxyz /= (double)centermask_.Nselected();
        // Extend min/max by buffer.
        xmin -= buffer_; 
        xmax += buffer_;
        ymin -= buffer_; 
        ymax += buffer_;
        zmin -= buffer_; 
        zmax += buffer_;
#     ifdef MPI
      }
      // Send values to children //TODO put in array instead
      trajComm_.MasterBcast( &xmin, 1, MPI_DOUBLE );
      trajComm_.MasterBcast( &xmax, 1, MPI_DOUBLE );
      trajComm_.MasterBcast( &ymin, 1, MPI_DOUBLE );
      trajComm_.MasterBcast( &ymax, 1, MPI_DOUBLE );
      trajComm_.MasterBcast( &zmin, 1, MPI_DOUBLE );
      trajComm_.MasterBcast( &zmax, 1, MPI_DOUBLE );
#     endif
      // Allocate grid of given size centered on mask.
      if (grid_->Allocate_N_O_D( (xmax-xmin)/dx_, (ymax-ymin)/dy_, (zmax-zmin)/dz_,
                                 Vec3(xmin, ymin, zmin), Vec3(dx_, dy_, dz_) ))
        return Action::ERR;
#     ifdef _OPENMP
      for (Garray::iterator gt = GRID_THREAD_.begin(); gt != GRID_THREAD_.end(); ++gt)
        gt->resize( grid_->NX(), grid_->NY(), grid_->NZ() );
#     endif
      xmin_ = xmin;
      ymin_ = ymin;
      zmin_ = zmin;
      setupGridOnMask_ = false;
    }
  }
  // Now calculate the density for every point
  // TODO: Convert everything to size_t?
  int nX = (int)grid_->NX();
  int nY = (int)grid_->NY();
  int nZ = (int)grid_->NZ();
  int maxidx = (int)Atoms_.size();
  int midx, atom;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(midx, atom, mythread)
  {
  mythread = omp_get_thread_num();
# pragma omp for
# endif
  for (midx = 0; midx < maxidx; midx++) {
    double rhalf = (double)halfradii_[midx];
    double rcut2;
    if (spheremode_)
      rcut2 = rhalf*rhalf;
    else
      rcut2 = 99999999.0;
    atom = Atoms_[midx];
    Vec3 pt = Vec3(frm.Frm().XYZ(atom));
    int ix = (int) ( floor( (pt[0]-xmin_) / dx_ + 0.5 ) );
    /* See how many steps in X dimension we smear out our Gaussian. This
     * formula is taken to be consistent with VMD's volmap tool.
     */
    int nxstep = (int) ceil(stepfac_ * rhalf / dx_);
    if (ix >= -nxstep && ix <= nX + nxstep) {
      int iy = (int) ( floor( (pt[1]-ymin_) / dy_ + 0.5 ) );
      // See how many steps in Y dim we smear out Gaussian.
      int nystep = (int) ceil(stepfac_ * rhalf / dy_);
      if (iy >= -nystep && iy <= nY + nystep) {
        int iz = (int) ( floor( (pt[2]-zmin_) / dz_ + 0.5 ) );
        // See how many steps in Z dim we smear out Gaussian.
        int nzstep = (int) ceil(stepfac_ * rhalf / dz_);
        if (iz >= -nzstep && iz <= nZ + nzstep) {
          // Calculate the gaussian normalization factor (in 3 dimensions with the
          // given half-radius)
          double norm = 1 / (sqrt_8_pi_cubed_ * rhalf*rhalf*rhalf);
          double exfac = -1.0 / (2.0 * rhalf * rhalf);
          //mprintf("DBG: Atom %i norm %g exfac %g\n", atom+1, norm, exfac);

          int xend = std::min(ix+nxstep, nX);
          int yend = std::min(iy+nystep, nY);
          int zend = std::min(iz+nzstep, nZ);
          for (int xval = std::max(ix-nxstep, 0); xval < xend; xval++)
            for (int yval = std::max(iy-nystep, 0); yval < yend; yval++)
              for (int zval = std::max(iz-nzstep, 0); zval < zend; zval++) {
                Vec3 gridpt = Vec3(xmin_+xval*dx_, ymin_+yval*dy_, zmin_+zval*dz_) - pt;
                double dist2 = gridpt.Magnitude2();
                if (dist2 < rcut2) {
                  //mprintf("DEBUG: rhalf= %g  dist2= %g  exfac= %g  exp= %g\n", rhalf, dist2, exfac, exfac*dist2);
#                 ifdef _OPENMP
                  // NOTE: It is OK to call table_.Yval() here because in OpenMP
                  //       local variables of called functions are private.
#                 if defined(VOLMAP_USEEXP)
                  GRID_THREAD_[mythread].incrementBy(xval, yval, zval, norm * exp(exfac * dist2));
//#                 elif defined(VOLMAP_USEFASTEXPS)
//                  GRID_THREAD[mythread].incrementBy(xval, yval, zval, norm * FASTEXPS(exfac * dist2));
//#                 elif defined(VOLMAP_USEFASTEXP64)
//                  GRID_THREAD[mythread].incrementBy(xval, yval, zval, norm * fast_exps_64(exfac * dist2));
//#                 elif defined(VOLMAP_USEFASTEXPIEEE)
//                  GRID_THREAD[mythread].incrementBy(xval, yval, zval, norm * fastexp::exp<double, fastexp::IEEE, 4>(exfac * dist2));
#                 elif defined(VOLMAP_USEACCURATE)
                  GRID_THREAD_[mythread].incrementBy(xval, yval, zval, norm * table_.Yval_accurate(exfac * dist2));
#                 elif defined(VOLMAP_USEXTABLE)
                  GRID_THREAD_[mythread].incrementBy(xval, yval, xval, norm * table_.Yval_xtable(exfac * dist2));
#                 else
                  GRID_THREAD_[mythread].incrementBy(xval, yval, zval, norm * table_.Yval(exfac * dist2));
#                 endif
#                 else /* OPENMP */
#                 if defined(VOLMAP_USEEXP)
                  grid_->Increment(xval, yval, zval, norm * exp(exfac * dist2));
//#                 elif defined(VOLMAP_USEFASTEXPS)
//                  grid_->Increment(xval, yval, zval, norm * FASTEXPS(exfac * dist2));
//#                 elif defined(VOLMAP_USEFASTEXP64)
//                  grid_->Increment(xval, yval, zval, norm * fast_exps_64(exfac * dist2));
//#                 elif defined(VOLMAP_USEFASTEXPIEEE)
//                  grid_->Increment(xval, yval, zval, norm * fastexp::exp<double, fastexp::IEEE, 4>(exfac * dist2));
#                 elif defined(VOLMAP_USEACCURATE)
                  grid_->Increment(xval, yval, zval, norm * table_.Yval_accurate(exfac * dist2));
#                 elif defined(VOLMAP_USEXTABLE)
                  grid_->Increment(xval, yval, zval, norm * table_.Yval_xtable(exfac * dist2));
#                 else
                  grid_->Increment(xval, yval, zval, norm * table_.Yval(exfac * dist2));
#                 endif
#                 endif /* OPENMP */
                }
              } // END loop over zval
        } // END z dim is valid
      } // END y dim is valid
    } // END x dim is valid
  } // END loop over atoms in densitymask_
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  
  // Increment frame counter
  Nframes_++;
  return Action::OK;
}

#ifdef MPI
/** Sync grid to the master process. */
int Action_Volmap::SyncAction() {
# ifdef _OPENMP
  CombineGridThreads();
# endif
  int total_frames = 0;
  trajComm_.ReduceMaster( &total_frames, &Nframes_, 1, MPI_INT, MPI_SUM );
  if (trajComm_.Master())
    Nframes_ = total_frames;
  return 0;
}
#endif

#ifdef _OPENMP
/** Combine results from each GRID_THREAD into main grid */
void Action_Volmap::CombineGridThreads() {
  if (!GRID_THREAD_.empty()) {
    for (Garray::const_iterator gt = GRID_THREAD_.begin(); gt != GRID_THREAD_.end(); ++gt)
      for (unsigned int idx = 0; idx != gt->size(); idx++)
        (*grid_)[idx] += (*gt)[idx];
    GRID_THREAD_.clear();
  }
}
#endif

// Need this instead of MAX since size_t can never be negative
inline size_t setStart(size_t xIn) {
  if (xIn == 0)
    return 0UL;
  else
    return xIn - 1L;
}

// Action_Volmap::Print()
void Action_Volmap::Print() {
  if (Nframes_ < 1) return;
# ifdef _OPENMP
  CombineGridThreads();
# endif
  // Divide our grid by the number of frames
  VOLMAP_T nf = (VOLMAP_T)Nframes_;
  for (VOLMAP_DS_T::iterator gval = grid_->begin(); gval != grid_->end(); ++gval)
    *gval /= nf;
  // Print volume estimate
  unsigned int nOccupiedVoxels = 0;
  for (VOLMAP_DS_T::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
    if (*gval > 0.0) {
      ++nOccupiedVoxels;
      //mprintf("DBG: %16.8e\n", *gval);
    }
  }
  double volume_estimate = (double)nOccupiedVoxels * grid_->Bin().VoxelVolume();
  total_volume_->Add(0, &volume_estimate);
  mprintf("\t%u occupied voxels, voxel volume= %f Ang^3, total volume %f Ang^3\n",
          nOccupiedVoxels, grid_->Bin().VoxelVolume(), volume_estimate);

//    grid_.PrintXplor( filename_, "This line is ignored", 
//                      "rdparm generated grid density" );
  
  // See if we need to write the peaks out somewhere
  if (calcpeaks_) {
    // Extract peaks from the current grid, setup another Grid instance. This
    // works by taking every grid point and analyzing all grid points adjacent
    // to it (including diagonals). If any of those grid points have a higher 
    // value (meaning there is a direction towards "increased" density) then 
    // that value is _not_ a maximum. Any density peaks less than the minimum
    // filter are discarded. The result is a Grid instance with all non-peak 
    // grid points zeroed-out.
    Grid<VOLMAP_T> peakgrid = grid_->InternalGrid();
    for (size_t i = 0; i < grid_->NX(); i++)
      for (size_t j = 0; j < grid_->NY(); j++)
        for (size_t k = 0; k < grid_->NZ(); k++) {
          VOLMAP_T val = grid_->GridVal(i, j, k);
          if (val < peakcut_) {
            peakgrid.setGrid(i, j, k, 0.0f);
            continue;
          }
          size_t i_end = std::min(i+2, grid_->NX());
          size_t j_end = std::min(j+2, grid_->NY());
          size_t k_end = std::min(k+2, grid_->NZ()); 
          for (size_t ii = setStart(i); ii < i_end; ii++)
            for (size_t jj = setStart(j); jj < j_end; jj++)
              for (size_t kk = setStart(k); kk < k_end; kk++) {
                if (ii==i && jj==j && kk==k) continue;
                if (grid_->GridVal(ii, jj, kk) > val)
                  peakgrid.setGrid(i,j,k,0.0f); // TODO: break after this?
              }
        }
    int npeaks = 0;
    //std::vector<double> peakdata;
    double pbuf[4];
    for (size_t i = 0; i < peakgrid.NX(); i++)
      for (size_t j = 0; j < peakgrid.NY(); j++)
        for (size_t k = 0; k < peakgrid.NZ(); k++) {
          double gval = peakgrid.element(i, j, k);
          if (gval > 0) {
            //npeaks++;
            //peakdata.push_back(xmin_+dx_*i);
            //peakdata.push_back(ymin_+dy_*j);
            //peakdata.push_back(zmin_+dz_*k);
            //peakdata.push_back(gval);
            pbuf[0] = xmin_+dx_*(double)i;
            pbuf[1] = ymin_+dy_*(double)j;
            pbuf[2] = zmin_+dz_*(double)k;
            pbuf[3] = gval;
            peakdata_->Add(npeaks++, pbuf);
          }
        }
    // Report on how many peaks found. 
    if (npeaks > 0) {
      mprintf("Volmap: %d density peaks found with higher density than %.4lf\n",
              npeaks, peakcut_);
    } else {
      mprintf("No peaks found with a density greater than %.3lf\n", peakcut_);
    }
  }
}
#undef VOLMAP_DS_T
#undef VOLMAP_T
