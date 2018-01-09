#include <cmath> // pow, exp, sqrt
#include <algorithm> // std::min, std::max
#include "Action_Volmap.h"
#include "Constants.h" // PI
#include "CpptrajStdio.h"
#ifdef _OPENMP
# include <omp.h>
#endif

const double Action_Volmap::sqrt_8_pi_cubed = sqrt(8.0*Constants::PI*Constants::PI*Constants::PI);
// CONSTRUCTOR
Action_Volmap::Action_Volmap() :
  dx_(0.0), dy_(0.0), dz_(0.0),
  xmin_(0.0), ymin_(0.0), zmin_(0.0),
  Nframes_(0),
  setupGridOnMask_(false),
  grid_(0),
  peakfile_(0),
  peakcut_(0.05),
  buffer_(3.0),
  radscale_(1.0)
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
  mprintf("\t[out <filename>] <mask> [radscale <factor>]\n"
          "\t{ data <existing set> |\n"
          "\t  name <setname> <dx> <dy> <dz>\n"
          "\t    { size <x,y,z> [center <x,y,z>] |\n"
          "\t      centermask <mask> [buffer <buffer>] }\n"
          "\t}\n"
          "\t[peakcut <cutoff>] [peakfile <xyzfile>]\n");
}

// Action_Volmap::Init()
Action::RetType Action_Volmap::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  // Get specific keywords
  peakcut_ = actionArgs.getKeyDouble("peakcut", 0.05);
  peakfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("peakfile"), "Volmap Peaks");
  radscale_ = 1.0 / actionArgs.getKeyDouble("radscale", 1.0);
  // Determine how to set up grid: previous data set, 'size'/'center', or 'centermask'
  std::string dsname = actionArgs.GetStringKey("data");
  std::string setname, sizestr, centerstr;
  setupGridOnMask_ = false;
  if (dsname.empty()) {
    setname = actionArgs.GetStringKey("name");
    sizestr = actionArgs.GetStringKey("size");
    if (sizestr.empty()) {
      std::string cmask = actionArgs.GetStringKey("centermask");
      if (cmask.empty()) {
        mprinterr("Error: To set up grid must specify either 'data', 'size', or 'centermask'.\n");
        return Action::ERR;
      }
      // Center mask specified. Get buffer size.
      setupGridOnMask_ = true;
      centermask_.SetMaskString( cmask );
      buffer_ = actionArgs.getKeyDouble("buffer", 3.0);
      if (buffer_ < 0) {
        mprintf("Error: Volmap: The buffer must be non-negative.\n");
        return Action::ERR;
      }
    } else {
      // 'size' specified. Get center.
      centerstr = actionArgs.GetStringKey("center");
    }
    // Get grid resolutions
    dx_ = actionArgs.getNextDouble(0.0);
    dy_ = actionArgs.getNextDouble(0.0);
    dz_ = actionArgs.getNextDouble(0.0);
  }
  //std::string density = actionArgs.GetStringKey("density"); // FIXME obsolete?
  // Get the required mask
  std::string reqmask = actionArgs.GetMaskNext();
  if (reqmask.empty()) {
     mprinterr("Error: No atom mask specified.\n");
     return Action::ERR;
  }
  densitymask_.SetMaskString(reqmask);
  // Get output filename
  std::string outfilename = actionArgs.GetStringKey("out");
  if (outfilename.empty())
    outfilename = actionArgs.GetStringNext(); // Backwards compat.
  DataFile* outfile = init.DFL().AddDataFile( outfilename, actionArgs );

# ifdef _OPENMP
  // Always need to allocate temp grid space for each thread.
  int numthreads = 0;
# pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
      numthreads = omp_get_num_threads();
  }
  GRID_THREAD_.resize( numthreads );
# endif
  // Create new grid or use existing. 
  if (!dsname.empty()) {
    // Get existing grid dataset
    grid_ = (DataSet_GridFlt*)init.DSL().FindSetOfType( dsname, DataSet::GRID_FLT );
    if (grid_ == 0) {
      mprinterr("Error: Could not find grid data set with name '%s'\n",
                dsname.c_str());
      return Action::ERR;
    }
    if (!grid_->Bin().IsOrthoGrid()) {
      mprinterr("Error: Currently only works with orthogonal grids.\n");
      return Action::ERR;
    }
    setname = dsname;
    GridBin_Ortho const& gbo = static_cast<GridBin_Ortho const&>( grid_->Bin() );
    dx_ = gbo.DX();
    dy_ = gbo.DY();
    dz_ = gbo.DZ();
    Vec3 const& oxyz = gbo.GridOrigin();
    xmin_ = oxyz[0];
    ymin_ = oxyz[1];
    zmin_ = oxyz[2];
  } else {
    // Create new grid.
    grid_ = (DataSet_GridFlt*)init.DSL().AddSet(DataSet::GRID_FLT, setname, "VOLMAP");
    if (grid_ == 0) return Action::ERR;
    if (!sizestr.empty()) {
      // Get grid sizes from the specified arguments
      ArgList sizeargs = ArgList(sizestr, ",");
      double xsize = sizeargs.getNextDouble(0.0);
      double ysize = sizeargs.getNextDouble(0.0);
      double zsize = sizeargs.getNextDouble(0.0);
      if (xsize <= 0 || ysize <= 0 || zsize <= 0) {
        mprinterr("Error: Volmap: Illegal grid sizes [%s]\n", sizestr.c_str());
        return Action::ERR;
      }
      ArgList centerargs = ArgList(centerstr, ",");
      double xcenter = centerargs.getNextDouble(0.0);
      double ycenter = centerargs.getNextDouble(0.0);
      double zcenter = centerargs.getNextDouble(0.0);
      // Allocate grid
      if ( grid_->Allocate_X_C_D(Vec3(xsize,ysize,zsize), 
                                 Vec3(xcenter,ycenter,zcenter), 
                                 Vec3(dx_,dy_,dz_)) ) return Action::ERR;
      Vec3 const& oxyz = grid_->Bin().GridOrigin();
      xmin_ = oxyz[0];
      ymin_ = oxyz[1];
      zmin_ = oxyz[2];
    } 
  }
# ifdef _OPENMP
  // If not setting up grid via mask later, allocate temp space now.
  if (!setupGridOnMask_)
    for (Garray::iterator gt = GRID_THREAD_.begin(); gt != GRID_THREAD_.end(); ++gt)
      gt->resize( grid_->NX(), grid_->NY(), grid_->NZ() );
# endif
  // Setup output file
  if (outfile != 0) outfile->AddDataSet( grid_ );
  // Create total volume set
  total_volume_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(setname, "totalvol"));
  if (total_volume_ == 0) return Action::ERR;

  // Info
  mprintf("    VOLMAP: Grid spacing will be %.2fx%.2fx%.2f Angstroms\n", dx_, dy_, dz_);
  if (setupGridOnMask_)
    mprintf("\tGrid centered around mask %s with %.2f Ang. clearance\n",
            centermask_.MaskExpression().c_str(), buffer_);
  else
    grid_->GridInfo();
  mprintf("\tGridding atoms in mask '%s'\n", densitymask_.MaskString());
  mprintf("\tDividing radii by %f\n", 1.0/radscale_);
  if (outfile != 0)
    mprintf("\tDensity will wrtten to '%s'\n", outfile->DataFilename().full());
  mprintf("\tGrid dataset name is '%s'\n", grid_->legend());
  mprintf("\tTotal grid volume dataset name is '%s'\n", total_volume_->legend());
  if (peakfile_ != 0)
    mprintf("\tDensity peaks above %.3f will be printed to %s in XYZ-format\n",
            peakcut_, peakfile_->Filename().full());
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
  halfradii_.reserve( setup.Top().Natom() );
  if (setup.Top().Nonbond().HasNonbond()) {
    for (int i = 0; i < setup.Top().Natom(); i++)
      halfradii_.push_back( (float)(setup.Top().GetVDWradius(i) * radscale_ / 2) );
  } else {
    for (Topology::atom_iterator it = setup.Top().begin();
            it != setup.Top().end(); it++) {
      halfradii_.push_back( (float)(it->ElementRadius() * radscale_ / 2) );
    }
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
  int maxidx = densitymask_.Nselected();
  int midx, atom;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(midx, atom, mythread)
  {
  mythread = omp_get_thread_num();
# pragma omp for
# endif
  for (midx = 0; midx < maxidx; midx++) {
    atom = densitymask_[midx];
    Vec3 pt = Vec3(frm.Frm().XYZ(atom));
    int ix = (int) ( floor( (pt[0]-xmin_) / dx_ + 0.5 ) );
    int iy = (int) ( floor( (pt[1]-ymin_) / dy_ + 0.5 ) );
    int iz = (int) ( floor( (pt[2]-zmin_) / dz_ + 0.5 ) );
    /* See how many steps in each dimension we smear out our Gaussian. This
     * formula is taken to be consistent with VMD's volmap tool
     */
    int nxstep = (int) ceil(4.1 * halfradii_[atom] / dx_);
    int nystep = (int) ceil(4.1 * halfradii_[atom] / dy_);
    int nzstep = (int) ceil(4.1 * halfradii_[atom] / dz_);
    // Calculate the gaussian normalization factor (in 3 dimensions with the
    // given half-radius)
    double norm = 1 / (sqrt_8_pi_cubed * 
                       halfradii_[atom]*halfradii_[atom]*halfradii_[atom]);
    double exfac = -1.0 / (2.0 * halfradii_[atom] * halfradii_[atom]);
    if (ix < -nxstep || ix > nX + nxstep ||
        iy < -nystep || iy > nY + nystep ||
        iz < -nzstep || iz > nZ + nzstep)
      continue;

    int xend = std::min(ix+nxstep, nX);
    int yend = std::min(iy+nystep, nY);
    int zend = std::min(iz+nzstep, nZ);
    for (int xval = std::max(ix-nxstep, 0); xval < xend; xval++)
      for (int yval = std::max(iy-nystep, 0); yval < yend; yval++)
        for (int zval = std::max(iz-nzstep, 0); zval < zend; zval++) {
          Vec3 gridpt = Vec3(xmin_+xval*dx_, ymin_+yval*dy_, zmin_+zval*dz_) - pt;
          double dist2 = gridpt.Magnitude2();
#         ifdef _OPENMP
          GRID_THREAD_[mythread].incrementBy(xval, yval, zval, norm * exp(exfac * dist2));
#         else
          grid_->Increment(xval, yval, zval, norm * exp(exfac * dist2));
#         endif
        }
  } // END loop over atoms in densitymask_
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  
  // Increment frame counter
  Nframes_++;
  return Action::OK;
}

#ifdef MPI
int Action_Volmap::SyncAction() {
# ifdef _OPENMP
  CombineGridThreads();
# endif
  int total_frames = 0;
  trajComm_.Reduce( &total_frames, &Nframes_, 1, MPI_INT, MPI_SUM );
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
  float nf = (float)Nframes_;
  for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval)
    *gval /= nf;
  // Print volume estimate
  unsigned int nOccupiedVoxels = 0;
  for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval)
    if (*gval > 0.0) ++nOccupiedVoxels;
  double volume_estimate = (double)nOccupiedVoxels * grid_->Bin().VoxelVolume();
  total_volume_->Add(0, &volume_estimate);
  mprintf("\t%u occupied voxels, voxel volume= %f Ang^3, total volume %f Ang^3\n",
          nOccupiedVoxels, grid_->Bin().VoxelVolume(), volume_estimate);

//    grid_.PrintXplor( filename_, "This line is ignored", 
//                      "rdparm generated grid density" );
  
  // See if we need to write the peaks out somewhere
  if (peakfile_ != 0) {
    // Extract peaks from the current grid, setup another Grid instance. This
    // works by taking every grid point and analyzing all grid points adjacent
    // to it (including diagonals). If any of those grid points have a higher 
    // value (meaning there is a direction towards "increased" density) then 
    // that value is _not_ a maximum. Any density peaks less than the minimum
    // filter are discarded. The result is a Grid instance with all non-peak 
    // grid points zeroed-out.
    Grid<float> peakgrid = grid_->InternalGrid();
    for (size_t i = 0; i < grid_->NX(); i++)
      for (size_t j = 0; j < grid_->NY(); j++)
        for (size_t k = 0; k < grid_->NZ(); k++) {
          float val = grid_->GridVal(i, j, k);
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
    std::vector<double> peakdata;
    for (size_t i = 0; i < peakgrid.NX(); i++)
      for (size_t j = 0; j < peakgrid.NY(); j++)
        for (size_t k = 0; k < peakgrid.NZ(); k++) {
          double gval = peakgrid.element(i, j, k);
          if (gval > 0) {
            npeaks++;
            peakdata.push_back(xmin_+dx_*i);
            peakdata.push_back(ymin_+dy_*j);
            peakdata.push_back(zmin_+dz_*k);
            peakdata.push_back(gval);
          }
        }
    // If we have peaks, open up our peak data and print it
    if (npeaks > 0) {
      peakfile_->Printf("%d\n\n", npeaks);
      for (int i = 0; i < npeaks; i++)
        peakfile_->Printf("C %16.8f %16.8f %16.8f %16.8f\n", peakdata[4*i],
                       peakdata[4*i+1], peakdata[4*i+2], peakdata[4*i+3]);
      mprintf("Volmap: %d density peaks found with higher density than %.4lf\n",
              npeaks, peakcut_);
    }else{
      mprintf("No peaks found with a density greater than %.3lf\n", peakcut_);
    }
  }
}
