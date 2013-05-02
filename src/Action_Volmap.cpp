#include <cstdio>
#include <cmath>
#include "Action_Volmap.h"
#include "ArgList.h"
#include "Atom.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_Volmap::Action_Volmap() :
  //xcenter_(0.0),
  //ycenter_(0.0),
  //zcenter_(0.0),
  Nframes_(0),
  halfradii_(NULL),
  radscale_(1.0),
  xsize_(0.0),
  ysize_(0.0),
  zsize_(0.0)
{}

// Destructor
Action_Volmap::~Action_Volmap() {
  if (halfradii_ != NULL) delete[] halfradii_;
}

void Action_Volmap::Help() {
  RawHelp();
  mprintf("filename   -- Output file name\n");
  mprintf("dx, dy, dz -- grid spacing in the x-, y-, and z-dimensions, respectively.\n\n");
  mprintf("The grid size can be determined either by the size (x,y,z in Angstroms)\n");
  mprintf("or by a rectangular prism enclosing a mask with <buffer> clearance\n");
  mprintf("in every dimension. The density is calculated from the atoms in the\n");
  mprintf("required <mask>. If a <buffer> is given, the grid is centered on the\n");
  mprintf("centermask if provided, or the required mask if not.\n");
}

void Action_Volmap::RawHelp() {
  mprintf("\tfilename dx dy dz <mask> [xplor] [radscale <factor>]\n");
  mprintf("\t[ [[buffer <buffer>] [centermask <mask>]] | [center <x,y,z>] [size <x,y,z>] ]\n");
  mprintf("\t[peakcut <cutoff>] [peakfile <xyzfile>]\n");
}

// Action_Volmap::init()
Action::RetType Action_Volmap::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get the required mask
  std::string reqmask = actionArgs.GetMaskNext();
  if (reqmask.empty()) {
     mprinterr("Error: Volmap: no density mask specified.\n");
     return Action::ERR;
  }
  densitymask_.SetMaskString(reqmask);
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: Volmap: no filename specified.\n");
    return Action::ERR;
  }
  // Get grid resolutions
  dx_ = actionArgs.getNextDouble(0.0);
  dy_ = actionArgs.getNextDouble(0.0);
  dz_ = actionArgs.getNextDouble(0.0);

  // Get extra options
  peakcut_ = actionArgs.getKeyDouble("peakcut", 0.05);
  peakfilename_ = actionArgs.GetStringKey("peakfile");
  buffer_ = actionArgs.getKeyDouble("buffer", 3.0);
  radscale_ = 1.0 / actionArgs.getKeyDouble("radscale", 1.0);
  std::string sizestr = actionArgs.GetStringKey("size");
  std::string center = actionArgs.GetStringKey("centermask");
  std::string density = actionArgs.GetStringKey("density");
  dxform_ = !actionArgs.hasKey("xplor");

  // See how we are going to be setting up our grid
  if (!dsname.empty()) {
    // Get existing grid dataset
    grid_ = (DataSet_GridFlt*)DSL.FindSetOfType( dsname, DataSet::GRID_FLT );
    if (grid_ == 0) {
      mprinterr("Error: volmap: Could not find grid data set with name %s\n",
                dsname.c_str());
      return Action::ERR;
    }
  } else if (!sizestr.empty()) {
    // Now get our size from the specified arguments
    ArgList sizeargs = ArgList(sizestr, ",");
    xsize_ = sizeargs.getNextDouble(0.0);
    ysize_ = sizeargs.getNextDouble(0.0);
    zsize_ = sizeargs.getNextDouble(0.0);
    if (xsize_ <= 0 || ysize_ <= 0 || zsize_ <= 0) {
      mprinterr("Error: Volmap: Illegal grid sizes [%s]\n", sizestr.c_str());
      return Action::ERR;
    }
    //double size[3] = {xsize_, ysize_, zsize_};
    //double res[3] = {dx_, dy_, dz_};
    std::string centerstr = actionArgs.GetStringKey("center");
    ArgList centerargs = ArgList(centerstr, ",");
    xcenter_ = centerargs.getNextDouble(0.0);
    ycenter_ = centerargs.getNextDouble(0.0);
    zcenter_ = centerargs.getNextDouble(0.0);
    //grid_.GridInitSizeRes("volmap", size, res, std::string("origin"));
    // Determine bin count, ensure its even.
    int nx = (int) (xsize_ / dx_);
    int ny = (int) (ysize_ / dy_);
    int nz = (int) (zsize_ / dz_);
    GridAction::CheckEven(nx, "volmap");
    GridAction::CheckEven(ny, "volmap");
    GridAction::CheckEven(nz, "volmap");
    grid_ = GridAction::AllocateGrid(DSL, actionArgs.GetStringKey("name"),
                                     nx, ny, nz, 0.0, 0.0, 0.0, dx_, dy_, dz_);
    if (grid_ == 0) return Action::ERR;
  } else {
    // Now we generate our grid around our mask. See if we have a center mask
    if (center.empty())
      centermask_.SetMaskString(reqmask);
    else
      centermask_.SetMaskString(center);
    if (buffer_ < 0) {
      mprintf("Error: Volmap: The buffer must be non-negative.\n");
      return Action::ERR;
    }
  }
  // Info
  mprintf("\tGrid spacing will be %.2lfx%.2lfx%.2lf Angstroms\n", dx_, dy_, dz_);
  if (sizestr.empty())
    mprintf("\tGrid centered around %s with %.2lf Ang. clearance\n",
            centermask_.MaskExpression().c_str(), buffer_);
  else
    mprintf("\tGrid centered at origin.\n");
  mprintf("\tDensity will be printed in %s in ", filename_.c_str());
  if (dxform_)
    mprintf("OpenDX format\n");
  else
    mprintf("Xplor format\n");
  if (!peakfilename_.empty())
    mprintf("\tDensity peaks above %.3lf will be printed to %s in XYZ-format\n",
            peakcut_, peakfilename_.c_str());

  return Action::OK;
}

// Action_Volmap::setup()
Action::RetType Action_Volmap::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up our density mask and make sure it's not empty
  if (currentParm->SetupIntegerMask( densitymask_ ))
    return Action::ERR;
  if (densitymask_.None()) {
    mprinterr("Error: Volmap: Density mask selection empty!\n");
    return Action::ERR;
  }
  mprintf("Volmap: Grid mask [%s] selects %d atoms.\n", densitymask_.MaskString(),
          densitymask_.Nselected());
  // If we did not specify a size, make sure we have a valid centermask_
  if (xsize_ == 0) {
    if (currentParm->SetupIntegerMask( centermask_ ))
      return Action::ERR;
    // The masks must be populated
    if (centermask_.None()) {
      mprinterr("Error: Volmap: mask selection(s) empty!\n");
      return Action::ERR;
    }
    mprintf("Volmap: Centered mask [%s] selects %d atoms.\n", centermask_.MaskString(),
            centermask_.Nselected());
  }
  // Set up our radii_
  halfradii_ = new float[currentParm->Natom()];
  for (int i = 0; i < currentParm->Natom(); i++)
    halfradii_[i] = (float) GetRadius_(*currentParm, i) * radscale_ / 2;

  // DEBUG
//for (AtomMask::const_iterator it = densitymask_.begin(); it != densitymask_.end(); it++)
//  mprintf("Radius of atom %d is %f\n", *it, 2 * halfradii_[*it]);
  
  // Setup mask
  return Action::OK;
}

// Action_Volmap::GetRadius_()
/** Takes an input topology and gives back the VDW radius
  * \return vdW radius of the requested atom number from a Topology instance
  */
double Action_Volmap::GetRadius_(Topology const& top, int atom) {
  int param = (top.Ntypes() * (top[atom].TypeIndex()-1)) + top[atom].TypeIndex() - 1;
  int idx = top.NB_index()[param] - 1;
  return 0.5 * pow(2 * top.LJA()[idx] / top.LJB()[idx], 1.0/6.0);
}

// Action_Volmap::DoAction()
Action::RetType Action_Volmap::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // If this is our first frame, then we may have to set up our grid from our masks
  if (Nframes_ == 0) {
    if (xsize_ == 0) {
      double res[3] = {dx_, dy_, dz_};
      double xmin, xmax, ymin, ymax, zmin, zmax;
      AtomMask::const_iterator it = centermask_.begin();
      Vec3 cxyz = Vec3(currentFrame->XYZ(*it));
      xmin = xmax = cxyz[0];
      ymin = ymax = cxyz[1];
      zmin = zmax = cxyz[2];
      ++it;
      for (; it != centermask_.end(); it++) {
        Vec3 pt = Vec3(currentFrame->XYZ(*it));
        cxyz += pt;
        xmin = std::min(xmin, pt[0]);
        xmax = std::max(xmax, pt[0]);
        ymin = std::min(ymin, pt[1]);
        ymax = std::max(ymax, pt[1]);
        zmin = std::min(zmin, pt[2]);
        zmax = std::max(zmax, pt[2]);
      }
      center /= (double)centermask_.Nselected();
      xmin -= buffer_; 
      xmax += buffer_;
      ymin -= buffer_; 
      ymax += buffer_;
      zmin -= buffer_; 
      zmax += buffer_;
      int nx = (int) ((xmax - xmin) / dx_);
      int ny = (int) ((ymax - ymin) / dy_);
      int nz = (int) ((zmax - zmin) / dz_);
      //grid_.GridInitSizeRes("volmap", size, res, "origin");
      grid_ = AllocateGrid(*masterDSL_, dsname_, nx, ny, nz,
                           cxyz[0], cxyz[1], cxyz[2], dx_, dy_, dz_);
      if (grid_ == 0) return Action::ERR;
      // Get the center from the mask and assign the 'origin'
      //Vec3 center = currentFrame->VGeometricCenter(centermask_);
      //xcenter_ = center[0]; ycenter_ = center[1]; zcenter_ = center[2];
      xmin_ = xmin; 
      ymin_ = ymin; 
      zmin_ = zmin;
    }
  }
  // Now calculate the density for every point
  for (AtomMask::const_iterator atom = densitymask_.begin();
       atom != densitymask_.end(); atom++) {
    Vec3 pt = Vec3(currentFrame->XYZ(*atom));
    int ix = (int) ( floor( (pt[0]-xmin_) / dx_ + 0.5 ) );
    int iy = (int) ( floor( (pt[1]-ymin_) / dy_ + 0.5 ) );
    int iz = (int) ( floor( (pt[2]-zmin_) / dz_ + 0.5 ) );
    /* See how many steps in each dimension we smear out our Gaussian. This
     * formula is taken to be consistent with VMD's volmap tool
     */
    int nxstep = (int) ceil(4.1 * halfradii_[*atom] / dx_);
    int nystep = (int) ceil(4.1 * halfradii_[*atom] / dy_);
    int nzstep = (int) ceil(4.1 * halfradii_[*atom] / dz_);

    // Calculate the gaussian normalization factor (in 3 dimensions with the
    // given half-radius)
    double norm = 1 / (sqrt( 8.0 * PI*PI*PI ) * 
                       halfradii_[*atom]*halfradii_[*atom]*halfradii_[*atom]);
    double exfac = -1.0 / (2.0 * halfradii_[*atom] * halfradii_[*atom]);
    bool skip_point = (ix < -nxstep || ix > grid_.NX() + nxstep ||
                       iy < -nystep || iy > grid_.NY() + nystep ||
                       iz < -nzstep || iz > grid_.NZ() + nzstep);
    if (skip_point)
      continue;

    for (int xval = MAX(ix-nxstep,0); xval < MIN(ix+nxstep, grid_.NX()); xval++)
    for (int yval = MAX(iy-nystep,0); yval < MIN(iy+nystep, grid_.NY()); yval++)
    for (int zval = MAX(iz-nzstep,0); zval < MIN(iz+nzstep, grid_.NZ()); zval++) {
      Vec3 gridpt = Vec3(xmin_ + xval * dx_, ymin_ + yval * dy_, zmin_ + zval * dz_);
      double dist2 = DIST2_NoImage(pt, gridpt);
      grid_.AddDensity(xval, yval, zval, norm * exp(exfac * dist2));
    }
  }
  
  // Increment frame counter
  Nframes_++;
  return Action::OK;
}

inline size_t setStart(size_t xIn) {
  if (xIn == 0)
    return 0UL;
  else
    return xIn - 1L;
}

inline size_t setEnd(size_t xIn, size_t End) {
  size_t e = xIn + 2L;
  if (e > End)
    return End;
  else
    return e;
}

// Action_Volmap::print()
void Action_Volmap::Print() {

  // Divide our grid by the number of frames
  grid_ /= Nframes_;
  // Write the output file
  if (dxform_)
    grid_.PrintDX(filename_, xmin_, ymin_, zmin_);
  else
    grid_.PrintXplor( filename_, "This line is ignored", 
                      "rdparm generated grid density" );
  
  // See if we need to write the peaks out somewhere
  if (!peakfilename_.empty()) {
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
          size_t i_end = setEnd(i, grid_->NX());
          size_t j_end = setEnd(j, grid_->NY());
          size_t k_end = setEnd(k, grid_->NZ()); 
          for (size_t ii = setStart(i); ii < i_end; ii++)
            for (size_t jj = setStart(j); jj < j_end; jj++)
              for (size_t kk = setStart(k); kk < k_end; kk++) {
                if (ii==i && jj=j && kk=k) continue;
                if (grid_->GridVal(ii, jj, kk) > val)
                  peakgrid.setGrid(i,j,k,0.0f); // TODO: break after this?
              }
        }
    //peakgrid_ = grid_.ExtractPeaks(peakcut_);
    int npeaks = 0;
    std::vector<double> peakdata;
    for (int i = 0; i < peakgrid.NX(); i++)
      for (int j = 0; j < peakgrid.NY(); j++)
        for (int k = 0; k < peakgrid.NZ(); k++) {
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
      CpptrajFile outfile;
      if(outfile.OpenWrite(peakfilename_)) {
        mprinterr("Error: Could not open %s for writing.\n", peakfilename_.c_str());
        return;
      }
      outfile.Printf("%d\n\n", npeaks);
      for (int i = 0; i < npeaks; i++)
        outfile.Printf("C %16.8f %16.8f %16.8f %16.8f\n", peakdata[4*i],
                       peakdata[4*i+1], peakdata[4*i+2], peakdata[4*i+3]);
      outfile.CloseFile();
      mprintf("VolMap: %d density peaks found with higher density than %.4lf\n",
              npeaks, peakcut_);
    }else{
      mprintf("No peaks found with a density greater than %.3lf\n", peakcut_);
    }
  }

}
