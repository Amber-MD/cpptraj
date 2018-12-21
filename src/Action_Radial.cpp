// Action_Radial
#include <cmath> // sqrt
#include "Action_Radial.h"
#include "CpptrajStdio.h"
#include "Constants.h" // FOURTHIRDSPI
#ifdef _OPENMP
#  include <omp.h>
#endif

// CONSTRUCTOR
Action_Radial::Action_Radial() :
  RDF_(0),
  rdf_thread_(0),
  rmode_(NORMAL),
  currentParm_(0),
  intramol_distances_(0),
  useVolume_(false),
  volume_(0),
  maximum2_(0),
  spacing_(-1),
  one_over_spacing_(-1),
  numBins_(0),
  numthreads_(1),
  numFrames_(0),
  // Default particle density (molecules/Ang^3) for water based on 1.0 g/mL
  density_(0.033456),
  Dset_(0),
  intrdf_(0),
  rawrdf_(0),
  debug_(0)
{} 

void Action_Radial::Help() const {
  mprintf("\t[out <outfilename>] <spacing> <maximum> <solvent mask1> [<solute mask2>] [noimage]\n"
          "\t[density <density> | volume] [<dataset name>] [intrdf <file>] [rawrdf <file>]\n"
          "\t[{{center1|center2|nointramol} | [byres1] [byres2] [bymol1] [bymol2]}]\n"
          "  Calculate the radial distribution function (RDF) of atoms in <solvent mask1>.\n"
          "  If <solute mask2> is given calculate RDF of all atoms in <solvent mask1>\n"
          "  to each atom in <solute mask2>.\n"
          "  center1|center2 will use the center of *all* atoms selected by masks 1 and 2 respectively.\n"
          "  nointramol will ignore distances when both atoms are part of the same molecule.\n"
          "  If byresX or bymolX are specified, distances will be between the centers of mass\n"
          "  of residues/molecules selected by mask1 or mask2.\n");
}

// DESTRUCTOR
Action_Radial::~Action_Radial() {
  //fprintf(stderr,"Radial Destructor.\n");
  if (RDF_!=0) delete[] RDF_;
  if (rdf_thread_!=0) {
    for (int i=0; i < numthreads_; i++)
      delete[] rdf_thread_[i];
    delete[] rdf_thread_;
  }
}

inline Action::RetType RDF_ERR(const char* msg) {
  mprinterr("Error: %s\n", msg);
  return Action::ERR;
}

// Action_Radial::Init()
Action::RetType Action_Radial::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  std::string outfilename = actionArgs.GetStringKey("out");
  // Default particle density (mols/Ang^3) for water based on 1.0 g/mL
  density_ = actionArgs.getKeyDouble("density",0.033456);
  // Determine mode, by site TODO better integrate with other modes
  siteMode1_ = OFF;
  siteMode2_ = OFF;
  if (actionArgs.hasKey("byres1")) siteMode1_ = BYRES;
  if (actionArgs.hasKey("bymol1")) siteMode1_ = BYMOL;
  if (actionArgs.hasKey("byres2")) siteMode2_ = BYRES; 
  if (actionArgs.hasKey("bymol2")) siteMode2_ = BYMOL; 
  // Determine mode, other
  if (actionArgs.hasKey("center1"))
    rmode_ = CENTER1;
  else if (actionArgs.hasKey("center2"))
    rmode_ = CENTER2;
  else if (actionArgs.hasKey("nointramol"))
    rmode_ = NO_INTRAMOL;
  else
    rmode_ = NORMAL;
  // Check for mode incompatibility
  if (siteMode1_ != OFF || siteMode2_ != OFF) {
    if (rmode_ != NORMAL) {
      mprinterr("Error: 'byres'/'bymol' mode cannot be active with other modes (center, nointramol).\n");
      return Action::ERR;
    }
    rmode_ = BYSITE;
  }
  useVolume_ = actionArgs.hasKey("volume");
  DataFile* intrdfFile = init.DFL().AddDataFile(actionArgs.GetStringKey("intrdf"));
  DataFile* rawrdfFile = init.DFL().AddDataFile(actionArgs.GetStringKey("rawrdf"));
  spacing_ = actionArgs.getNextDouble(-1.0);
  if (spacing_ < 0) {
    mprinterr("Error: Radial: No spacing argument or arg < 0.\n");
    Help();
    return Action::ERR;
  }
  double maximum = actionArgs.getNextDouble(-1.0);
  if (maximum < 0) {
    mprinterr("Error: Radial: No maximum argument or arg < 0.\n");
    Help();
    return Action::ERR;
  }
  // Store max^2, distances^2 greater than max^2 do not need to be
  // binned and therefore do not need a sqrt calc.
  maximum2_ = maximum * maximum;

  // Get First Mask
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: Radial: No mask given.\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);

  // Check for second mask - if none specified use first mask
  std::string mask2 = actionArgs.GetMaskNext();
  if (!mask2.empty()) 
    Mask2_.SetMaskString(mask2);
  else
    Mask2_.SetMaskString(mask1);
  // If filename not yet specified check for backwards compat.
  if (outfilename.empty() && actionArgs.Nargs() > 1 && !actionArgs.Marked(1))
    outfilename = actionArgs.GetStringNext();

  // Set up output dataset.
  Dset_ = init.DSL().AddSet( DataSet::DOUBLE, MetaData(actionArgs.GetStringNext(), "",
                                                       MetaData::NOT_TS), "g(r)");
  if (Dset_ == 0) return RDF_ERR("Could not allocate RDF data set.");
  DataFile* outfile = init.DFL().AddDataFile(outfilename, actionArgs);
  if (outfile != 0) outfile->AddDataSet( Dset_ );
  // Make default precision a little higher than normal
  Dset_->SetupFormat().SetFormatWidthPrecision(12,6);
  // Set DataSet legend from mask strings.
  Dset_->SetLegend(Mask1_.MaskExpression() + " => " + Mask2_.MaskExpression());
  // TODO: Set Yaxis label in DataFile
  // Calculate number of bins
  one_over_spacing_ = 1 / spacing_;
  double temp_numbins = maximum * one_over_spacing_;
  temp_numbins = ceil(temp_numbins);
  numBins_ = (int) temp_numbins;
  // Setup output datafile. Align on bin centers instead of left.
  // TODO: Use Rdim for binning?
  Dimension Rdim( spacing_ / 2.0, spacing_, "Distance (Ang)" ); 
  Dset_->SetDim(Dimension::X, Rdim);
  // Set up output for integral of mask2 if specified.
  if (intrdfFile != 0) {
    intrdf_ = init.DSL().AddSet( DataSet::DOUBLE, MetaData(Dset_->Meta().Name(), "int",
                                                           MetaData::NOT_TS) );
    if (intrdf_ == 0) return RDF_ERR("Could not allocate RDF integral data set.");
    intrdf_->SetupFormat().SetFormatWidthPrecision(12,6);
    intrdf_->SetLegend("Int[" + Mask2_.MaskExpression() + "]");
    intrdf_->SetDim(Dimension::X, Rdim);
    intrdfFile->AddDataSet( intrdf_ );
  } else
    intrdf_ = 0;
  // Set up output for raw rdf
  if (rawrdfFile != 0) {
    rawrdf_ = init.DSL().AddSet( DataSet::DOUBLE, MetaData(Dset_->Meta().Name(), "raw",
                                                           MetaData::NOT_TS) );
    if (rawrdf_ == 0) return RDF_ERR("Could not allocate raw RDF data set.");
    rawrdf_->SetupFormat().SetFormatWidthPrecision(12,6);
    rawrdf_->SetLegend("Raw[" + Dset_->Meta().Legend() + "]");
    rawrdf_->SetDim(Dimension::X, Rdim);
    rawrdfFile->AddDataSet( rawrdf_ );
  } else
    rawrdf_ = 0;
# ifdef MPI
  // These do not need to be synced since they are not time series
  Dset_->SetNeedsSync( false );
  if (intrdf_ != 0) intrdf_->SetNeedsSync( false );
  if (rawrdf_ != 0) rawrdf_->SetNeedsSync( false );
# endif
  // Set up histogram
  RDF_ = new int[ numBins_ ];
  std::fill(RDF_, RDF_ + numBins_, 0);
# ifdef _OPENMP
  // Since RDF is shared by all threads and we cant guarantee that a given
  // bin in RDF wont be accessed at the same time by the same thread,
  // each thread needs its own bin space.
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    numthreads_ = omp_get_num_threads();
}
  rdf_thread_ = new int*[ numthreads_ ];
  for (int i=0; i < numthreads_; i++) {
    rdf_thread_[i] = new int[ numBins_ ];
    std::fill(rdf_thread_[i], rdf_thread_[i] + numBins_, 0);
  }
# endif
  
  mprintf("    RADIAL: Calculating RDF for atoms in mask [%s]",Mask1_.MaskString());
  if (!mask2.empty()) 
    mprintf(" to atoms in mask [%s]",Mask2_.MaskString());
  mprintf("\n");
  if (outfile != 0)
    mprintf("\tOutput to %s.\n", outfile->DataFilename().full());
  if (intrdf_ != 0)
    mprintf("\tIntegral of mask2 atoms will be output to %s\n",
            intrdfFile->DataFilename().full());
  if (rawrdf_ != 0)
    mprintf("\tRaw RDF bin values will be output to %s\n",
            rawrdfFile->DataFilename().full());
  if (rmode_ == BYSITE) {
    if (siteMode1_ == BYRES)
      mprintf("\tUsing center of residues selected by mask1 '%s'\n", Mask1_.MaskString());
    else if (siteMode1_ == BYMOL)
      mprintf("\tUsing center of molecules selected by mask1 '%s'\n", Mask1_.MaskString());
    if (siteMode2_ == BYRES)
      mprintf("\tUsing center of residues selected by mask2 '%s'\n", Mask2_.MaskString());
    else if (siteMode2_ == BYMOL)
      mprintf("\tUsing center of molecules selected by mask2 '%s'\n", Mask2_.MaskString());
  } else {
    if (rmode_==CENTER1)
      mprintf("\tUsing center of all atoms selected by mask1.\n");
    else if (rmode_==CENTER2)
      mprintf("\tUsing center of all atoms selected by mask2.\n");
    else if (rmode_==NO_INTRAMOL)
      mprintf("\tIgnoring intramolecular distances.\n");
  }
  mprintf("\tHistogram max %f, spacing %f, bins %i.\n",maximum,
          spacing_,numBins_);
  if (useVolume_)
    mprintf("\tNormalizing based on cell volume.\n");
  else
    mprintf("\tNormalizing using particle density of %f molecules/Ang^3.\n",density_);
  if (!image_.UseImage()) 
    mprintf("\tImaging disabled.\n");
  if (numthreads_ > 1)
    mprintf("\tParallelizing RDF calculation with %i threads.\n",numthreads_);

  return Action::OK;
}

/** Set up site array by atom. */
int Action_Radial::SetupSiteArrayByAtom(Marray& sites, AtomMask const& mask)
const
{
  sites.clear();
  sites.reserve(mask.Nselected());
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
    sites.push_back( AtomMask(*at) );
  return 0;
}

/** Set up site array by residue. */
int Action_Radial::SetupSiteArrayByRes(Marray& sites, Topology const& top, AtomMask const& mask)
const
{
  if (mask.Nselected() < 1) return 1;
  sites.clear();
  int lastRes = top[ mask[0] ].ResNum();
  sites.push_back( AtomMask() );
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
  {
    int currentRes = top[ *at ].ResNum();
    if (currentRes != lastRes) {
      sites.push_back( AtomMask() );
      lastRes = currentRes;
    }
    sites.back().AddSelectedAtom( *at );
  }
  // DEBUG
  if (debug_ > 1) {
    mprintf("DEBUG: Sites selected by residue for '%s'\n", mask.MaskString());
    for (Marray::const_iterator m = sites.begin(); m != sites.end(); ++m) {
      mprintf("%8u :", m - sites.begin());
      for (AtomMask::const_iterator at = m->begin(); at != m->end(); at++)
        mprintf(" %i", *at);
      mprintf("\n");
    }
  }
  return 0;
}

/** Set up site array by molecule. */
int Action_Radial::SetupSiteArrayByMol(Marray& sites, Topology const& top, AtomMask const& mask)
const
{
  if (mask.Nselected() < 1) return 1;
  if (top.Nmol() < 1) {
    mprinterr("Error: No topology info for '%s', cannot set up sites by molecule.\n");
    return -1;
  }
  sites.clear();
  int lastMol = top[ mask[0] ].MolNum();
  sites.push_back( AtomMask() );
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
  {
    int currentMol = top[ *at ].MolNum();
    if (currentMol != lastMol) {
      sites.push_back( AtomMask() );
      lastMol = currentMol;
    }
    sites.back().AddSelectedAtom( *at );
  }
  // DEBUG
  if (debug_ > 1) {
    mprintf("DEBUG: Sites selected by molecule for '%s'\n", mask.MaskString());
    for (Marray::const_iterator m = sites.begin(); m != sites.end(); ++m) {
      mprintf("%8u :", m - sites.begin());
      for (AtomMask::const_iterator at = m->begin(); at != m->end(); at++)
        mprintf(" %i", *at);
      mprintf("\n");
    }
  }
  return 0;
}

// Action_Radial::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
Action::RetType Action_Radial::Setup(ActionSetup& setup) {

  if ( setup.Top().SetupIntegerMask( Mask1_ ) ) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("Warning: First mask has no atoms.\n");
    return Action::SKIP;
  }
  if (setup.Top().SetupIntegerMask( Mask2_ ) ) return Action::ERR;
  if (Mask2_.None()) {
    mprintf("Warning: Second mask has no atoms.\n");
    return Action::SKIP;
  }
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );

  // If not computing center for mask 1 or 2, make the outer loop for distance
  // calculation correspond to the mask with the most atoms.
  if (rmode_ == NORMAL || rmode_ == NO_INTRAMOL) {
    if (Mask1_.Nselected() > Mask2_.Nselected()) {
      OuterMask_ = Mask1_;
      InnerMask_ = Mask2_;
    } else {
      OuterMask_ = Mask2_;
      InnerMask_ = Mask1_;
    }
  } else if (rmode_ == CENTER1) {
    OuterMask_ = Mask1_;
    InnerMask_ = Mask2_;
  } else if (rmode_ == CENTER2) {
    OuterMask_ = Mask2_;
    InnerMask_ = Mask1_;
  } else if (rmode_ == BYSITE) {
    // One or both masks will be by residue.
    int err = 0;
    if (siteMode1_ == BYRES)
      err = SetupSiteArrayByRes(Sites1_, setup.Top(), Mask1_);
    else if (siteMode1_ == BYMOL)
      err = SetupSiteArrayByMol(Sites1_, setup.Top(), Mask1_);
    else
      err = SetupSiteArrayByAtom(Sites1_, Mask1_);
    if (err != 0) return Action::ERR;
    if (siteMode2_ == BYRES)
      err = SetupSiteArrayByRes(Sites2_, setup.Top(), Mask2_);
    else if (siteMode2_ == BYMOL)
      err = SetupSiteArrayByMol(Sites2_, setup.Top(), Mask2_);
    else
      err = SetupSiteArrayByAtom(Sites2_, Mask2_);
    if (err != 0) return Action::ERR;
  } else {
    // SANITY CHECK
    mprinterr("Internal Error: Action_Radial: No mode set!\n");
    return Action::ERR;
  }
  // If ignoring intra-molecular distances, need to count how many we
  // are ignoring.
  if (rmode_ == NO_INTRAMOL) {
    int ndist = 0;
    for (AtomMask::const_iterator atom1 = OuterMask_.begin(); 
                                  atom1 != OuterMask_.end(); ++atom1)
      for (AtomMask::const_iterator atom2 = InnerMask_.begin();
                                    atom2 != InnerMask_.end(); ++atom2)
        if ( setup.Top()[*atom1].MolNum() == setup.Top()[*atom2].MolNum() )
          ++ndist;
    if (currentParm_ != 0 && ndist != intramol_distances_)
      mprintf("Warning: # of intramolecular distances (%i) has changed from the last"
              " topology (%i).\nWarning: Normalization will not be correct.\n",
              ndist, intramol_distances_);
    intramol_distances_ = ndist;
    currentParm_ = setup.TopAddress();
    mprintf("\tIgnoring %i intra-molecular distances.\n", intramol_distances_);
  }

  // Check volume information
  if (useVolume_ && setup.CoordInfo().TrajBox().Type()==Box::NOBOX) {
    mprintf("Warning: 'volume' specified but no box information for %s, skipping.\n",
            setup.Top().c_str());
    return Action::SKIP;
  }

  // Print mask and imaging info for this parm
  if (rmode_ == BYSITE) {
    mprintf("\t%zu sites selected by Mask1 (%i atoms), %zu sites selected by Mask2 (%i atoms)\n",
            Sites1_.size(), Mask1_.Nselected(), Sites2_.size(), Mask2_.Nselected());
  } else {
    mprintf("\t%i atoms in Mask1, %i atoms in Mask2\n",
            Mask1_.Nselected(), Mask2_.Nselected());
  }
  if (image_.ImagingEnabled())
    mprintf("\tImaging on.\n");
  else
    mprintf("\tImaging off.\n");
  return Action::OK;  
}

// Action_Radial::DoAction()
/** Calculate distances from atoms in mask1 to atoms in mask 2 and
  * bin them.
  */
// NOTE: Because of maximum2 not essential to check idx>numBins?
Action::RetType Action_Radial::DoAction(int frameNum, ActionFrame& frm) {
  double D;
  Matrix_3x3 ucell, recip;
  int atom1, atom2;
  int nmask1, nmask2;
  int idx;
# ifdef _OPENMP
  int mythread;
# endif

  // Set imaging information and store volume if specified
  // NOTE: Ucell and recip only needed for non-orthogonal boxes.
  if (image_.ImagingEnabled() || useVolume_) {
    D = frm.Frm().BoxCrd().ToRecip(ucell,recip);
    if (useVolume_)  volume_ += D;
  }
  // ---------------------------------------------
  if ( rmode_ == NORMAL ) { 
    // Calculation of all atoms in Mask1 to all atoms in Mask2
    int outer_max = OuterMask_.Nselected();
    int inner_max = InnerMask_.Nselected();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask1,nmask2,atom1,atom2,D,idx,mythread)
    {
    //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
      atom1 = OuterMask_[nmask1];
      for (nmask2 = 0; nmask2 < inner_max; nmask2++) {
        atom2 = InnerMask_[nmask2];
        if (atom1 != atom2) {
          D = DIST2( frm.Frm().XYZ(atom1), frm.Frm().XYZ(atom2),
                     image_.ImageType(), frm.Frm().BoxCrd(), ucell, recip);
          if (D <= maximum2_) {
            // NOTE: Can we modify the histogram to store D^2?
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing_);
            if (idx > -1 && idx < numBins_)
#             ifdef _OPENMP
              ++rdf_thread_[mythread][idx];
#             else
              ++RDF_[idx];
#             endif
          }
        }
      } // END loop over 2nd mask
    } // END loop over 1st mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  // ---------------------------------------------
  } else if ( rmode_ == NO_INTRAMOL ) {
    // Calculation of all atoms in Mask1 to all atoms in Mask2, ignoring
    // intra-molecular distances.
    int outer_max = OuterMask_.Nselected();
    int inner_max = InnerMask_.Nselected();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask1,nmask2,atom1,atom2,D,idx,mythread)
    {
    //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
      atom1 = OuterMask_[nmask1];
      for (nmask2 = 0; nmask2 < inner_max; nmask2++) {
        atom2 = InnerMask_[nmask2];
        if ( (*currentParm_)[atom1].MolNum() != (*currentParm_)[atom2].MolNum() ) {
          D = DIST2( frm.Frm().XYZ(atom1), frm.Frm().XYZ(atom2),
                     image_.ImageType(), frm.Frm().BoxCrd(), ucell, recip);
          if (D <= maximum2_) {
            // NOTE: Can we modify the histogram to store D^2?
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing_);
            if (idx > -1 && idx < numBins_)
#             ifdef _OPENMP
              ++rdf_thread_[mythread][idx];
#             else
              ++RDF_[idx];
#             endif
          }
        }
      } // END loop over 2nd mask
    } // END loop over 1st mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  // ---------------------------------------------
  } else if (rmode_ == BYSITE) {
    // Calculation of center of masks in Sites1 to center of masks in Sites2
    int mask1_max = (int)Sites1_.size();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask1,D,idx,mythread)
    {
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask1 = 0; nmask1 < mask1_max; nmask1++)
    {
      AtomMask const& site1 = Sites1_[nmask1];
      Vec3 com1 = frm.Frm().VGeometricCenter( site1 );
      for (Marray::const_iterator site2 = Sites2_.begin(); site2 != Sites2_.end(); ++site2)
      {
        if (site1 != *site2) {
          Vec3 com2 = frm.Frm().VGeometricCenter( *site2 );
          D = DIST2(com1.Dptr(), com2.Dptr(), image_.ImageType(),
                    frm.Frm().BoxCrd(), ucell, recip);
          if (D <= maximum2_) {
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing_);
            if (idx > -1 && idx < numBins_)
#             ifdef _OPENMP
              ++rdf_thread_[mythread][idx];
#             else
              ++RDF_[idx];
#             endif
          }
        } // END site1 != site2
      } // END inner loop over Sites2
    } // END outer loop over Sites1
#   ifdef _OPENMP
    }
#   endif
  // ---------------------------------------------
  } else { // CENTER1 || CENTER2
    // Calculation of center of one Mask to all atoms in other Mask
    Vec3 coord_center = frm.Frm().VGeometricCenter(OuterMask_);
    int mask2_max = InnerMask_.Nselected();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask2,atom2,D,idx,mythread)
    {
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask2 = 0; nmask2 < mask2_max; nmask2++) {
      atom2 = InnerMask_[nmask2];
      D = DIST2(coord_center.Dptr(), frm.Frm().XYZ(atom2), image_.ImageType(),
                frm.Frm().BoxCrd(), ucell, recip);
      if (D <= maximum2_) {
        // NOTE: Can we modify the histogram to store D^2?
        D = sqrt(D);
        //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
        idx = (int) (D * one_over_spacing_);
        if (idx > -1 && idx < numBins_)
#         ifdef _OPENMP
          ++rdf_thread_[mythread][idx];
#         else
          ++RDF_[idx];
#         endif
      }
    } // END loop over 2nd mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif 
  }
  ++numFrames_;

  return Action::OK;
} 

#ifdef MPI
int Action_Radial::SyncAction() {
# ifdef _OPENMP
  CombineRdfThreads();
# endif
  int total_frames = 0;
  trajComm_.ReduceMaster( &total_frames, &numFrames_, 1, MPI_INT, MPI_SUM );
  if (trajComm_.Master()) {
    numFrames_ = total_frames;
    int* sum_bins = new int[ numBins_ ];
    trajComm_.ReduceMaster( sum_bins, RDF_, numBins_, MPI_INT, MPI_SUM );
    std::copy( sum_bins, sum_bins + numBins_, RDF_ );
    delete[] sum_bins;
  } else
    trajComm_.ReduceMaster( 0,        RDF_, numBins_, MPI_INT, MPI_SUM );
  return 0;
}
#endif

#ifdef _OPENMP
/** Combine results from each rdf_thread into rdf. */
void Action_Radial::CombineRdfThreads() {
  if (rdf_thread_ == 0) return;
  for (int thread = 0; thread < numthreads_; thread++) { 
    for (int bin = 0; bin < numBins_; bin++)
      RDF_[bin] += rdf_thread_[thread][bin];
    delete[] rdf_thread_[thread];
  }
  delete[] rdf_thread_;
  rdf_thread_ = 0;
}
#endif

// Action_Radial::Print()
/** Convert the histogram to a dataset, normalize, create datafile.
  */
// NOTE: Currently the normalization is based on number of atoms in each mask;
//       if multiple prmtops are loaded and number of atoms changes from 
//       prmtop to prmtop this will throw off normalization.
void Action_Radial::Print() {
  if (numFrames_==0) return;
# ifdef _OPENMP
  CombineRdfThreads(); 
# endif
  mprintf("    RADIAL: %i frames,", numFrames_);
  double nmask1 = (double)Mask1_.Nselected();
  double nmask2 = (double)Mask2_.Nselected();
  int numSameAtoms = 0;
  if (rmode_ == NORMAL) {
    // If Mask1 and Mask2 have any atoms in common distances were not calcd
    // between them (because they are 0.0 of course); need to correct for this.
    numSameAtoms = Mask1_.NumAtomsInCommon( Mask2_ );
  } else if (rmode_ == NO_INTRAMOL) {
    // # of intra-molecular distances already counted. Any atoms in common
    // between mask1 and mask2 will have been included as an intramol dist.
    numSameAtoms = intramol_distances_ ;
  } else if (rmode_ == CENTER1) {
    // If the center1 option was specified only one distance was calcd
    // from mask 1. Assume COM of mask 1 != atom(s) in mask2.
    nmask1 = 1.0;
    numSameAtoms = 0;
  } else if (rmode_ == CENTER2) {
    // If the center2 option was specified only one distance was calcd
    // from mask 2. Assume COM of mask 2 != atom(s) in mask1.
    nmask2 = 1.0;
    numSameAtoms = 0;
  } else if (rmode_ == BYSITE) {
    // Count sites in common
    nmask1 = (double)Sites1_.size();
    nmask2 = (double)Sites2_.size();
    numSameAtoms = 0;
    for (Marray::const_iterator site1 = Sites1_.begin(); site1 != Sites1_.end(); ++site1)
      for (Marray::const_iterator site2 = Sites2_.begin(); site2 != Sites2_.end(); ++site2)
        if (*site1 == *site2)
          numSameAtoms++;
  }
  mprintf(" # in mask1= %.0f, # in mask2 = %.0f, # in common = %i\n",
          nmask1, nmask2, numSameAtoms);
  
  // If useVolume, calculate the density from the average volume
  if (useVolume_) {
    double avgVol = volume_ / numFrames_;
    mprintf("\tAverage volume is %f Ang^3.\n",avgVol);
    density_ = (nmask1 * nmask2 - (double)numSameAtoms) / avgVol;
    mprintf("\tAverage density is %f distances / Ang^3.\n",density_);
  } else {
    density_ = density_ * 
               (nmask1 * nmask2 - (double)numSameAtoms) / nmask1;
    mprintf("\tDensity is %f distances / Ang^3.\n",density_);
  }
  // Need to normalize each bin, which holds the particle count at that
  // distance. Calculate the expected number of molecules for that 
  // volume slice. Expected # of molecules is particle density times volume 
  // of each slice:
  // Density * ( [(4/3)*PI*(R+dr)^3] - [(4/3)*PI*(dr)^3] )
  double sum = 0.0;
  for (int bin = 0; bin < numBins_; bin++) {
    //mprintf("DBG:\tNumBins= %i\n",rdf[bin]); 
    // Number of particles in this volume slice over all frames.
    double N = (double) RDF_[bin];
    if (rawrdf_ != 0)
      rawrdf_->Add(bin, &N);
    // r, r + dr
    double R = spacing_ * (double)bin;
    double Rdr = R + spacing_;
    // Volume of slice: 4/3_pi * [(r+dr)^3 - (dr)^3]
    double dv = Constants::FOURTHIRDSPI * ( (Rdr * Rdr * Rdr) - (R * R * R) );
    // Expected # distances in this volume slice
    double expectedD = dv * density_;
    if (debug_>0)
      mprintf("    \tBin %f->%f <Pop>=%f, V=%f, D=%f, norm %f distances.\n",
              R,Rdr,N/numFrames_,dv,density_,expectedD);
    // Divide by # frames
    double norm = expectedD * (double)numFrames_;
    N /= norm;
    Dset_->Add(bin, &N);
    // If specified, calc integral of # mask2 atoms as fn of distance
    if (intrdf_ != 0) {
      sum += N * expectedD / nmask2;
      intrdf_->Add(bin, &sum);
    }
  }
}
