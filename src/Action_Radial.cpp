// Action_Radial
#include <cmath> // sqrt
#include <cstring> // memset
#include "Action_Radial.h"
#include "CpptrajStdio.h"
#include "Constants.h" // FOURTHIRDSPI
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_Radial::Action_Radial() :
  RDF_(NULL),
  rdf_thread_(NULL),
  center1_(false),
  useVolume_(false),
  volume_(0),
  outfilename_(NULL),
  maximum2_(0),
  spacing_(-1),
  one_over_spacing_(-1),
  numBins_(0),
  numthreads_(1),
  numFrames_(0),
  numDistances_(0),
  // Default particle density (molecules/Ang^3) for water based on 1.0 g/mL
  density_(0.033456)
{
  //fprintf(stderr,"Radial Con\n");
  useImage_ = true;
} 

// DESTRUCTOR
Action_Radial::~Action_Radial() {
  //fprintf(stderr,"Radial Destructor.\n");
  if (RDF_!=NULL) delete[] RDF_;
  if (rdf_thread_!=NULL) {
    for (int i=0; i < numthreads_; i++)
      delete[] rdf_thread_[i];
    delete[] rdf_thread_;
  }
}

// Action_Radial::init()
/** Expected call: radial <outfilename> <spacing> <maximum> <mask1> [<mask2>] [noimage]
  *                       [density <density> | volume] [center1] 
  */
int Action_Radial::init() {
  char *mask1, *mask2;

  // Get Keywords
  useImage_ = !(actionArgs.hasKey("noimage"));
  // Default particle density (mols/Ang^3) for water based on 1.0 g/mL
  density_ = actionArgs.getKeyDouble("density",0.033456);
  center1_ = actionArgs.hasKey("center1");
  useVolume_ = actionArgs.hasKey("volume");

  // Get required args
  outfilename_ = actionArgs.getNextString();
  if (outfilename_==NULL) {
    mprinterr("Error: Radial: No output filename given.\n");
    return 1;
  }
  spacing_ = actionArgs.getNextDouble(-1.0);
  if (spacing_ < 0) {
    mprinterr("Error: Radial: No spacing argument or arg < 0.\n");
    return 1;
  }
  double maximum = actionArgs.getNextDouble(-1.0);
  if (maximum < 0) {
    mprinterr("Error: Radial: No maximum argument or arg < 0.\n");
    return 1;
  }
  // Store max^2, distances^2 greater than max^2 do not need to be
  // binned and therefore do not need a sqrt calc.
  maximum2_ = maximum * maximum;

  // Get First Mask
  mask1 = actionArgs.getNextMask();
  if (mask1==NULL) {
    mprinterr("Error: Radial: No mask given.\n");
    return 1;
  }
  Mask1_.SetMaskString(mask1);

  // Check for second mask - if none specified use first mask
  mask2 = actionArgs.getNextMask();
  if (mask2!=NULL) 
    Mask2_.SetMaskString(mask2);
  else
    Mask2_.SetMaskString(mask1);

  // Set up histogram
  one_over_spacing_ = 1 / spacing_;
  double temp_numbins = maximum * one_over_spacing_;
  temp_numbins = ceil(temp_numbins);
  numBins_ = (int) temp_numbins;
  RDF_ = new int[ numBins_ ];
  memset(RDF_, 0, numBins_ * sizeof(int));
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
    memset(rdf_thread_[i], 0, numBins_ * sizeof(int));
  }
# endif
  
  mprintf("    RADIAL: Calculating RDF for atoms in mask [%s]",Mask1_.MaskString());
  if (mask2!=NULL) 
    mprintf(" to atoms in mask [%s]",Mask2_.MaskString());
  mprintf("\n            Output to %s.\n",outfilename_);
  if (center1_)
    mprintf("            Using center of atoms in mask1.\n");
  //mprintf("            Histogram max %lf, spacing %lf, bins %i.\n",rdf.Max(0),
  //        rdf.Step(0),rdf.NBins(0));
  mprintf("            Histogram max %lf, spacing %lf, bins %i.\n",maximum,
          spacing_,numBins_);
  if (useVolume_)
    mprintf("            Normalizing based on cell volume.\n");
  else
    mprintf("            Normalizing using particle density of %lf molecules/Ang^3.\n",density_);
  if (!useImage_) 
    mprintf("            Imaging disabled.\n");
  if (numthreads_ > 1)
    mprintf("            Parallelizing RDF calculation with %i threads.\n",numthreads_);

  return 0;
}

// Action_Radial::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
int Action_Radial::setup() {

  if ( currentParm->SetupIntegerMask( Mask1_ ) ) return 1;
  if (Mask1_.None()) {
    mprintf("    Error: Radial::setup: Masks has no atoms.\n");
    return 1;
  }
  if (currentParm->SetupIntegerMask( Mask2_ ) ) return 1;
  if (Mask2_.None()) {
    mprintf("    Error: Radial::setup: Second mask has no atoms.\n");
    return 1;
  }

  // If not computing a center for mask 1, make the outer loop for distance 
  // calculation correspond to the mask with the most atoms.
  if (!center1_) {
    if (Mask1_.Nselected() > Mask2_.Nselected()) {
      OuterMask_ = Mask1_;
      InnerMask_ = Mask2_;
    } else {
      OuterMask_ = Mask2_;
      InnerMask_ = Mask1_;
    }
  }

  // Check volume information
  if (useVolume_ && currentParm->BoxType()==Box::NOBOX) {
    mprintf("    Warning: Radial: 'volume' specified but no box information for %s, skipping.\n",
            currentParm->c_str());
    return 1;
  }

  // Print mask and imaging info for this parm
  mprintf("    RADIAL: %i atoms in Mask1, %i atoms in Mask2, ",
          Mask1_.Nselected(), Mask2_.Nselected());
  if (imageType_ != Frame::NOIMAGE)
    mprintf("Imaging on.\n");
  else
    mprintf("Imaging off.\n");
        
  return 0;  
}

// Action_Radial::action()
/** Calculate distances from atoms in mask1 to atoms in mask 2 and
  * bin them.
  */
// NOTE: Because of maximum2 not essential to check idx>numBins?
int Action_Radial::action() {
  double D, ucell[9], recip[9], coord_center[3];
  int atom1, atom2;
  int nmask1, nmask2;
  int idx, mydistances;
# ifdef _OPENMP
  int mythread;
# endif

  // Set imaging information and store volume if specified
  // NOTE: Ucell and recip only needed for non-orthogonal boxes.
  if (imageType_ != Frame::NOIMAGE) {
    D = currentFrame->BoxToRecip(ucell,recip);
    if (useVolume_)  volume_ += D;
  }

  mydistances = 0;
  // Calculation of center of Mask1 to all atoms in Mask2
  if (center1_) {
    currentFrame->GeometricCenter(&Mask1_,coord_center);
    int mask2_max = Mask2_.Nselected();
#ifdef _OPENMP
#pragma omp parallel private(nmask2,atom2,D,idx,mythread) reduction(+:mydistances)
{
  mythread = omp_get_thread_num();
#pragma omp for
#endif
    for (nmask2 = 0; nmask2 < mask2_max; nmask2++) {
      atom2 = Mask2_[nmask2];
      D = currentFrame->DIST2(coord_center,atom2,imageType_,ucell,recip);
      if (D <= maximum2_) {
        // NOTE: Can we modify the histogram to store D^2?
        D = sqrt(D);
        //mprintf("MASKLOOP: %10i %10i %10.4lf\n",atom1,atom2,D);
        idx = (int) (D * one_over_spacing_);
        if (idx > -1 && idx < numBins_)
#         ifdef _OPENMP
          ++rdf_thread_[mythread][idx];
#         else
          ++RDF_[idx];
#         endif
        ++mydistances;
      }
    } // END loop over 2nd mask
#ifdef _OPENMP
} // END pragma omp parallel
#endif 
  // Calculation of all atoms in Mask1 to all atoms in Mask2
  } else {
    int outer_max = OuterMask_.Nselected();
    int inner_max = InnerMask_.Nselected();
#ifdef _OPENMP
#pragma omp parallel private(nmask1,nmask2,atom1,atom2,D,idx,mythread) reduction(+:mydistances) 
{
  //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
  mythread = omp_get_thread_num();
#pragma omp for
#endif
    for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
      atom1 = OuterMask_[nmask1];
      for (nmask2 = 0; nmask2 < inner_max; nmask2++) {
        atom2 = InnerMask_[nmask2];
        if (atom1 != atom2) {
          D = currentFrame->DIST2(atom1,atom2,imageType_,ucell,recip);
          if (D <= maximum2_) {
            // NOTE: Can we modify the histogram to store D^2?
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4lf\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing_);
            if (idx > -1 && idx < numBins_)
#             ifdef _OPENMP
              ++rdf_thread_[mythread][idx];
#             else
              ++RDF_[idx];
#             endif
            ++mydistances;
          }
        }
      } // END loop over 2nd mask
    } // END loop over 1st mask
#ifdef _OPENMP
} // END pragma omp parallel
#endif 
  } // END if center1
  
  numDistances_ += mydistances;
  ++numFrames_;

  return 0;
} 

// Action_Radial::print()
/** Convert the histogram to a dataset, normalize, create datafile.
  */
// NOTE: Currently the normalization is based on number of atoms in each mask;
//       if multiple prmtops are loaded and number of atoms changes from 
//       prmtop to prmtop this will throw off normalization.
void Action_Radial::print() {
  DataFile *outfile;
  DataSet *Dset;
  double N, R, Rdr, dv, norm;
 
  if (numFrames_==0) return;
#  ifdef _OPENMP 
  // Combine results from each rdf_thread into rdf
  for (int thread=0; thread < numthreads_; thread++) 
    for (int bin = 0; bin < numBins_; bin++) 
      RDF_[bin] += rdf_thread_[thread][bin];
#  endif

  // Set up output dataset. 
  Dset = DSL->Add( DataSet::DOUBLE, NULL, "g(r)");
  outfile = DFL->Add(outfilename_, Dset);
  if (outfile==NULL) {
    mprinterr("Error: Radial: Could not setup output file %s\n",outfilename_);
    return;
  }
  // Make default precision a little higher than normal
  Dset->SetPrecision(12,6);

  mprintf("    RADIAL: %i frames, %i distances.\n",numFrames_,numDistances_);

  // If Mask1 and Mask2 have any atoms in common distances were not calcd
  // between them (because they are 0.0 of course); need to correct for this.
  int numSameAtoms = Mask1_.NumAtomsInCommon( Mask2_ );
  
  // If useVolume, calculate the density from the average volume
  if (useVolume_) {
    dv = volume_ / numFrames_;
    mprintf("            Average volume is %lf Ang^3.\n",dv);
    density_ = ((double)Mask1_.Nselected() * (double)Mask2_.Nselected() - (double)numSameAtoms) / 
               dv;
    mprintf("            Average density is %lf distances / Ang^3.\n",density_);
  } else {
    density_ = density_ * 
              ((double)Mask1_.Nselected() * (double)Mask2_.Nselected() - (double)numSameAtoms) /
              ((double)Mask1_.Nselected());
    mprintf("            Density is %lf distances / Ang^3.\n",density_);
  }

  // Need to normalize each bin, which holds the particle count at that
  // distance. Calculate the expected number of molecules for that 
  // volume slice. Expected # of molecules is particle density times volume 
  // of each slice:
  // Density * ( [(4/3)*PI*(R+dr)^3] - [(4/3)*PI*(dr)^3] )
  for (int bin = 0; bin < numBins_; bin++) {
    //mprintf("DBG:\tNumBins= %i\n",rdf[bin]); 
    // Number of particles in this volume slice over all frames.
    N = (double) RDF_[bin];
    // r
    R = spacing_ * (double)bin;
    // r + dr
    Rdr = spacing_ * (double)(bin+1);
    // Volume of slice: 4/3_pi * [(r+dr)^3 - (dr)^3]
    dv = FOURTHIRDSPI * ( (Rdr * Rdr * Rdr) - (R * R * R) );
    // Expected # distances in this volume slice
    norm = dv * density_;
    if (debug>0)
      mprintf("    \tBin %lf->%lf <Pop>=%lf, V=%lf, D=%lf, norm %lf distances.\n",
              R,Rdr,N/numFrames_,dv,density_,norm);
    // Divide by # frames
    norm *= numFrames_;
    N /= norm;

    Dset->Add(bin,&N);
  }
  // Setup output datafile.
  outfile->ProcessArgs("xmin 0.0 xstep " + doubleToString(spacing_));
  // Create label from mask strings. Enclose in quotes so the label is 1 arg.
  outfile->ProcessArgs("xlabel \"[" + Mask1_.MaskExpression() + "] => [" + 
                                      Mask2_.MaskExpression() + "]\"" );
  outfile->ProcessArgs("ylabel g(r)");
}

