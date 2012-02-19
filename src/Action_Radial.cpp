// Radial
#include <cmath> // sqrt
#include <cstring> // memset
#include "Action_Radial.h"
#include "CpptrajStdio.h"
#include "Constants.h" // FOURTHIRDSPI
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Radial::Radial() {
  //fprintf(stderr,"Radial Con\n");
  rdf=NULL;
  rdf_thread=NULL;
  numthreads=1;
  useImage=true;
  center1=false;
  useVolume=false;
  volume=0;
  outfilename=NULL;
  maximum=0;
  maximum2=0;
  spacing=-1;
  one_over_spacing=-1;
  numBins=0;
  numFrames=0;
  numDistances=0;
  // Default particle density (molecules/Ang^3) for water based on 1.0 g/mL
  density = 0.033456;
} 

// DESTRUCTOR
Radial::~Radial() {
  //fprintf(stderr,"Radial Destructor.\n");
  if (rdf!=NULL) delete[] rdf;
  if (rdf_thread!=NULL) {
    for (int i=0; i < numthreads; i++)
      delete[] rdf_thread[i];
    delete[] rdf_thread;
  }
}

// Radial::init()
/** Expected call: radial <outfilename> <spacing> <maximum> <mask1> [<mask2>] [noimage]
  *                       [density <density> | volume] [center1] 
  */
int Radial::init() {
  char *mask1, *mask2;

  // Get Keywords
  useImage = !(actionArgs.hasKey("noimage"));
  // Default particle density (mols/Ang^3) for water based on 1.0 g/mL
  density = actionArgs.getKeyDouble("density",0.033456);
  center1 = actionArgs.hasKey("center1");
  useVolume = actionArgs.hasKey("volume");

  // Get required args
  outfilename = actionArgs.getNextString();
  if (outfilename==NULL) {
    mprinterr("Error: Radial: No output filename given.\n");
    return 1;
  }
  spacing = actionArgs.getNextDouble(-1.0);
  if (spacing < 0) {
    mprinterr("Error: Radial: No spacing argument or arg < 0.\n");
    return 1;
  }
  maximum = actionArgs.getNextDouble(-1.0);
  if (maximum < 0) {
    mprinterr("Error: Radial: No maximum argument or arg < 0.\n");
    return 1;
  }
  // Store max^2, distances^2 greater than max^2 do not need to be
  // binned and therefore do not need a sqrt calc.
  maximum2 = maximum * maximum;

  // Get First Mask
  mask1 = actionArgs.getNextMask();
  if (mask1==NULL) {
    mprinterr("Error: Radial: No mask given.\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);

  // Check for second mask - if none specified use first mask
  mask2 = actionArgs.getNextMask();
  if (mask2!=NULL) 
    Mask2.SetMaskString(mask2);
  else
    Mask2.SetMaskString(mask1);

  // Set up histogram
  one_over_spacing = 1 / spacing;
  double temp_numbins = maximum * one_over_spacing;
  temp_numbins = ceil(temp_numbins);
  numBins = (int) temp_numbins;
  rdf = new int[ numBins ];
  memset(rdf, 0, numBins * sizeof(int));
# ifdef _OPENMP
  // Since rdf is shared by all threads and we cant guarantee that a given
  // bin in rdf wont be accessed at the same time by the same thread,
  // each thread needs its own bin space.
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    numthreads = omp_get_num_threads();
}
  rdf_thread = new int*[ numthreads ];
  for (int i=0; i < numthreads; i++) {
    rdf_thread[i] = new int[ numBins ];
    memset(rdf_thread[i], 0, numBins * sizeof(int));
  }
# endif
  
  mprintf("    RADIAL: Calculating RDF for atoms in mask [%s]",Mask1.MaskString());
  if (mask2!=NULL) 
    mprintf(" to atoms in mask [%s]",Mask2.MaskString());
  mprintf("\n            Output to %s.\n",outfilename);
  if (center1)
    mprintf("            Using center of atoms in mask1.\n");
  //mprintf("            Histogram max %lf, spacing %lf, bins %i.\n",rdf.Max(0),
  //        rdf.Step(0),rdf.NBins(0));
  mprintf("            Histogram max %lf, spacing %lf, bins %i.\n",maximum,
          spacing,numBins);
  if (useVolume)
    mprintf("            Normalizing based on cell volume.\n");
  else
    mprintf("            Normalizing using particle density of %lf molecules/Ang^3.\n",density);
  if (!useImage) 
    mprintf("            Imaging disabled.\n");
  if (numthreads > 1)
    mprintf("            Parallelizing RDF calculation with %i threads.\n",numthreads);

  return 0;
}

// Radial::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
int Radial::setup() {

  if ( currentParm->SetupIntegerMask( Mask1, activeReference) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Radial::setup: Masks has no atoms.\n");
    return 1;
  }
  if (currentParm->SetupIntegerMask( Mask2, activeReference) ) return 1;
  if (Mask2.None()) {
    mprintf("    Error: Radial::setup: Second mask has no atoms.\n");
    return 1;
  }

  // If not computing a center for mask 1, make the outer loop for distance 
  // calculation correspond to the mask with the most atoms.
  if (!center1) {
    if (Mask1.Nselected > Mask2.Nselected) {
      OuterMask = Mask1;
      InnerMask = Mask2;
    } else {
      OuterMask = Mask2;
      InnerMask = Mask1;
    }
  }

  // Check volume information
  if (useVolume && currentParm->boxType==NOBOX) {
    mprintf("    Warning: Radial: 'volume' specified but no box information for %s, skipping.\n",
            currentParm->parmName);
    return 1;
  }

  // Print mask and imaging info for this parm
  mprintf("    RADIAL: %i atoms in Mask1, %i atoms in Mask2, ",Mask1.Nselected,Mask2.Nselected);
  if (imageType !=NOBOX)
    mprintf("Imaging on.\n");
  else
    mprintf("Imaging off.\n");
        
  return 0;  
}

// Radial::action()
/** Calculate distances from atoms in mask1 to atoms in mask 2 and
  * bin them.
  */
// NOTE: Because of maximum2 not essential to check idx>numBins?
int Radial::action() {
  double D, ucell[9], recip[9], coord_center[3];
  int atom1, atom2;
  int nmask1, nmask2;
  int idx, mydistances;
# ifdef _OPENMP
  int mythread;
# endif

  // Set imaging information and store volume if specified
  // NOTE: Ucell and recip only needed for non-orthogonal boxes.
  if (imageType!=NOBOX) {
    D = currentFrame->BoxToRecip(ucell,recip);
    if (useVolume)  volume += D;
  }

  mydistances = 0;
  // Calculation of center of Mask1 to all atoms in Mask2
  if (center1) {
    currentFrame->GeometricCenter(&Mask1,coord_center);
#ifdef _OPENMP
#pragma omp parallel private(nmask2,atom2,D,idx,mythread) reduction(+:mydistances)
{
  mythread = omp_get_thread_num();
#pragma omp for
#endif
    for (nmask2 = 0; nmask2 < Mask2.Nselected; nmask2++) {
      atom2 = Mask2.Selected[nmask2];
      D = currentFrame->DIST2(coord_center,atom2,imageType,ucell,recip);
      if (D <= maximum2) {
        // NOTE: Can we modify the histogram to store D^2?
        D = sqrt(D);
        //mprintf("MASKLOOP: %10i %10i %10.4lf\n",atom1,atom2,D);
        idx = (int) (D * one_over_spacing);
        if (idx > -1 && idx < numBins)
#         ifdef _OPENMP
          ++rdf_thread[mythread][idx];
#         else
          ++rdf[idx];
#         endif
        ++mydistances;
      }
    } // END loop over 2nd mask
#ifdef _OPENMP
} // END pragma omp parallel
#endif 
  // Calculation of all atoms in Mask1 to all atoms in Mask2
  } else {
#ifdef _OPENMP
#pragma omp parallel private(nmask1,nmask2,atom1,atom2,D,idx,mythread) reduction(+:mydistances) 
{
  //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
  mythread = omp_get_thread_num();
#pragma omp for
#endif
    for (nmask1 = 0; nmask1 < OuterMask.Nselected; nmask1++) {
      atom1 = OuterMask.Selected[nmask1];
      for (nmask2 = 0; nmask2 < InnerMask.Nselected; nmask2++) {
        atom2 = InnerMask.Selected[nmask2];
        if (atom1 != atom2) {
          D = currentFrame->DIST2(atom1,atom2,imageType,ucell,recip);
          if (D <= maximum2) {
            // NOTE: Can we modify the histogram to store D^2?
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4lf\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing);
            if (idx > -1 && idx < numBins)
#             ifdef _OPENMP
              ++rdf_thread[mythread][idx];
#             else
              ++rdf[idx];
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
  
  numDistances += mydistances;
  ++numFrames;

  return 0;
} 

// Radial::print()
/** Convert the histogram to a dataset, normalize, create datafile.
  */
// NOTE: Currently the normalization is based on number of atoms in each mask;
//       if multiple prmtops are loaded and number of atoms changes from 
//       prmtop to prmtop this will throw off normalization.
void Radial::print() {
  DataFile *outfile;
  DataSet *Dset;
  std::string xlabel;
  double N, R, Rdr, dv, norm;
 
  if (numFrames==0) return;
#  ifdef _OPENMP 
  // Combine results from each rdf_thread into rdf
  for (int thread=0; thread < numthreads; thread++) 
    for (int bin = 0; bin < numBins; bin++) 
      rdf[bin] += rdf_thread[thread][bin];
#  endif

  // Set up output dataset. 
  Dset = DSL->Add( DOUBLE, NULL, "g(r)");
  outfile = DFL->Add(outfilename, Dset);
  if (outfile==NULL) {
    mprinterr("Error: Radial: Could not setup output file %s\n",outfilename);
    return;
  }
  // Make default precision a little higher than normal
  Dset->SetPrecision(12,6);

  mprintf("    RADIAL: %i frames, %i distances.\n",numFrames,numDistances);

  // If Mask1 and Mask2 have any atoms in common distances were not calcd
  // between them (because they are 0.0 of course); need to correct for this.
  int numSameAtoms = Mask1.NumAtomsInCommon( Mask2 );
  
  // If useVolume, calculate the density from the average volume
  if (useVolume) {
    dv = volume / numFrames;
    mprintf("            Average volume is %lf Ang^3.\n",dv);
    density = ((double)Mask1.Nselected * (double)Mask2.Nselected - (double)numSameAtoms) / dv;
    mprintf("            Average density is %lf distances / Ang^3.\n",density);
  } else {
    density = density * ((double)Mask1.Nselected * (double)Mask2.Nselected - (double)numSameAtoms) /
              ((double)Mask1.Nselected);
    mprintf("            Density is %lf distances / Ang^3.\n",density);
  }

  // Need to normalize each bin, which holds the particle count at that
  // distance. Calculate the expected number of molecules for that 
  // volume slice. Expected # of molecules is particle density times volume 
  // of each slice:
  // Density * ( [(4/3)*PI*(R+dr)^3] - [(4/3)*PI*(dr)^3] )
  for (int bin = 0; bin < numBins; bin++) {
    //mprintf("DBG:\tNumBins= %i\n",rdf[bin]); 
    // Number of particles in this volume slice over all frames.
    N = (double) rdf[bin];
    // r
    R = spacing * (double)bin;
    // r + dr
    Rdr = spacing * (double)(bin+1);
    // Volume of slice: 4/3_pi * [(r+dr)^3 - (dr)^3]
    dv = FOURTHIRDSPI * ( (Rdr * Rdr * Rdr) - (R * R * R) );
    // Expected # distances in this volume slice
    norm = dv * density;
    if (debug>0)
      mprintf("    \tBin %lf->%lf <Pop>=%lf, V=%lf, D=%lf, norm %lf distances.\n",
              R,Rdr,N/numFrames,dv,density,norm);
    // Divide by # frames
    norm *= numFrames;
    N /= norm;

    Dset->Add(bin,&N);
  }
  // Setup output datafile. Create label from mask strings
  xlabel = '[' + Mask1.MaskExpression() + "] => [" + Mask2.MaskExpression() + ']';
  outfile->SetCoordMinStep(0.0,spacing,0,-1);
  outfile->SetXlabel((char*)xlabel.c_str());
  outfile->SetYlabel((char*)"g(r)");
}

