// Radial
#include "Action_Radial.h"
#include "CpptrajStdio.h"
#include <cmath> // sqrt
#include <cstdio> // sprintf
#include "Constants.h" // FOURTHIRDSPI

// CONSTRUCTOR
Radial::Radial() {
  //fprintf(stderr,"Radial Con\n");
  noimage=false;
  imageType=0;
  center1=false;
  useVolume=false;
  volume=0;
  outfilename=NULL;
  maximum=0;
  maximum2=0;
  spacing=-1;
  numFrames=0;
  numDistances=0;
  // Default particle density (mols/Ang^3) for water based on 1.0 g/mL
  density = 0.033456;
} 

// DESTRUCTOR
Radial::~Radial() {
  //fprintf(stderr,"Radial Destructor.\n");
}

/* Radial::init()
 * Expected call: radial <outfilename> <spacing> <maximum> <mask1> [<mask2>] [noimage]
 *                       [density <density> | volume] [center1] 
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Radial::init() {
  char *mask1, *mask2;

  // Get Keywords
  noimage = actionArgs.hasKey("noimage");
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
  rdf.SetDebug(debug);
  if (rdf.AddDimension((char*)"RDF", 0.0, maximum, spacing)) {
    mprinterr("Error: Radial: Could not set up histogram for RDF.\n");
    return 1;
  }

  mprintf("    RADIAL: Calculating RDF for atoms in mask [%s]",Mask1.MaskString());
  if (mask2!=NULL) 
    mprintf(" to atoms in mask [%s]",Mask2.MaskString());
  mprintf("\n            Output to %s.\n",outfilename);
  if (center1)
    mprintf("            Using center of atoms in mask1.\n");
  mprintf("            Histogram max %lf, spacing %lf, bins %i.\n",rdf.Max(0),
          rdf.Step(0),rdf.NBins(0));
  if (useVolume)
    mprintf("            Normalizing based on cell volume.\n");
  else
    mprintf("            Normalizing using particle density of %lf mols/Ang^3.\n",density);
  if (noimage) 
    mprintf("            Imaging disabled.\n");

  return 0;
}

/* Radial::setup()
 * Determine what atoms each mask pertains to for the current parm file.
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

  // Check imaging - check box based on prmtop box
  imageType = 0;
  if (!noimage) {
    imageType = (int)currentParm->boxType;
    if (currentParm->boxType==NOBOX && debug>0) {
      mprintf("    Warning: No box info in %s, disabling imaging.\n",currentParm->parmName);
    }
  }

  // Check volume information
  if (useVolume && currentParm->boxType==NOBOX) {
    mprintf("    Warning: Radial: 'volume' specified but no box information for %s, skipping.\n",
            currentParm->parmName);
    return 1;
  }

  // Print imaging info for this parm
  mprintf("    RADIAL: %i atoms in Mask1, %i atoms in Mask2, ",Mask1.Nselected,Mask2.Nselected);
  if (imageType > 0)
    mprintf("Imaging on.\n");
  else
    mprintf("Imaging off.\n");
        
  return 0;  
}

/* Radial::action()
 * Calculate distances from atoms in mask1 to atoms in mask 2 and
 * bin them.
 */
int Radial::action() {
  double D, ucell[9], recip[9], coord_center[3];
  int atom1, atom2;
  int nmask1, nmask2;

  // Set imaging information and store volume if specified
  // NOTE: Ucell and recip only needed for non-orthogonal boxes.
  if (imageType>0) {
    D = currentFrame->BoxToRecip(ucell,recip);
    if (useVolume)  volume += D;
  }

  if (center1) {
    currentFrame->GeometricCenter(&Mask1,coord_center);
    for (nmask2 = 0; nmask2 < Mask2.Nselected; nmask2++) {
      atom2 = Mask2.Selected[nmask2];
      D = currentFrame->DIST2(coord_center,atom2,imageType,ucell,recip);
      if (D > maximum2) continue;
      // NOTE: Can we modify the histogram to store D^2?
      D = sqrt(D);
      //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
      rdf.BinData(&D);
      numDistances++;
    } // END loop over 2nd mask

  } else {
//#ifdef _OPENMP
//#pragma omp parallel private(nmask1,nmask2,atom1,atom2,D)
//{
//  //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
//#pragma omp for
//#endif
  for (nmask1 = 0; nmask1 < Mask1.Nselected; nmask1++) {
    atom1 = Mask1.Selected[nmask1];
    for (nmask2 = 0; nmask2 < Mask2.Nselected; nmask2++) {
      atom2 = Mask2.Selected[nmask2];
      if (atom1 != atom2) {
        D = currentFrame->DIST2(atom1,atom2,imageType,ucell,recip);
        if (D > maximum2) continue;
        // NOTE: Can we modify the histogram to store D^2?
        D = sqrt(D);
        //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
        rdf.BinData(&D);
        numDistances++;
      }
    } // END loop over 2nd mask
  } // END loop over 1st mask
//#ifdef _OPENMP
//} // END pragma omp parallel
//#endif
  }

  numFrames++;

  return 0;
} 

/* Radial::print()
 * Convert the histogram to a dataset, normalize, create datafile.
 */
void Radial::print() {
  DataFile *outfile;
  DataSet *Dset;
  bool histloop=true;
  int bin;
  double R, Rdr, dv, norm;
  double N;
  char temp[128];
 
  if (numFrames==0) return;
 
  // Create label from mask strings
  sprintf(temp,"[%s] => [%s]",Mask1.MaskString(),Mask2.MaskString());

  Dset = rdfdata.Add( DOUBLE, (char*)"RDF", "RDF");
  outfile = DFL->Add(outfilename, Dset);
  if (outfile==NULL) {
    mprinterr("Error: Radial: Could not setup output file %s\n",outfilename);
    return;
  }

  mprintf("    RADIAL: %i frames, %i distances.\n",numFrames,numDistances);
  //mprintf("            Histogram has %.0lf values.\n",rdf.BinTotal());
  
  // If useVolume, calculate the density from the average volume
  if (useVolume) {
    dv = volume / numFrames;
    mprintf("            Average volume is %lf Ang^3.\n",dv);
    density = (Mask1.Nselected * Mask2.Nselected) / dv;
    mprintf("            Average density is %lf distances / Ang^3.\n",density);
  }

  // Need to normalize each bin, which holds the particle count at that
  // distance. Calculate the expected number of molecules for that 
  // volume slice. Expected # of molecules is particle density times volume 
  // of each slice:
  // Density * ( [(4/3)*PI*(R+dr)^3] - [(4/3)*PI*(dr)^3] )
  // If there is more than 1 atom in the first mask we have essentially
  // created Mask1.Nselected RDFs and are averaging them, so further divide
  // by Mask1.Nselected. 
  rdf.BinStart(false);
  bin = 0;
  while (histloop) {
    // Number of particles in this volume slice over all frames.
    N = rdf.CurrentBinData();
    // r
    rdf.CurrentBinCoord(&R);
    // r + dr
    Rdr = R + rdf.Step(0);
    // Volume of slice: 4/3_pi * [(r+dr)^3 - (dr)^3]
    dv = FOURTHIRDSPI * ( (Rdr * Rdr * Rdr) - (R * R * R) );
    // Expected # molecules in this volume slice
    norm = dv * density;
    if (debug>0)
      mprintf("    \tBin %lf->%lf %lf V %lf, D %lf, expect %lf molecules.\n",
              R,Rdr,N/numFrames,dv,density,norm);
    // DEBUG: Hack for selecting both Hs and O.
    //norm *= 3;
    // Divide by # frames, expected # of molecules, and # of RDFs
    norm *= numFrames;
    //norm *= Mask1.Nselected;
    N /= norm;

    Dset->Add(bin,&N);
    bin++;
    if (rdf.NextBin()) histloop=false;
  }
  outfile->SetCoordMinStep(rdf.Min(0),rdf.Step(0),0,-1);
  outfile->SetXlabel(temp);
  outfile->SetYlabel((char*)"g(r)");
}

  
