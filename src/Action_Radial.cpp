// Radial
#include "Action_Radial.h"
#include "CpptrajStdio.h"
#include <cmath> // sqrt
#include <cstdio> // sprintf
#include "vectormath.h" // PI

// CONSTRUCTOR
Radial::Radial() {
  //fprintf(stderr,"Radial Con\n");
  noimage=false;
  imageType=0;
  hasSecondMask=false;
  spacing=-1;
  maximum=0;
  maximum2=0;
  numFrames=0;
  numDistances=0;
} 

// DESTRUCTOR
Radial::~Radial() {
  //fprintf(stderr,"Radial Destructor.\n");
}

/* Radial::init()
 * Expected call: radial <name> <spacing> <maximum> <mask1> [<mask2>] [noimage]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Radial::init() {
  char *mask1, *mask2;

  // Get Keywords
  noimage = A->hasKey("noimage");

  // Get required args
  outfilename = A->getNextString();
  if (outfilename==NULL) {
    mprinterr("Error: Radial: No output filename given.\n");
    return 1;
  }
  spacing = A->getNextDouble(-1.0);
  if (spacing < 0) {
    mprinterr("Error: Radial: No spacing argument or arg < 0.\n");
    return 1;
  }
  maximum = A->getNextDouble(-1.0);
  if (maximum < 0) {
    mprinterr("Error: Radial: No maximum argument or arg < 0.\n");
    return 1;
  }
  maximum2 = maximum * maximum;

  // Get First Mask
  mask1 = A->getNextMask();
  if (mask1==NULL) {
    mprinterr("Error: Radial: No mask given.\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);

  // Check for second mask - if none specified use first mask
  mask2 = A->getNextMask();
  if (mask2!=NULL) {
    hasSecondMask=true;
    Mask2.SetMaskString(mask2);
  } else
    Mask2.SetMaskString(mask1);

  // Set up histogram
  rdf.SetDebug(debug);
  if (rdf.AddDimension((char*)"RDF", 0.0, maximum, spacing)) {
    mprinterr("Error: Radial: Could not set up histogram for RDF.\n");
    return 1;
  }

  mprintf("    RADIAL: Calculating RDF for atoms in mask [%s]",Mask1.maskString);
  if (hasSecondMask) 
    mprintf(" to atoms in mask [%s]",Mask2.maskString);
  mprintf(", output to %s.\n",outfilename);
  mprintf("            Histogram max %lf, spacing %lf, bins %i.\n",rdf.Max(0),
          rdf.Step(0),rdf.NBins(0));
  if (noimage) 
    mprintf("            Imaging disabled.\n");

  return 0;
}

/* Radial::setup()
 * Determine what atoms each mask pertains to for the current parm file.
 * Also determine whether imaging should be performed.
 */
int Radial::setup() {

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Radial::setup: Masks has no atoms.\n");
    return 1;
  }
  if ( Mask2.SetupMask(P,debug) ) return 1;
  if (Mask2.None()) {
    mprintf("    Error: Radial::setup: Second mask has no atoms.\n");
    return 1;
  }

  // Check imaging - check box based on prmtop box
  imageType = 0;
  if (!noimage) {
    imageType = (int)P->boxType;
    if (P->boxType==NOBOX && debug>0) {
      mprintf("    Warning: No box info in %s, disabling imaging.\n",P->parmName);
    }
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
  double D, ucell[9], recip[9];
  int atom1, atom2;

  if (imageType>0) F->BoxToRecip(ucell,recip);
  for (int nmask1 = 0; nmask1 < Mask1.Nselected; nmask1++) {
    atom1 = Mask1.Selected[nmask1];
    for (int nmask2 = 0; nmask2 < Mask2.Nselected; nmask2++) {
      atom2 = Mask2.Selected[nmask2];
      if (atom1 != atom2) {
        D = F->DIST2(atom1,atom2,imageType,ucell,recip);
        if (D > maximum2) continue;
        // NOTE: Can we modify the histogram to store D^2?
        D = sqrt(D);
        //fprintf(outfile,"%10i %10.4lf\n",currentFrame,D);
        rdf.BinData(&D);
        numDistances++;
      }
    } // END loop over 2nd mask
  } // END loop over 1st mask

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
  double four_thirds_pi = (4/3) * PI;
  double R, Rdr, dv, norm;
  double N;
  double avgDistances;
  double density = 0.033456;
  char temp[128];
  
  // Create label from mask strings
  sprintf(temp,"[%s] => [%s]",Mask1.maskString,Mask2.maskString);

  Dset = rdfdata.Add( DOUBLE, (char*)"RDF", "RDF");
  outfile = DFL->Add(outfilename, Dset);
  if (outfile==NULL) {
    mprinterr("Error: Radial: Could not setup output file %s\n",outfilename);
    return;
  }

  mprintf("    RADIAL: %i frames, %i distances.\n",numFrames,numDistances);
  avgDistances = (double) numDistances;
  avgDistances /= (double) numFrames;
  mprintf("    RADIAL: Avg number of distances = %lf\n",avgDistances);
  mprintf("    RADIAL: Histogram has %.0lf values.\n",rdf.BinTotal());

  rdf.BinStart(false);
  bin = 0;
  while (histloop) {
    N = rdf.CurrentBinData();
    // R
    rdf.CurrentBinCoord(&R);
    // R + dr
    Rdr = R + rdf.Step(0);
    // 4/3 pi * density * [(r+dr)^3 - r^3]
    dv = four_thirds_pi * density * ( (Rdr * Rdr * Rdr) - (R * R * R) );
    // normalized w.r.t. the average # of particles binned per frame
    norm = dv * numDistances;
    N /= norm;

    Dset->Add(bin,&N);
    bin++;
    if (rdf.NextBin()) histloop=false;
  }
  outfile->SetCoordMinStep(rdf.Min(0),rdf.Step(0),0,-1);
  outfile->SetXlabel(temp);
  outfile->SetYlabel((char*)"g(r)");
}

  
