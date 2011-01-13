#include <cstdio> // sprintf
#include "Action_Hbond.h"
#include "CpptrajStdio.h"
// Hbond 

// CONSTRUCTOR
Hbond::Hbond() {
  Nframes=0; 
}

// DESTRUCTOR
Hbond::~Hbond() {
}

/*
 * Hbond::init()
 * Expected call: hbond [out <filename>] <mask> [angle <cut>] [dist <cut>] [avgout <filename>]
 * Search for Hbonding atoms in region specified by mask. 
 * Arg. check order is:
 *    1) Keywords
 *    2) Masks
 */
int Hbond::init() {
  char *mask, *outfilename;

  // Get keywords
  outfilename = A->getKeyString("out",NULL);
  avgout = A->getKeyString("avgout",NULL);
  acut = A->getKeyDouble("angle",135.0);
  dcut = A->getKeyDouble("dist",3.0);
  // Get masks
  mask = A->getNextMask();
  Mask.SetMaskString(mask);

  // Setup datasets
  NumHbonds = DSL->Add(INT, A->getNextString(),"NumHB");
  if (NumHbonds==NULL) return 1;
  DFL->Add(outfilename,NumHbonds);

  mprintf( "  HBOND: Calculating Hbonds in region specified by %s\n",Mask.maskString);
  mprintf( "         Distance cutoff = %8.3lf, Angle Cutoff = %8.3lf\n",dcut,acut);
  if (outfilename!=NULL) 
    mprintf( "         Dumping # Hbond v time results to %s\n", outfilename);
  if (avgout!=NULL)
    mprintf( "         Dumping Hbond avgs to %s\n",avgout);

  return 0;
}

/*
 * Hbond::setup()
 * Search for hbond donors and acceptors. 
 */
int Hbond::setup() {
  int atom, selected, a2;

  // Set up mask
  if ( Mask.SetupMask(P,debug) ) return 1;
  if ( Mask.None() ) {
    mprintf("    Error: Hbond::setup: Mask has no atoms.\n");
    return 1;
  }

  // Set up acceptors: F, O, N
  // NOTE: Attempt to determine electronegative carbons?
  for (selected=0; selected < Mask.Nselected; selected++) {
    atom = Mask.Selected[selected];
    if (P->names[atom][0]=='F' ||
        P->names[atom][0]=='O' ||
        P->names[atom][0]=='N'   )
      Acceptor.push_back(atom);
  }
  mprintf("      HBOND: Set up %i acceptors:\n",(int)Acceptor.size());
  for (accept = Acceptor.begin(); accept!=Acceptor.end(); accept++)
    mprintf("        %8i: %4s\n",*accept,P->names[*accept]);

  // Set up donors: O-H, N-H
  // NOTE: Scan donor list and determine which ones have H?
  for (accept = Acceptor.begin(); accept!=Acceptor.end(); accept++) {
    // Is this atom in the bondsh array? If so the other atom must be hydrogen
    for ( selected=0; selected < P->Nbonh()*3; selected+=3) {
      // Actual atom #s in bondsh array = x / 3
      atom = P->bondsh[selected  ] / 3;
      a2   = P->bondsh[selected+1] / 3;
      //mprintf("DEBUG: HBOND: Donor Setup: Accept=%i, selected=%i, atom=%i, a2=%i\n",
      //        *accept, selected, atom, a2);
      if (*accept == atom) {
        Donor.push_back(atom);
        Donor.push_back(a2);
      } else if (*accept == a2) {
        Donor.push_back(a2);
        Donor.push_back(atom);
      }
    }
  }
  mprintf("      HBOND: Set up %i donors:\n",((int)Donor.size())/2);
  for (donor = Donor.begin(); donor!=Donor.end(); donor++) {
    atom = (*donor);
    donor++;
    a2   = (*donor);
    mprintf("        %8i:%4s - %8i:%4s\n",atom,P->names[atom],a2,P->names[a2]);  
  }


  return 0;
}

/*
 * Hbond::action()
 * Calculate distance between all donors and acceptors. Store Hbond info.
 */    
int Hbond::action() {
  // accept ... H-D
  int D, H, Nhb, numHB;
  double dist, angle;
  std::map<int,HbondType>::iterator it;
  HbondType HB;

  Nhb = 0; numHB=0;
  for (donor = Donor.begin(); donor!=Donor.end(); donor++) {
    D = (*donor);
    donor++;
    H = (*donor);
    for (accept = Acceptor.begin(); accept!=Acceptor.end(); accept++, Nhb++) {
      if (*accept == D) continue;
      dist = F->DIST(*accept, D);
      if (dist > dcut) continue;
      angle = F->ANGLE(*accept, H, D);
      if (angle < acut) continue;
//      mprintf( "HBOND[%i]: %i:%s ... %i:%s-%i:%s Dist=%lf Angle=%lf\n", 
//              Nhb, *accept, P->names[*accept],
//              H, P->names[H], D, P->names[D], dist, angle);
      numHB++;
      // Find hbond in map
      it = HbondMap.find( Nhb );
      if (it == HbondMap.end() ) {
        // New Hbond
        HB.A=*accept;
        HB.D=D;
        HB.H=H;
        HB.Frames = 1;
        HB.dist=dist;
        HB.angle=angle;
        HbondMap.insert( it, std::pair<int,HbondType>(Nhb, HB) );
      } else {
        (*it).second.Frames++;
        (*it).second.dist+=dist;
        (*it).second.angle+=angle;
      }
    }
  }
  NumHbonds->Add(currentFrame, &numHB);
//  mprintf("HBOND: Scanned %i hbonds.\n",Nhb);
  Nframes++;

  return 0;
}

/*
 * Hbond::print()
 * Print average occupancies over all frames for all detected Hbonds
 */
void Hbond::print() {
  std::map<int,HbondType>::iterator it;
  std::list<HbondType> HbondList;
  std::list<HbondType>::iterator hbond;
  double avg, dist, angle;
  int Ares, Hres, Dres;
  char Aname[20], Hname[20], Dname[20], temp[20];
  PtrajFile OutFile;

  if (OutFile.SetupFile(avgout, WRITE, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug)) return;
  OutFile.OpenFile();

  OutFile.IO->Printf("HBONDS:\n");
  OutFile.IO->Printf("  %-15s %-15s %-15s %6s %6s %8s %8s\n",
                 "Acceptor","DonorH","Donor","Frames","Frac","AvgDist","AvgAng");
  for (it = HbondMap.begin(); it!=HbondMap.end(); it++) {
//    mprintf( "      %8i:%s ... %8i:%s-%8i:%s %8i Frames.\n",
//              (*it).second.A, P->names[(*it).second.A],
//              (*it).second.H, P->names[(*it).second.H],
//              (*it).second.D, P->names[(*it).second.D],
//              (*it).second.Frames);
    HbondList.push_back( (*it).second );
  }

  HbondList.sort( hbond_cmp() );
  for ( hbond = HbondList.begin(); hbond!=HbondList.end(); hbond++ ) {
    avg = (double) (*hbond).Frames;
    avg = avg / ((double) Nframes);
    dist = (double) (*hbond).dist;
    dist = dist / ((double) (*hbond).Frames);
    angle = (double) (*hbond).angle;
    angle = angle / ((double) (*hbond).Frames);
    Ares = P->atomToResidue((*hbond).A);
    Hres = P->atomToResidue((*hbond).H);
    Dres = P->atomToResidue((*hbond).D);
    P->ResName(temp,Ares);
    sprintf(Aname,"%s%i@%s",temp,Ares,P->names[(*hbond).A]);
    P->ResName(temp,Hres);
    sprintf(Hname,"%s%i@%s",temp,Hres,P->names[(*hbond).H]);
    P->ResName(temp,Dres);
    sprintf(Dname,"%s%i@%s",temp,Dres,P->names[(*hbond).D]);
    OutFile.IO->Printf("%-15s %-15s %-15s %6i %6.2lf %8.3lf %8.3lf\n",
                       Aname,Hname,Dname, (*hbond).Frames,avg,dist,angle);
//    OutFile.IO->Printf("  %6i:%4s %6i:%4s %6i:%4s %6i %6.2lf %8.3lf %8.3lf\n",
//              (*hbond).A, P->names[(*hbond).A],
//              (*hbond).H, P->names[(*hbond).H],
//              (*hbond).D, P->names[(*hbond).D],
//              (*hbond).Frames,avg,dist,angle);
  }
  OutFile.CloseFile();
}
