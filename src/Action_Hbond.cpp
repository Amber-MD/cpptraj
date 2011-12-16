// Hbond 
#include <cmath> // sqrt
#include "Action_Hbond.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG, DEGRAD

// CONSTRUCTOR
Hbond::Hbond() {
  Nframes=0;
  HBavg=NULL;
  hasDonorMask=false;
  hasAcceptorMask=false; 
}

// DESTRUCTOR
Hbond::~Hbond() {
  if (HBavg!=NULL) delete HBavg;
}

// Hbond::init()
/** Expected call: hbond [out <filename>] <mask> [angle <cut>] [dist <cut>] [avgout <filename>]
  *                      [donormask <mask>] [acceptormask <mask>]
  * Search for Hbonding atoms in region specified by mask. 
  * Arg. check order is:
  * - Keywords
  * - Masks
  * If just <mask> is specified donors and acceptors will be automatically
  * searched for.
  * If donormask is specified but not acceptormask, acceptors will be 
  * automatically searched for in <mask>.
  * If acceptormask is specified but not donormask, donors will be automatically
  * searched for in <mask>.
  * If both donormask and acceptor mask are specified no searching will occur.
  */
int Hbond::init() {
  char *mask, *outfilename;

  // Get keywords
  outfilename = actionArgs.getKeyString("out",NULL);
  avgout = actionArgs.getKeyString("avgout",NULL);
  acut = actionArgs.getKeyDouble("angle",135.0);
  // Convert angle cutoff to radians
  acut *= DEGRAD;
  dcut = actionArgs.getKeyDouble("dist",3.0);
  dcut2 = dcut * dcut;
  // Get donor mask
  mask = actionArgs.getKeyString("donormask",NULL);
  if (mask!=NULL) {
    DonorMask.SetMaskString(mask);
    hasDonorMask=true;
  }
  // Get acceptor mask
  mask = actionArgs.getKeyString("acceptormask",NULL);
  if (mask!=NULL) {
    AcceptorMask.SetMaskString(mask);
    hasAcceptorMask=true;
  }
  // Get generic mask
  mask = actionArgs.getNextMask();
  Mask.SetMaskString(mask);

  // Setup datasets
  NumHbonds = DSL->Add(INT, actionArgs.getNextString(),"NumHB");
  if (NumHbonds==NULL) return 1;
  DFL->Add(outfilename,NumHbonds);

  mprintf( "  HBOND: ");
  if (!hasDonorMask && !hasAcceptorMask)
    mprintf("Searching for Hbond donors/acceptors in region specified by %s\n",Mask.MaskString());
  else if (hasDonorMask && !hasAcceptorMask)
    mprintf("Donor mask is %s, acceptors will be searched for in region specified by %s\n",
            DonorMask.MaskString(),Mask.MaskString());
  else if (hasAcceptorMask && !hasDonorMask)
    mprintf("Acceptor mask is %s, donors will be searched for in a region specified by %s\n",
            AcceptorMask.MaskString(),Mask.MaskString());
  else
    mprintf("Donor mask is %s, Acceptor mask is %s\n",
            DonorMask.MaskString(),AcceptorMask.MaskString());
  mprintf( "         Distance cutoff = %8.3lf, Angle Cutoff = %8.3lf\n",dcut,acut*RADDEG);
  if (outfilename!=NULL) 
    mprintf( "         Dumping # Hbond v time results to %s\n", outfilename);
  if (avgout!=NULL)
    mprintf( "         Dumping Hbond avgs to %s\n",avgout);

  return 0;
}

// Hbond::SearchAcceptor()
/** Search for hbond acceptors X in the region specified by amask.
  * If Auto is true select acceptors based on the rule that "Hydrogen 
  * bonds are FON"
  */
void Hbond::SearchAcceptor(AtomMask *amask, bool Auto) {
  int atom;
  bool isAcceptor;
  // Set up acceptors: F, O, N
  // NOTE: Attempt to determine electronegative carbons?
  for (int selected=0; selected < amask->Nselected; selected++) {
    atom = amask->Selected[selected];
    isAcceptor=true;
    // If auto searching, only consider acceptor atoms as F, O, N
    if (Auto) {
      isAcceptor=false;
      if (currentParm->AtomElementIs(atom,FLUORINE) ||
          currentParm->AtomElementIs(atom,OXYGEN)   ||
          currentParm->AtomElementIs(atom,NITROGEN)   )
        isAcceptor=true;
    }
    if (isAcceptor)
      Acceptor.push_back(atom);
  }
}

// Hbond::SearchDonor()
/** Search for hydrogen bond donors X-H in the region specified by dmask.
  * If Auto is true select donors based on the rule that "Hydrogen bonds 
  * are FON"
  */
void Hbond::SearchDonor(AtomMask *dmask, bool Auto) {
  int donoratom;//, atom1, atom2;
  bool isDonor;
  std::vector<int> HatomList;
  // Set up donors: F-H, O-H, N-H
  for (int selected=0; selected < dmask->Nselected; selected++) {
    donoratom = dmask->Selected[selected];
    // If this is already an H atom continue
    if (currentParm->AtomElementIs(donoratom,HYDROGEN)) continue;
    isDonor = true;
    // If auto searching, only consider donor atoms as F, O, N
    if (Auto) {
      isDonor=false;
      if (currentParm->AtomElementIs(donoratom,FLUORINE) ||
          currentParm->AtomElementIs(donoratom,OXYGEN)   ||
          currentParm->AtomElementIs(donoratom,NITROGEN)   )
        isDonor=true;
    }
    if (isDonor) {
      // Get list of hydrogen atoms bonded to this atom
      if ( currentParm->GetBondedHatoms(donoratom, HatomList) ) {
        for (std::vector<int>::iterator h_atom = HatomList.begin();
                                        h_atom != HatomList.end();
                                        h_atom++)
        {
          //mprintf("BOND TO H: %i@%s -- %i@%s\n",donoratom+1,currentParm->AtomName(donoratom),
          //        *h_atom+1,currentParm->AtomName(*h_atom));
          Donor.push_back(donoratom);
          Donor.push_back(*h_atom);
        }
      }
    } // END atom is potential donor
  } // END loop over selected atoms
}

// Hbond::setup()
/** Search for hbond donors and acceptors. 
  */
int Hbond::setup() {
  int atom, a2;

  // Set up bond information for parm
  if (currentParm->SetupBondInfo()) return 1;

  // Set up mask
  if (!hasDonorMask && !hasAcceptorMask) {
    if ( currentParm->SetupIntegerMask( Mask, activeReference) ) return 1;
    if ( Mask.None() ) {
      mprintf("    Error: Hbond::setup: Mask has no atoms.\n");
      return 1;
    }
  }
  // Set up donor mask
  if (hasDonorMask) {
    if ( currentParm->SetupIntegerMask( DonorMask, activeReference) ) return 1;
    if (DonorMask.None()) {
      mprintf("    Error: Hbond: DonorMask has no atoms.\n");
      return 1;
    }
  }
  // Set up acceptor mask
  if (hasAcceptorMask) {
    if ( currentParm->SetupIntegerMask( AcceptorMask, activeReference) ) return 1;
    if (AcceptorMask.None()) {
      mprintf("    Error: Hbond: AcceptorMask has no atoms.\n");
      return 1;
    }
  }

  // Four cases:
  // 1) DonorMask and AcceptorMask NULL: donors and acceptors automatically searched for.
  if (!hasDonorMask && !hasAcceptorMask) {
    SearchAcceptor(&Mask,true);
    SearchDonor(&Mask,true);
  
  // 2) DonorMask only: acceptors automatically searched for in Mask
  } else if (hasDonorMask && !hasAcceptorMask) {
    SearchAcceptor(&Mask,true);
    SearchDonor(&DonorMask, false);

  // 3) AcceptorMask only: donors automatically searched for in Mask
  } else if (!hasDonorMask && hasAcceptorMask) {
    SearchAcceptor(&AcceptorMask, false);
    SearchDonor(&Mask,true);

  // 4) Both DonorMask and AcceptorMask: No automatic search.
  } else {
    SearchAcceptor(&AcceptorMask, false);
    SearchDonor(&DonorMask, false);
  }

  // Print acceptor/donor information
  mprintf("      HBOND: Set up %i acceptors:\n",(int)Acceptor.size());
  if (debug>0) {
    for (accept = Acceptor.begin(); accept!=Acceptor.end(); accept++)
      mprintf("        %8i: %4s\n",*accept+1,currentParm->AtomName(*accept));
  }
  mprintf("      HBOND: Set up %i donors:\n",((int)Donor.size())/2);
  if (debug>0) {
    for (donor = Donor.begin(); donor!=Donor.end(); donor++) {
      atom = (*donor);
      donor++;
      a2   = (*donor);
      mprintf("        %8i:%4s - %8i:%4s\n",atom+1,currentParm->AtomName(atom),
              a2+1,currentParm->AtomName(a2)); 
    } 
  }

  return 0;
}

// Hbond::action()
/** Calculate distance between all donors and acceptors. Store Hbond info.
  */    
int Hbond::action() {
  // accept ... H-D
  int D, H, Nhb, numHB;
  double dist, dist2, angle;//, ucell[9], recip[9];
  std::map<int,HbondType>::iterator it;
  HbondType HB;

  Nhb = 0; numHB=0;
  for (donor = Donor.begin(); donor!=Donor.end(); donor++) {
    D = (*donor);
    donor++;
    H = (*donor);
    for (accept = Acceptor.begin(); accept!=Acceptor.end(); accept++, Nhb++) {
      if (*accept == D) continue;
      dist2 = currentFrame->DIST2(*accept, D);
      //dist2 = currentFrame->DIST2(*accept, D, (int)P->boxType, ucell, recip);
      if (dist2 > dcut2) continue;
      angle = currentFrame->ANGLE(*accept, H, D);
      if (angle < acut) continue;
//      mprintf( "HBOND[%i]: %i:%s ... %i:%s-%i:%s Dist=%lf Angle=%lf\n", 
//              Nhb, *accept, P->names[*accept],
//              H, P->names[H], D, P->names[D], dist, angle);
      numHB++;
      dist = sqrt(dist2);
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
  NumHbonds->Add(frameNum, &numHB);
//  mprintf("HBOND: Scanned %i hbonds.\n",Nhb);
  Nframes++;

  return 0;
}

// Hbond::print()
/** Print average occupancies over all frames for all detected Hbonds
  */
void Hbond::print() {
  std::map<int,HbondType>::iterator it;
  std::list<HbondType> HbondList;
  std::list<HbondType>::iterator hbond;
  double avg, dist, angle;
  char Aname[32], Hname[32], Dname[32];
  DataFile *hbavgFile;

  // If avgout is NULL no averaging.
  if (avgout==NULL) return;

  // Set up data set list for all avg-related data.
  HBavg = new DataSetList(); 
  DFL->Add(avgout, HBavg->Add(STRING, (char*)"Acceptor", "Acceptor"));
  DFL->Add(avgout, HBavg->Add(STRING, (char*)"DonorH", "DonorH"));
  DFL->Add(avgout, HBavg->Add(STRING, (char*)"Donor", "Donor"));
  DFL->Add(avgout, HBavg->Add(INT, (char*)"Frames", "Frames"));
  DFL->Add(avgout, HBavg->Add(DOUBLE, (char*)"Frac", "Frac"));
  DFL->Add(avgout, HBavg->Add(DOUBLE, (char*)"AvgDist", "AvgDist"));
  hbavgFile = DFL->Add(avgout, HBavg->Add(DOUBLE, (char*)"AvgAng", "AvgAng"));

  //if (OutFile.SetupFile(avgout, WRITE, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug)) return;
  //OutFile.OpenFile();

  //OutFile.IO->Printf("HBONDS:\n");
  //OutFile.IO->Printf("  %-15s %-15s %-15s %6s %6s %8s %8s\n",
  //               "Acceptor","DonorH","Donor","Frames","Frac","AvgDist","AvgAng");
  for (it = HbondMap.begin(); it!=HbondMap.end(); it++) {
//    mprintf( "      %8i:%s ... %8i:%s-%8i:%s %8i Frames.\n",
//              (*it).second.A, P->names[(*it).second.A],
//              (*it).second.H, P->names[(*it).second.H],
//              (*it).second.D, P->names[(*it).second.D],
//              (*it).second.Frames);
    HbondList.push_back( (*it).second );
  }

  HbondList.sort( hbond_cmp() );
  int hbondnum=0;
  for ( hbond = HbondList.begin(); hbond!=HbondList.end(); hbond++ ) {
    avg = (double) (*hbond).Frames;
    avg = avg / ((double) Nframes);
    dist = (double) (*hbond).dist;
    dist = dist / ((double) (*hbond).Frames);
    angle = (double) (*hbond).angle;
    angle = angle / ((double) (*hbond).Frames);
    angle *= RADDEG;

    currentParm->ResAtomName(Aname, (*hbond).A);
    currentParm->ResAtomName(Hname, (*hbond).H);
    currentParm->ResAtomName(Dname, (*hbond).D);

    HBavg->AddData(hbondnum, Aname, 0);
    HBavg->AddData(hbondnum, Hname, 1);
    HBavg->AddData(hbondnum, Dname, 2);
    HBavg->AddData(hbondnum, &((*hbond).Frames), 3);
    HBavg->AddData(hbondnum, &avg, 4);
    HBavg->AddData(hbondnum, &dist, 5);
    HBavg->AddData(hbondnum, &angle, 6);
    hbondnum++;
    //OutFile.IO->Printf("%-15s %-15s %-15s %6i %6.2lf %8.3lf %8.3lf\n",
    //                   Aname,Hname,Dname, (*hbond).Frames,avg,dist,angle);
  }
  //OutFile.CloseFile();
  hbavgFile->SetNoXcol();
}
