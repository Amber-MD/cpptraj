#include "Action_Contacts.h"
#include "CpptrajStdio.h"
// TODO: Separate byResidue stuff

// CONSTRUCTOR
Action_Contacts::Action_Contacts() :
  byResidue_(false),
  distance_(49.0), // 7.0^2
  first_(false)
{
  // Do not allow reference traj
  //ForbidRefMode( REFTRAJ );
}

// DESTRUCTOR
Action_Contacts::~Action_Contacts() {
  outfile_.CloseFile();
  if (byResidue_)
    outfile2_.CloseFile();
}

/** Set up native contacts based on reference structure. */
int Action_Contacts::SetupContacts(Frame* refframe, Topology* refparm) {
  // Determine which pairs of atoms satisfy cutoff, build contact list
  AtomMask::const_iterator atom1end = Mask_.end() - 1;
  for (AtomMask::const_iterator atom1 = Mask_.begin();
                                atom1 != atom1end; ++atom1)
  {
    for (AtomMask::const_iterator atom2 = atom1 + 1;
                                  atom2 != Mask_.end(); ++atom2)
    {
      double d2 = refframe->DIST2(*atom1, *atom2);
      if (d2 < distance_)
        nativecontacts_.insert( contactType(*atom1, *atom2) );
    }
  }
  // DEBUG - Print contacts
  mprintf("\tSetup %zu native contacts:\n", nativecontacts_.size());
  for (contactListType::iterator contact = nativecontacts_.begin();
                                 contact != nativecontacts_.end(); ++contact)
  {
    int a1 = (*contact).first;
    int a2 = (*contact).second;
    mprintf("\t\tAtom %i[%s] to %i[%s]\n", a1+1, (*refparm)[a1].c_str(),
            a2+1, (*refparm)[a2].c_str());
  }
  return 0;
}

/** USAGE:
  *
  *  contacts [ first | reference | ref <ref> | refindex <#> ] [byresidue]
  *      [out <filename>] [time <interval>] [distance <cutoff>] [<mask>]
  *
  *  action argument usage:
  *
  *  byresidue: calculate number of contacts for every specified atom and save result per residue
  */
int Action_Contacts::init() {

  byResidue_ = actionArgs.hasKey("byresidue");
  double dist = actionArgs.getKeyDouble("distance", 7.0);
  // Square the cutoff
  distance_ = dist * dist;
  first_ = actionArgs.hasKey("first");
  char* referenceName = actionArgs.getKeyString("ref",0);
  int refindex = actionArgs.getKeyInt("refindex",-1);
  // For compatibility with ptraj, keyword 'reference' == 'refindex 0'
  if (actionArgs.hasKey("reference")) refindex = 0;
  char* outfilename = actionArgs.getKeyString("out",0); // NOTE: Should be NULL?
  if (outfile_.SetupWrite(outfilename, debug))
    return 1;
  if (outfile_.OpenFile())
    return 1;
  if (byResidue_) {
    if (outfilename==0) {
      mprinterr("Error: Contacts 'byresidue' requires output filename.\n");
      return 1;
    }
    std::string file2name( outfilename );
    file2name += ".native";
    if (outfile2_.SetupWrite( file2name.c_str(), debug ))
      return 1;
    if (outfile2_.OpenFile())
      return 1;
  }

  // Get Mask
  char* mask0 = actionArgs.getNextMask();
  if (mask0==0 && byResidue_)
    Mask_.SetMaskString("@CA");
  else
    Mask_.SetMaskString( mask0 );
  
  // Initialize reference. If no reference mask is given mask0 will be used.
  // First arg 'nofit' set to true, no fitting with contacts. Allows last arg
  // 'RefTrans' to be NULL.
  //if (RefInit(true, false, Mask_.MaskString(), actionArgs, FL, PFL, NULL)!=0)
  //  return 1;
  if (!first_ && referenceName==0 && refindex==-1) {
    mprintf("\tNo reference structure specified. Defaulting to first.\n");
    first_ = true;
  }

  if (!first_) {
    // Attempt to get the reference index by name/tag
    if (referenceName != 0)
      refindex = FL->FindName(referenceName);
    // Get reference frame by index
    // TODO: Convert FrameList to return frame reference?
    Frame* RefFrame = FL->GetFrame(refindex);
    if (RefFrame == 0) {
      mprinterr("Error: Could not get reference index %i\n",refindex);
      return 1;
    }
    // Get reference parm for frame
    Topology* RefParm = FL->GetFrameParm(refindex);
    if (RefParm==0) {
      mprinterr("Error: Could not get parm for frame %s\n", FL->FrameName(refindex));
      return 1;
    }
    // Set up atom mask for reference frame
    if (RefParm->SetupIntegerMask(Mask_, *RefFrame)) return 1;
    // Set up reference contacts 
    SetupContacts(RefFrame, RefParm);
  }

  // Output file header - only if not byresidue
  if (!byResidue_) {
    outfile_.Printf("#time\tContacts\tnative Contacts ");
    if (!first_)
      outfile_.Printf("(number of natives: %zu)", nativecontacts_.size());
    outfile_.Printf("\n");
  }

  mprintf("    CONTACTS: [%s] Calculating current contacts and comparing results to",
          Mask_.MaskString());
  if (first_)
    mprintf(" first frame.\n");
  else
    mprintf(" reference structure.\n");
  mprintf("                   Distance cutoff is %lf angstroms.\n", dist);
  if (outfilename==0)
    mprintf("              Results will be written to stdout");
  else
    mprintf("              Writing results to %s", outfilename);
  mprintf("\n");
  if (byResidue_)
    mprintf("              Results are output on a per-residue basis.\n");

  return 0;
}

int Action_Contacts::setup() {
  //if (first_) 
  //  RefParm_ = currentParm;
  // Set up atom mask 
  if (currentParm->SetupIntegerMask(Mask_)) return 1;

  // Determine which residues are active based on the mask
  activeResidues_.clear();
  for (AtomMask::const_iterator atom = Mask_.begin();
                                atom != Mask_.end(); ++atom)
  {
    int resnum = (*currentParm)[*atom].ResNum();
    activeResidues_.insert( resnum );
  }

  // byresidue header - only on first time through
  if (residueContacts_.empty() && byResidue_) {
    outfile_.Printf("#time");
    outfile2_.Printf("#time");
    for (std::set<int>::iterator res = activeResidues_.begin();
                                 res != activeResidues_.end(); ++res)
    {
      outfile_.Printf("\tresidue %i", *res);
      outfile2_.Printf("\tresidue %i", *res);
    }
    outfile_.Printf("\tTotal\n");
    outfile2_.Printf("\tTotal\n");
  }

  // Reserve space for residue contact counts
  residueContacts_.reserve( currentParm->Nres() );
  residueNative_.reserve( currentParm->Nres() );


  return 0;
}

int Action_Contacts::action() {
  if (first_) {
    SetupContacts( currentFrame, currentParm );
    first_ = false;
  }

  // Determine how many contacts and how many native contacts are present.
  //contactListType::iterator native = nativecontacts_.begin();
  residueContacts_.assign( currentParm->Nres(), 0 );
  residueNative_.assign( currentParm->Nres(), 0 );

  // Determine which pairs of atoms satisfy cutoff
  int numcontacts = 0;
  int numnative = 0;
  AtomMask::const_iterator atom1end = Mask_.end() - 1;
  for (AtomMask::const_iterator atom1 = Mask_.begin();
                                atom1 != atom1end; ++atom1)
  {
    int res1 = (*currentParm)[*atom1].ResNum();
    for (AtomMask::const_iterator atom2 = atom1 + 1;
                                  atom2 != Mask_.end(); ++atom2)
    {
      double d2 = currentFrame->DIST2(*atom1, *atom2);
      // Contact?
      if (d2 < distance_) {
        ++numcontacts;
        int res2 = (*currentParm)[*atom2].ResNum();
        mprintf("CONTACT: %i res %i to %i res %i [%i]",*atom1,res1,*atom2,res2,numcontacts);
        ++residueContacts_[res1];
        ++residueContacts_[res2];
        // Is this a native contact?
        contactListType::iterator nativebegin = nativecontacts_.lower_bound( *atom1 );
        if (nativebegin != nativecontacts_.end()) {
          contactListType::iterator nativeend = nativecontacts_.upper_bound( *atom1 );
          for (contactListType::iterator native = nativebegin; native!=nativeend; ++native) {
            if ( *atom2 == (*native).second ) {
              ++numnative;
              mprintf(" NATIVE [%i]",numnative);
              ++residueNative_[res1];
              ++residueNative_[res2];
              break;
            }
          }
        }
        mprintf("\n");
      }
/*
      if (native != nativecontacts_.end()) {
        if ( *atom1 == (*native).first && *atom2 == (*native).second ) {
          if (d2 < distance_)
            ++numnative;
          ++native;
        }
      }
*/
    }
  }

  // The total # of contacts is multiplied by 2 since each contact
  // is bidirectional, i.e. atom1->atom2 and atom2->atom1 each
  // count as a separate contact.
  if (!byResidue_) {
    outfile_.Printf("%10.2lf\t%i\t%i\n", (double)frameNum+OUTPUTFRAMESHIFT,
                    numcontacts*2, numnative*2);
  } else {
    outfile_.Printf("%10.2lf", (double)frameNum+OUTPUTFRAMESHIFT);
    outfile2_.Printf("%10.2lf", (double)frameNum+OUTPUTFRAMESHIFT);
    for (std::set<int>::iterator res = activeResidues_.begin();
                                 res != activeResidues_.end(); ++res)
    {
      outfile_.Printf("\t%i", residueContacts_[ *res ]);
      outfile2_.Printf("\t%i", residueNative_[ *res ]);
    }
    outfile_.Printf("\t%i\n", numcontacts*2);
    outfile2_.Printf("\t%i\n", numnative*2);
  }

  return 0;
}
     
