#include "Action_Contacts.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Contacts::Action_Contacts() :
  byResidue_(false),
  distance_(49.0), // 7.0^2
  RefParm_(0),
  first_(false)
{
  // Do not allow reference traj
  //ForbidRefMode( REFTRAJ );
}

// DESTRUCTOR
Action_Contacts::~Action_Contacts() {
  outfile_.CloseFile();
}

/** Set up native contacts based on reference structure. */
int Action_Contacts::SetupContacts(Frame *refframe) {
  if (RefParm_ == 0) return 1;
  // Set up atom mask for reference frame
  if (RefParm_->SetupIntegerMask(Mask_, *refframe)) return 1;
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
        contactlist_.push_back( std::pair<int,int>(*atom1, *atom2) );
    }
  }
  // DEBUG - Print contacts
  mprintf("\tSetup %zu contacts:\n", contactlist_.size());
  for (std::vector< std::pair<int,int> >::iterator contact = contactlist_.begin();
                                                   contact != contactlist_.end(); ++contact)
  {
    int a1 = (*contact).first;
    int a2 = (*contact).second;
    mprintf("\t\tAtom %i[%s] to %i[%s]\n", a1+1, (*RefParm_)[a1].c_str(),
            a2+1, (*RefParm_)[a2].c_str());
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
    RefParm_ = FL->GetFrameParm(refindex);
    if (RefParm_==0) {
      mprinterr("Error: Could not get parm for frame %s\n", FL->FrameName(refindex));
      return 1;
    }
    // Set up reference contacts 
    SetupContacts(RefFrame);
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
  if (first_) 
    RefParm_ = currentParm;

  return 0;
}

int Action_Contacts::action() {
  if (first_) {
    SetupContacts( currentFrame );
    first_ = false;
  }

  return 0;
}
     
