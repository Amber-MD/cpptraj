// FrameList
#include "FrameList.h"
#include "Trajin_Single.h"
#include "CpptrajStdio.h"

const ReferenceFrame FrameList::ErrorFrame_(0,0,"",-1);

// CONSTRUCTOR
FrameList::FrameList() : 
  refFrameNum_(0)
{}

// DESTRUCTOR
FrameList::~FrameList() {
  Clear();
}

/** Clear the FrameList. */
void FrameList::Clear() {
  for (std::vector<ReferenceFrame>::iterator ref = frames_.begin();
                                             ref != frames_.end(); ++ref)
    delete (*ref).Coord();
  frames_.clear();
  for (std::vector<Topology*>::iterator parm = StrippedRefParms_.begin();
                                        parm != StrippedRefParms_.end(); parm++)
    delete *parm;
  StrippedRefParms_.clear();
  refFrameNum_ = 0;
  FileList::Clear();
}

// -----------------------------------------------------------------------------
// FrameList::ActiveReference()
/** Return the address of the frame pointed to by refFrameNum_.
  */
Frame* FrameList::ActiveReference() {
  if (frames_.empty()) return 0;
  return frames_[refFrameNum_].Coord();
}

// FrameList::SetActiveRef()
/** Set the given frame list number as the active reference.
  */
void FrameList::SetActiveRef(int numIn) {
  if (numIn < 0 || numIn >= (int)frames_.size()) {
    mprintf("Warning: FrameList::SetActiveRef: Ref # %i out of bounds.\n",numIn);
    return;
  }
  refFrameNum_ = numIn;
}
// -----------------------------------------------------------------------------
/** Add Frame from the given trajectory file to the FrameList. Store trajectory 
  * name and frame number that this frame came from in frameNames and frameNums 
  * respectively. Store the associated parm in FrameParm. 
  */
int FrameList::AddReference(ArgList& argIn, TopologyList& topListIn) {
  Trajin_Single traj;

  traj.SetDebug(debug_);
  // Get associated parmtop
  Topology* parmIn = topListIn.GetParm( argIn );
  // Check if we want to obtain the average structure
  bool average = argIn.hasKey("average");

  // Set up trajectory - false = do not modify box info
  if ( traj.SetupTrajRead( argIn.GetStringNext(), &argIn, parmIn, false ) ) {
    mprinterr("Error: reference: Could not set up trajectory.\n");
    return 1;
  }

  // Check for mask expression
  // TODO: This is done after SetupTrajRead because forward slash is a valid
  //       mask operand for element, so getNextMask will unfortunately pick up 
  //       filenames as well as masks.
  std::string maskexpr = argIn.GetMaskNext();

  // Check for tag - done after SetupTrajRead so traj can process args
  std::string reftag = argIn.getNextTag();

  // Check and warn if filename/reftag currently in use
  if (FindName( traj.TrajFilename().Full() )!=-1) {
    mprintf("Warning: Reference with name %s already exists!\n",traj.TrajFilename().full());
    //return 1;
  }
  if (FindName(reftag)!=-1) {
    mprintf("Warning: Reference with tag %s already exists!\n",reftag.c_str());
    //return 1;
  }

  // If not obtaining average structure, tell trajectory to only process
  // the start frame.
  if (!average)
    traj.SingleFrame();

  // Check number of frames to be read
  int trajFrames = traj.TotalReadFrames();
  if (trajFrames < 1) {
    mprinterr("Error: No frames could be read for reference %s\n", traj.TrajFilename().base());
    return 1;
  }
  // Start trajectory read
  if ( traj.BeginTraj(false) ) {
    mprinterr("Error: Could not open reference %s\n.", traj.TrajFilename().base());
    return 1;
  }
  traj.PrintInfo(1);
  // Set up input frame
  // NOTE: If ever need ref velocity change this alloc
  Frame *CurrentFrame = new Frame( parmIn->Atoms() );
  // If averaging requested, loop over specified frames and avg coords.
  if (average) {
    mprintf("    Averaging over %i frames.\n",trajFrames);
    Frame *AvgFrame = new Frame( parmIn->Atoms() );
    AvgFrame->ZeroCoords();
    while ( traj.GetNextFrame( *CurrentFrame ) ) {
      *AvgFrame += *CurrentFrame;
    }
    AvgFrame->Divide( (double)trajFrames );
    delete CurrentFrame;
    CurrentFrame = AvgFrame;
  } else {
    // Not averaging, get the one frame from traj
    traj.GetNextFrame( *CurrentFrame );
  }
  traj.EndTraj();

  // If a mask expression was specified, strip to match the expression.
  Topology *CurrentParm = parmIn;
  if (!maskexpr.empty()) {
    AtomMask stripMask( maskexpr );
    mprintf("    reference: Keeping atoms in mask [%s]\n",stripMask.MaskString());
    if (parmIn->SetupIntegerMask(stripMask, *CurrentFrame)) {
      delete CurrentFrame;
      return 1;
    }
    if (stripMask.None()) {
      mprinterr("Error: No atoms kept for reference.\n");
      delete CurrentFrame;
      return 1;
    }
    // Create new stripped frame
    Frame *strippedRefFrame = new Frame( *CurrentFrame, stripMask );
    mprintf("\tKept %i atoms.\n", strippedRefFrame->Natom());
    // Create new stripped parm
    Topology *strippedRefParm = CurrentParm->modifyStateByMask( stripMask );
    if (strippedRefParm==0) {
      mprinterr("Error: could not strip reference.\n");
      return 1;
    }
    strippedRefParm->Summary();
    // Store the new stripped parm in this class so it can be freed later
    StrippedRefParms_.push_back( strippedRefParm );
    delete CurrentFrame;
    CurrentFrame = strippedRefFrame;
    // No need to free CurrentParm since it exists in the parm file list.
    CurrentParm = strippedRefParm;
  }

  frames_.push_back( ReferenceFrame( CurrentFrame, CurrentParm, 
                                     traj.TrajFilename().Base(), traj.Start() ) );
  AddNameWithTag( traj.TrajFilename(), reftag );
  // If the top currently has no reference coords, set them now
  if (CurrentParm->NoRefCoords()) CurrentParm->SetReferenceCoords( CurrentFrame ); 
  return 0;
}

// FrameList::GetFrameFromArgs()
/** \return ReferenceFrame based on args in argIn.
  * The keywords in order of precedence are:
  *   - 'ref <name>'  : Get reference frame by name/tag.
  *   - 'reference'   : First reference frame in list.
  *   - 'refindex <#>': Reference frame at position.
  */
ReferenceFrame FrameList::GetFrameFromArgs(ArgList& argIn) const {
  int refindex;
  // By name/tag
  std::string refname = argIn.GetStringKey("ref");
  if (!refname.empty()) {
    refindex = FindName( refname );
    if (refindex == -1) {
      mprinterr("Error: Could not get reference with name %s\n", refname.c_str());
      return ErrorFrame_;
    }
    return frames_[refindex];
  }
  // First defined reference
  if (argIn.hasKey("reference")) {
    if (frames_.empty()) {
      mprinterr("Error: No reference frames defined.\n");
      return ErrorFrame_;
    }
    return frames_[0];
  }
  // By index
  refindex = argIn.getKeyInt("refindex", -1);
  if (refindex != -1) {
    if (refindex < 0 || refindex >= (int)frames_.size()) {
      mprinterr("Error: reference index %i is out of bounds.\n", refindex);
      return ErrorFrame_;
    }
    return frames_[refindex];
  }
  // No frame desired, return empty.
  return ReferenceFrame();
}

// FrameList::GetFrameByName()
ReferenceFrame FrameList::GetFrameByName(std::string const& refName) const{
  int refIndex = FindName( refName );
  if (refIndex < 0)
    return ReferenceFrame();
  return frames_[refIndex];
}

// FrameList::ReplaceFrame()
/** Replace the frame/parm at the given position with the given frame/parm.
  * The old frame is deleted. 
  */
int FrameList::ReplaceFrame(ReferenceFrame const& refIn, Frame *newFrame, Topology *newParm) {
  if (newFrame==0 || newParm==0) return 1;
  for (std::vector<ReferenceFrame>::iterator ref = frames_.begin(); 
                                             ref != frames_.end(); ++ref) 
  {
    if ( *ref == refIn ) {
      delete (*ref).Coord();
      (*ref).SetRef( newFrame, newParm );
      return 0;
    }
  }
  return 1;
}

// FrameList::List()
/** Print a list of trajectory names that frames have been taken from.
  */
void FrameList::List() const {
  if (frames_.empty()) {
    mprintf("  No frames defined.\n");
    return;
  }
  if (HasNames()) {
    mprintf("  The following %zu frames have been defined:\n",frames_.size());
    for (int fn=0; fn < (int)frames_.size(); fn++) { 
      if (!Tag(fn).empty())
        mprintf("    %i: %s frame %i\n", fn, Tag(fn).c_str(),
                frames_[fn].Num()+1);
      else
        mprintf("    %i: %s frame %i\n",fn,Name(fn).base(),
                frames_[fn].Num()+1);
    }
  } else {
    mprintf("  %zu frames have been defined.\n",frames_.size());
  }
  mprintf("\tActive reference frame for masks is %i\n",refFrameNum_);
}
