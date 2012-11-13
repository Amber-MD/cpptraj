// FrameList
#include "FrameList.h"
#include "Trajin_Single.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
FrameList::FrameList() : 
  refFrameNum_(0)
{}

// DESTRUCTOR
FrameList::~FrameList() {
  for (std::vector<Frame*>::iterator frame = frames_.begin();
                                     frame != frames_.end(); frame++)
    delete *frame;
  for (std::vector<Topology*>::iterator parm = StrippedRefParms_.begin();
                                        parm != StrippedRefParms_.end(); parm++)
    delete *parm;
}

// -----------------------------------------------------------------------------
// FrameList::ActiveReference()
/** Return the address of the frame pointed to by refFrameNum_.
  */
Frame *FrameList::ActiveReference() {
  if (frames_.empty()) return NULL;
  return frames_[refFrameNum_];
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

  // Set up trajectory
  if ( traj.SetupTrajRead( argIn.GetStringNext(), &argIn, parmIn ) ) {
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
  if (FindName( traj.TrajName().Full() )!=-1) {
    mprintf("Warning: Reference with name %s already exists!\n",traj.FullTrajStr());
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
    mprinterr("Error: No frames could be read for reference %s\n", traj.BaseTrajStr());
    return 1;
  }
  // Start trajectory read
  if ( traj.BeginTraj(false) ) {
    mprinterr("Error: Could not open reference %s\n.", traj.BaseTrajStr());
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
    if (strippedRefParm==NULL) {
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

  frames_.push_back( CurrentFrame );
  AddNameWithTag( traj.TrajName().Full(), traj.TrajName().Base(), reftag );
  nums_.push_back( traj.Start() );
  parms_.push_back( CurrentParm );
  return 0;
}

// FrameList::AddFrame()
/** Add given Frame to the FrameList. Store the associated parm in FrameParm.
  */
int FrameList::AddFrame(Frame *F, Topology *P) {
  if (F==NULL || P==NULL) return 1;
  frames_.push_back(F);
  parms_.push_back(P);
  return 0;
}

// FrameList::GetFrameParm()
/** Given index of frame, return parm in FrameParm
  */
Topology *FrameList::GetFrameParm(int idx) {
  if (idx < 0 || idx >= (int)parms_.size()) return NULL;
  return parms_[idx];
}

// FrameList::GetFrame()
/** Return the frame in the frame list specified by index.
  */
Frame *FrameList::GetFrame(int idx) {
  if (idx<0 || idx>=(int)frames_.size()) return NULL;
  return frames_[idx];
}

// FrameList::ReplaceFrame()
/** Replace the frame/parm at the given position with the given frame/parm.
  * The old frame is deleted. 
  */
int FrameList::ReplaceFrame(int idx, Frame *newFrame, Topology *newParm) {
  if (newFrame==NULL || newParm==NULL) return 1;
  if (idx<0 || idx>=(int)frames_.size()) return 1;
  delete frames_[idx];
  frames_[idx] = newFrame;
  parms_[idx] = newParm;
  return 0;
}

// FrameList::List()
/** Print a list of trajectory names that frames have been taken from.
  */
void FrameList::List() {
  if (frames_.empty()) {
    mprintf("  No frames defined.\n");
    return;
  }
  if (HasNames()) {
    mprintf("  The following %zu frames have been defined:\n",frames_.size());
    for (int fn=0; fn < (int)frames_.size(); fn++) { 
      if (!Tag(fn).empty())
        mprintf("    %i: %s frame %i\n",fn,Tag(fn).c_str(),
                nums_[fn]+OUTPUTFRAMESHIFT);
      else
        mprintf("    %i: %s frame %i\n",fn,Name(fn).c_str(),
                nums_[fn]+OUTPUTFRAMESHIFT);
    }
  } else {
    mprintf("  %zu frames have been defined.\n",frames_.size());
  }
  mprintf("\tActive reference frame for masks is %i\n",refFrameNum_);
}

// FrameList::FrameName()
/** Return name of given frame.
  */
const char *FrameList::FrameName(int idx) {
  if (idx<0 || idx>=(int)frames_.size()) return NULL;
  if (Name(idx).empty()) return NULL;
  return Name(idx).c_str();
}

