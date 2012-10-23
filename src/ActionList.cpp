// ActionList
#include "ActionList.h"
#include "CpptrajStdio.h"
// All action classes go here
// TODO: Re-enable actions
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_Angle.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_DSSP.h"
#include "Action_Center.h"
#include "Action_Hbond.h"
#include "Action_Image.h"
#include "Action_Surf.h"
#include "Action_Radgyr.h"
#include "Action_Mask.h"
#include "Action_Closest.h"
#include "Action_NAstruct.h"
#include "Action_Pucker.h"
#include "Action_Outtraj.h"
#include "Action_Rms2d.h"
#include "Action_Average.h"
#include "Action_Radial.h"
#include "Action_DistRmsd.h"
#include "Action_Jcoupling.h"
#include "Action_Clustering.h"
//#include "Action_Pairwise.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
//#include "Action_DihedralScan.h"
#include "Action_Rotdif.h"
#include "Action_RunningAvg.h"
#include "Action_RmsAvgCorr.h"
#include "Action_AtomicFluct.h"
#include "Action_Watershell.h"
#include "Action_AvgCoord.h"
#include "Action_Contacts.h"
#include "Action_Vector.h"
#include "Action_Principal.h"
#include "Action_Matrix.h"
#include "Action_Grid.h"
#include "Action_GridFreeEnergy.h"
#include "Action_Dipole.h"
#include "Action_Projection.h"
#include "Action_ClusterDihedral.h"
#include "Action_Unwrap.h"
#include "Action_Diffusion.h"
#include "Action_DNAionTracker.h"
#include "Action_Scale.h"
#include "Action_RandomizeIons.h"
#include "Action_AutoImage.h"
#include "Action_STFC_Diffusion.h"
#include "Action_AtomicCorr.h"
#include "Action_Bounds.h"
#include "Action_Rotate.h"
#include "Action_Translate.h"

const DispatchObject::Token ActionList::DispatchArray[] = {
  { DispatchObject::ACTION, "2drms", Action_Rms2d::Alloc, Action_Rms2d::Help, 0 },
  { DispatchObject::ACTION, "angle", Action_Angle::Alloc, Action_Angle::Help, 0 },
  { DispatchObject::ACTION, "atomiccorr", Action_AtomicCorr::Alloc, Action_AtomicCorr::Help, 0 },
  { DispatchObject::ACTION, "atomicfluct", Action_AtomicFluct::Alloc, Action_AtomicFluct::Help, 0 },
    { DispatchObject::ACTION, "atommap", Action_AtomMap::Alloc, Action_AtomMap::Help, 0 },
  { DispatchObject::ACTION, "autoimage", Action_AutoImage::Alloc, Action_AutoImage::Help, 0 },
  { DispatchObject::ACTION, "average", Action_Average::Alloc, Action_Average::Help, 0 },
  { DispatchObject::ACTION, "avgcoord", Action_AvgCoord::Alloc, Action_AvgCoord::Help, 0 },
  { DispatchObject::ACTION, "bounds", Action_Bounds::Alloc, Action_Bounds::Help, 0 },
  { DispatchObject::ACTION, "center", Action_Center::Alloc, Action_Center::Help, 0 },
  { DispatchObject::ACTION, "check", Action_CheckStructure::Alloc, Action_CheckStructure::Help, 0 },
  { DispatchObject::ACTION, "checkstructure", Action_CheckStructure::Alloc, Action_CheckStructure::Help, 0 },
  { DispatchObject::ACTION, "closest", Action_Closest::Alloc, Action_Closest::Help, 0 },
  { DispatchObject::ACTION, "cluster", Action_Clustering::Alloc, Action_Clustering::Help, 0 },
  { DispatchObject::ACTION, "clusterdihedral", Action_ClusterDihedral::Alloc, Action_ClusterDihedral::Help, 0 },
  { DispatchObject::ACTION, "contacts", Action_Contacts::Alloc, Action_Contacts::Help, 0 },
  { DispatchObject::ACTION, "diffusion", Action_Diffusion::Alloc, Action_Diffusion::Help, 0 },
  { DispatchObject::ACTION, "dihedral", Action_Dihedral::Alloc, Action_Dihedral::Help, 0 },
//  { DispatchObject::ACTION, "dihedralscan", DihedralScan::Alloc, DihedralScan::Help, 0 },
  { DispatchObject::ACTION, "dipole", Action_Dipole::Alloc, Action_Dipole::Help, 0 },
  { DispatchObject::ACTION, "distance", Action_Distance::Alloc, Action_Distance::Help, 0 },
  { DispatchObject::ACTION, "dnaiontracker", Action_DNAionTracker::Alloc, Action_DNAionTracker::Help, 0 },
  { DispatchObject::ACTION, "drms", Action_DistRmsd::Alloc, Action_DistRmsd::Help, 0 },
  { DispatchObject::ACTION, "drmsd", Action_DistRmsd::Alloc, Action_DistRmsd::Help, 0 },
  { DispatchObject::ACTION, "gfe", Action_GridFreeEnergy::Alloc, Action_GridFreeEnergy::Help, 0 },
  { DispatchObject::ACTION, "grid", Action_Grid::Alloc, Action_Grid::Help, 0 },
  { DispatchObject::ACTION, "hbond", Action_Hbond::Alloc, Action_Hbond::Help, 0 },
  { DispatchObject::ACTION, "image", Action_Image::Alloc, Action_Image::Help, 0 },
  { DispatchObject::ACTION, "jcoupling", Action_Jcoupling::Alloc, Action_Jcoupling::Help, 0 },
  { DispatchObject::ACTION, "mask", Action_Mask::Alloc, Action_Mask::Help, 0 },
  { DispatchObject::ACTION, "matrix", Action_Matrix::Alloc, Action_Matrix::Help, 0 },
  { DispatchObject::ACTION, "molsurf", Action_Molsurf::Alloc, Action_Molsurf::Help, 0 },
  { DispatchObject::ACTION, "nastruct", Action_NAstruct::Alloc, Action_NAstruct::Help, 0 },
  { DispatchObject::ACTION, "outtraj", Action_Outtraj::Alloc, Action_Outtraj::Help, 0 },
//  { DispatchObject::ACTION, "pairwise", Pairwise::Alloc, Pairwise::Help, 0 },
  { DispatchObject::ACTION, "principal", Action_Principal::Alloc, Action_Principal::Help, 0 },
  { DispatchObject::ACTION, "projection", Action_Projection::Alloc, Action_Projection::Help, 0 },
  { DispatchObject::ACTION, "pucker", Action_Pucker::Alloc, Action_Pucker::Help, 0 },
  { DispatchObject::ACTION, "radgyr", Action_Radgyr::Alloc, Action_Radgyr::Help, 0 },
  { DispatchObject::ACTION, "radial", Action_Radial::Alloc, Action_Radial::Help, 0 },
  { DispatchObject::ACTION, "randomizeions", Action_RandomizeIons::Alloc, Action_RandomizeIons::Help, 0 },
  { DispatchObject::ACTION, "rms2d", Action_Rms2d::Alloc, Action_Rms2d::Help, 0 },
  { DispatchObject::ACTION, "rms", Action_Rmsd::Alloc, Action_Rmsd::Help, 0 },
  { DispatchObject::ACTION, "rmsd", Action_Rmsd::Alloc, Action_Rmsd::Help, 0 },
  { DispatchObject::ACTION, "rog", Action_Radgyr::Alloc, Action_Radgyr::Help, 0 },
  { DispatchObject::ACTION, "rotate", Action_Rotate::Alloc, Action_Rotate::Help, 0 },
  { DispatchObject::ACTION, "rotdif", Action_Rotdif::Alloc, Action_Rotdif::Help, 0 },
  { DispatchObject::ACTION, "runavg", Action_RunningAvg::Alloc, Action_RunningAvg::Help, 0 },
  { DispatchObject::ACTION, "runningaverage", Action_RunningAvg::Alloc, Action_RunningAvg::Help, 0 },
  { DispatchObject::ACTION, "scale", Action_Scale::Alloc, Action_Scale::Help, 0 },
  { DispatchObject::ACTION, "secstruct", Action_DSSP::Alloc, Action_DSSP::Help, 0 },
  { DispatchObject::ACTION, "stfcdiffusion", Action_STFC_Diffusion::Alloc, Action_STFC_Diffusion::Help, 0 },
  { DispatchObject::ACTION, "strip", Action_Strip::Alloc, Action_Strip::Help, 0 },
  { DispatchObject::ACTION, "surf", Action_Surf::Alloc, Action_Surf::Help, 0 },
  { DispatchObject::ACTION, "trans", Action_Translate::Alloc, Action_Translate::Help, 0 },
  { DispatchObject::ACTION, "unstrip", Action_Unstrip::Alloc, Action_Unstrip::Help, 0 },
  { DispatchObject::ACTION, "unwrap", Action_Unwrap::Alloc, Action_Unwrap::Help, 0 },
  { DispatchObject::ACTION, "vector", Action_Vector::Alloc, Action_Vector::Help, 0 },
  { DispatchObject::ACTION, "watershell", Action_Watershell::Alloc, Action_Watershell::Help, 0 },
  { DispatchObject::NONE,        0,                  0,                 0, 0 }
};

// CONSTRUCTOR
ActionList::ActionList() :
  debug_(0)
{}

// DESTRUCTOR
ActionList::~ActionList() {
    // No need to cast back to whatever action was allocd since Action destructor is virtual
    for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act)
      delete *act; 
}

// ActionList::SetDebug()
void ActionList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("ActionList DEBUG LEVEL SET TO %i\n",debug_);
}

int ActionList::AddAction(DispatchObject::DispatchAllocatorType Alloc, ArgList& argIn,
                          TopologyList* PFL, FrameList* FL, DataSetList* DSL,
                          DataFileList* DFL)
{
  Action* act = (Action*)Alloc();
  // Attempt to initialize action
  if ( act->Init( argIn, PFL, FL, DSL, DFL, debug_ ) != Action::OK ) {
    mprinterr("Error: Could not initialize action [%s]\n", argIn.Command());
    return 1;
  }
  actionlist_.push_back( act );
  actioncmd_.push_back( argIn.ArgLine() );
  actionstatus_.push_back( INIT );
  return 0;
}

// ActionList::Init()
/** Initialize non-parm-specific data for each action (like datasets). If an 
  * action cannot be initialized deactivate it. Also set action debug level.
  */
/*int ActionList::Init( DataSetList *DSL, FrameList *FL, DataFileList *DFL, 
                      TopologyList *PFL, bool exitOnError) 
{
  mprintf("\nACTIONS: Initializing %zu actions:\n",actionlist_.size());
  unsigned int actnum = 0;
  for (action_it act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    mprintf("  %u: [%s]\n", actnum++, (*act)->CmdLine());
    if ((*act)->Status() == Action::INACTIVE) {
      mprintf("Warning: Action %s is not active.\n", (*act)->ActionCommand());
    } else {
      if ( (*act)->Init( DSL, FL, DFL, PFL, debug_ ) ) {
        if (exitOnError) {
          mprinterr("Error: Init failed for [%s].\n", (*act)->CmdLine());
          return 1;
        } else {
          mprintf("Warning: Init failed for [%s]: DEACTIVATING\n",
                  (*act)->CmdLine());
        }
      }
    }
    mprintf("\n");
  }

  return 0;
}*/

// ActionList::SetupActions()
/** Attempt to set up all actions in the action list with the given parm
  * If an action cannot be set up skip it.
  */
int ActionList::SetupActions(Topology **ParmAddress) {
  Topology *OriginalParm = *ParmAddress;

  mprintf(".....................................................\n");
  mprintf("PARM [%s]: Setting up %zu actions.\n",(*ParmAddress)->c_str(),actionlist_.size());
  unsigned int actnum = 0;
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    // Only attempt to set up action if active 
    if (actionstatus_[actnum] != INACTIVE) {
      mprintf("  %u: [%s]\n", actnum, actioncmd_[actnum].c_str());
      actionstatus_[actnum] = SETUP;
      Action::RetType err = (*act)->Setup(*ParmAddress, ParmAddress);
      if (err == Action::ERR) {
        mprintf("Warning: Setup failed for [%s]: Skipping\n",
                actioncmd_[actnum].c_str());
        // Reset action status to INIT (pre-setup)
        actionstatus_[actnum] = INIT;
        //return 1;
      } else if (err == Action::USEORIGINALFRAME) {
        // Return value of 2 requests return to original parm
        *ParmAddress = OriginalParm;
      }
      //fprintf(stdout,"DEBUG: After Action %i Setup parmName is %s\n",act,P->parmName);
    }
    ++actnum;
  }
  //mprintf(".....................................................\n");

  return 0;
}

// ActionList::DoActions()
/** Perform actions in the action list on the given Frame. Skip actions not 
  * initialized or not setup. 
  * \param FrameAddress Memory address of the current frame.
  * \param frameNumIn The current frame number.
  * \return true if coordinate output should be suppressed.
  * \return false if coordinate output should be performed.
  */
bool ActionList::DoActions(Frame **FrameAddress, int frameNumIn) {
  Frame *OriginalFrame = *FrameAddress;

  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameNumIn);
  unsigned int actnum = 0;
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act) 
  {
    // Only do actions which were properly set up
    if (actionstatus_[actnum] == SETUP) { 
      // Perform action on frame
      Action::RetType err = (*act)->DoAction(frameNumIn, *FrameAddress, FrameAddress);
      // Check for action special conditions/errors
      if (err != Action::OK) {
        if (err == Action::USEORIGINALFRAME) {
          // Return value of 2 requests return to original frame
          *FrameAddress = OriginalFrame;
        } else if (err == Action::SUPPRESSCOORDOUTPUT) {
          // Skip the rest of the actions and suppress output. Necessary when
          // e.g. performing a running average over coords.
          return true;
        } else {
          // If here return type is ACTION_ERR.
          // Treat actions that fail as if they could not be set up
          mprintf("Warning: Action [%s] failed, frame %i.\n", actioncmd_[actnum].c_str(),
                frameNumIn);
          actionstatus_[actnum] = INIT;
        }
      }
    }
    ++actnum;
  }
  return false;
}

// ActionList::Print()
void ActionList::Print() {
  mprintf("\nACTION OUTPUT:\n");
  unsigned int actnum = 0;
  for (Aarray::iterator act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    // Skip deactivated actions
    if (actionstatus_[actnum++] != INACTIVE)
      (*act)->Print();
  }
}

void ActionList::List() {
  mprintf("ACTIONS:\n");
  for (unsigned int actnum = 0; actnum < actionlist_.size(); ++actnum)
    mprintf("  %u: [%s]\n", actnum, actioncmd_[actnum].c_str());   
}
