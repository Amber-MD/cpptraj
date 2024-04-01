#include "Action_MinMaxDist.h"
#include "CharMask.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Action_MinMaxDist::Action_MinMaxDist() :
  mode_(NO_MODE),
  distType_(NO_DIST)
{}

const char* Action_MinMaxDist::modeStr_[] = {
  "atoms",
  "residues",
  "molecules",
  0
};

const char* Action_MinMaxDist::distTypeStr_[] = {
  "minimum",
  "maximum",
  "minimum and maximum",
  0
};

// Action_MinMaxDist::Help()
void Action_MinMaxDist::Help() const {
  mprintf("\tmask1 <mask1> [mask2 <mask2>] [{byatom|byres|bymol}]\n"
          "\t[mindist] [maxdist] [bothdist] [noimage]\n"
          "  Record the min/max distance between atoms/residues/molecules.\n"
         );
}

// Action_MinMaxDist::Init()
Action::RetType Action_MinMaxDist::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  imageOpt_.InitImaging( !(actionArgs.hasKey("noimage")) );
  // Mask Keywords
  std::string mask1str = actionArgs.GetStringKey("mask1");
  if (mask1str.empty()) {
    mprinterr("Error: Must specify at least 'mask1'\n");
    return Action::ERR;
  }
  if (mask1_.SetMaskString( mask1str )) {
    mprinterr("Error: Could not set mask1 '%s'\n", mask1str.c_str());
    return Action::ERR;
  }
  std::string mask2str = actionArgs.GetStringKey("mask2");
  if (!mask2str.empty()) {
    if (mask2_.SetMaskString( mask2str )) {
      mprinterr("Error: Could not set mask2 '%s'\n", mask2str.c_str());
      return Action::ERR;
    }
  }
  // Mode args
  if (actionArgs.hasKey("byatom"))
    mode_ = BY_ATOM;
  else if (actionArgs.hasKey("byres"))
    mode_ = BY_RES;
  else if (actionArgs.hasKey("bymol"))
    mode_ = BY_MOL;
  // Distance calc type args
  bool calc_mindist = (actionArgs.hasKey("mindist") || (actionArgs[0] == "mindist"));
  bool calc_maxdist = (actionArgs.hasKey("maxdist") || (actionArgs[0] == "maxdist"));
  if (actionArgs.hasKey("bothdist"))
    distType_ = BOTH_DIST;
  // Default mode
  if (mode_ == NO_MODE)
    mode_ = BY_ATOM;
  // Default distance calc type
  if (distType_ == NO_DIST) {
    if (calc_mindist && calc_maxdist)
      distType_ = BOTH_DIST;
    else if (calc_mindist)
      distType_ = MIN_DIST;
    else if (calc_maxdist)
      distType_ = MAX_DIST;
    else {
      mprintf("Warning: No distance type specified and command name '%s' unrecognized. Using default.\n");
      distType_ = MIN_DIST;
    }
  }

  mprintf("    MINMAXDIST: Calculating %s distance for selected %s.\n",
          distTypeStr_[distType_], modeStr_[mode_]);
  mprintf("\tMask1: %s\n", mask1_.MaskString());
  if (mask2_.MaskStringSet()) {
    mprintf("\tMask2: %s\n", mask2_.MaskString());
  }
  if (imageOpt_.UseImage())
    mprintf("\tDistances will use minimum image convention if box info present.\n");
  else
    mprintf("\tDistances will not be imaged.\n");

  return Action::OK;
}

/** For DEBUG, print selected atoms in entity array., */
void Action_MinMaxDist::printEntities(Earray const& entities, AtomMask const& maskIn) const {
  mprintf("DEBUG: Selected %s in mask %s\n",
          modeStr_[mode_], maskIn.MaskString());
  for (Earray::const_iterator it = entities.begin(); it != entities.end(); ++it) {
    if (!it->name_.empty()) {
      mprintf("\t%s :", it->name_.c_str());
      for (AtomMask::const_iterator at = it->emask_.begin(); at != it->emask_.end(); ++at)
        mprintf(" %i", *at + 1);
      mprintf("\n");
    }
  }
}

/** Set up entities for atoms selected in given mask. */
int Action_MinMaxDist::setupEntityArray(Earray& entities, AtomMask const& maskIn,
                                        Topology const& topIn)
const
{
  entities.clear();
  CharMask cmask(maskIn.ConvertToCharMask(), maskIn.Nselected());

  if (mode_ == BY_RES) {
    entities.reserve(topIn.Nres());
    for (unsigned int ires = 0; ires < (unsigned int)topIn.Nres(); ires++) {
      Residue const& Res = topIn.Res(ires);
      if (ires >= entities.size()) {
        entities.resize(ires+1);
      }
      Entity& currentEntity = entities[ires];
      for (int at = Res.FirstAtom(); at != Res.LastAtom(); at++) {
        if (cmask.AtomInCharMask( at )) {
          if (currentEntity.name_.empty()) {
            currentEntity.name_ = Res.Name().Truncated();
            mprintf("NAME %s\n", currentEntity.name_.c_str());
          }
          currentEntity.emask_.AddSelectedAtom( at );
        }
      }
    } // END loop over residues
  } else {
    mprinterr("Internal Error: Action_MinMaxDist::setupEntityArray() Unhandled mode\n");
    return 1;
  }
  // DEBUG - print entities
  printEntities(entities, maskIn);
  return 0;
}
      

// Action_MinMaxDist::Setup()
Action::RetType Action_MinMaxDist::Setup(ActionSetup& setup)
{
  // Set up imaging info for this topology 
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );
  if (imageOpt_.ImagingEnabled())
    mprintf("\tDistance imaging on.\n");
  else
    mprintf("\tDistance imaging off.\n");
  // Set up masks
  if (setup.Top().SetupIntegerMask( mask1_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", mask1_.MaskString());
    return Action::OK;
  }
  mask1_.MaskInfo();
  if (mask1_.None()) {
    mprintf("Warning: Nothing selected by mask '%s'\n", mask1_.MaskString());
    return Action::SKIP;
  }
  if (mask2_.MaskStringSet()) {
    if (setup.Top().SetupIntegerMask( mask2_ )) {
      mprinterr("Error: Could not set up mask '%s'\n", mask2_.MaskString());
      return Action::OK;
    }
    mask2_.MaskInfo();
    if (mask2_.None()) {
      mprintf("Warning: Nothing selected by mask '%s'\n", mask2_.MaskString());
      return Action::SKIP;
    }
  }
  if (mode_ != BY_ATOM) {
    if (setupEntityArray(entities1_, mask1_, setup.Top())) {
      mprinterr("Error: Could not set up %s for mask '%s'\n",
                modeStr_[mode_], mask1_.MaskString());
      return Action::ERR;
    }
    if (mask2_.MaskStringSet()) {
      if (setupEntityArray(entities2_, mask2_, setup.Top())) {
        mprinterr("Error: Could not set up %s for mask '%s'\n",
                  modeStr_[mode_], mask2_.MaskString());
        return Action::ERR;
      }
    }
  }

  return Action::OK;
}

// Action_MinMaxDist::DoAction()
Action::RetType Action_MinMaxDist::DoAction(int frameNum, ActionFrame& frm)
{
  return Action::OK;
}
