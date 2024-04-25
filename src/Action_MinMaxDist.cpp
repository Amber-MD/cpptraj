#include "Action_MinMaxDist.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "StringRoutines.h" // integerToString
#include <limits> // Max double
#include <cmath> // sqrt

/** CONSTRUCTOR */
Action_MinMaxDist::Action_MinMaxDist() :
  mode_(NO_MODE),
  distType_(NO_DIST),
  outfile_(0),
  byAtomSet_(0),
  masterDSL_(0),
  debug_(0),
  resOffset_(0)
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
          "\t[noimage] [name <setname>] [out <file>] [resoffset <#>]\n"
          "  Record the min/max distance between atoms/residues/molecules.\n"
         );
}

// Action_MinMaxDist::Init()
Action::RetType Action_MinMaxDist::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  imageOpt_.InitImaging( !(actionArgs.hasKey("noimage")) );
  resOffset_ = actionArgs.getKeyInt("resoffset", 1);
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
  // File args
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Mode args
  if (actionArgs.hasKey("byatom"))
    mode_ = BY_ATOM;
  else if (actionArgs.hasKey("byres"))
    mode_ = BY_RES;
  else if (actionArgs.hasKey("bymol"))
    mode_ = BY_MOL;
  // Distance calc type args
  // TODO This action is a little different in that it looks at the 
  //      command name to figure out if user wants minimum or maximum
  //      distance. May break if command name is updated and this isnt.
  //bool calc_mindist = (actionArgs.hasKey("mindist") || (actionArgs[0] == "mindist"));
  //bool calc_maxdist = (actionArgs.hasKey("maxdist") || (actionArgs[0] == "maxdist"));
  bool calc_mindist = (actionArgs[0][0] == 'm' && actionArgs[0][1] == 'i' && actionArgs[0][2] == 'n');
  bool calc_maxdist = (actionArgs[0][0] == 'm' && actionArgs[0][1] == 'a' && actionArgs[0][2] == 'x');
  //if (actionArgs.hasKey("bothdist"))
  //  distType_ = BOTH_DIST;
  // DataSet Name
  dsname_ = actionArgs.GetStringKey("name");
  // Default name
  if (dsname_.empty())
    dsname_ = init.DSL().GenerateDefaultName("MINMAXDIST");
  // Default mode
  if (mode_ == NO_MODE)
    mode_ = BY_ATOM;
  // Default distance calc type
  if (distType_ == NO_DIST) {
    if (calc_mindist && calc_maxdist) {
      //distType_ = BOTH_DIST;
      mprinterr("Error: Can only have 'mindist' or 'maxdist', not both.\n");
      return Action::ERR;
    } else if (calc_mindist)
      distType_ = MIN_DIST;
    else if (calc_maxdist)
      distType_ = MAX_DIST;
    else {
      mprinterr("Internal Error: No distance type specified and command name '%s' unrecognized.\n",
                actionArgs[0].c_str());
      return Action::ERR;
    }
  }
  // Checks
  if (mode_ == BY_ATOM && distType_ == MIN_DIST && !mask2_.MaskStringSet())
    mprintf("Warning: Minimum by atom distance for a single mask will always be the shortest bond distance.\n");
  // Save master data set list
  masterDSL_ = init.DslPtr();
  // Allocate by atom set
  if (mode_ == BY_ATOM) {
    byAtomSet_ = init.DSL().AddSet(DataSet::FLOAT, MetaData(dsname_));
    if (byAtomSet_ == 0) {
      mprinterr("Error: Could not allocate set '%s'\n", dsname_.c_str());
      return Action::ERR;
    }
    if (outfile_ != 0)
      outfile_->AddDataSet( byAtomSet_ );
  }

  mprintf("    MINMAXDIST: Calculating %s distance for selected %s.\n",
          distTypeStr_[distType_], modeStr_[mode_]);
  if (mode_ == BY_RES)
    mprintf("\tResidue offset: %i\n", resOffset_);
  mprintf("\tMask1: %s\n", mask1_.MaskString());
  if (mask2_.MaskStringSet()) {
    mprintf("\tMask2: %s\n", mask2_.MaskString());
  }
  if (imageOpt_.UseImage())
    mprintf("\tDistances will use minimum image convention if box info present.\n");
  else
    mprintf("\tDistances will not be imaged.\n");
  mprintf("\tData set name: %s\n", dsname_.c_str());
  if (outfile_ != 0)
    mprintf("\tData file name: %s\n", outfile_->DataFilename().full());

  return Action::OK;
}

/** For DEBUG, print selected atoms in entity array., */
void Action_MinMaxDist::printEntities(Earray const& entities, AtomMask const& maskIn) const {
  mprintf("DEBUG: Selected %s in mask %s\n",
          modeStr_[mode_], maskIn.MaskString());
  for (Earray::const_iterator it = entities.begin(); it != entities.end(); ++it) {
    mprintf("\t%s %i :", it->name_.c_str(), it->num_ + 1);
    for (AtomMask::const_iterator at = it->emask_.begin(); at != it->emask_.end(); ++at)
      mprintf(" %i", *at + 1);
    mprintf("\n");
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
    //entities.reserve(topIn.Nres());
    for (int ires = 0; ires < topIn.Nres(); ires++) {
      Residue const& Res = topIn.Res(ires);
      bool needs_alloc = true;
      for (int at = Res.FirstAtom(); at != Res.LastAtom(); at++) {
        if (cmask.AtomInCharMask( at )) {
          if (needs_alloc) {
            entities.push_back( Entity(Res.Name().Truncated(), ires) );
            needs_alloc = false;
          }
          Entity& currentEntity = entities.back();
          currentEntity.emask_.AddSelectedAtom( at );
        }
      }
    } // END loop over residues
  } else if (mode_ == BY_MOL) {
    for (int imol = 0; imol < topIn.Nmol(); imol++) {
      Molecule const& Mol = topIn.Mol(imol);
      bool needs_alloc = true;
      for (Unit::const_iterator seg = Mol.MolUnit().segBegin();
                                seg != Mol.MolUnit().segEnd(); ++seg)
      {
        for (int at = seg->Begin(); at < seg->End(); ++at) {
          if (cmask.AtomInCharMask( at )) {
            if (needs_alloc) {
              int rnum = topIn[at].ResNum();
              entities.push_back( Entity(topIn.Res(rnum).Name().Truncated(), imol) );
              needs_alloc = false;
            }
            Entity& currentEntity = entities.back();
            currentEntity.emask_.AddSelectedAtom( at );
          }
        }
      } // END loop over molecule segments
    } // END loop over molecules
  } else {
    mprinterr("Internal Error: Action_MinMaxDist::setupEntityArray() Unhandled mode\n");
    return 1;
  }
  // DEBUG - print entities
  if (debug_ > 0)
    printEntities(entities, maskIn);
  mprintf("\t%zu %s selected by mask '%s'\n", entities.size(), modeStr_[mode_], maskIn.MaskString());
  return 0;
}

/// Aspect for BY_RES/BY_MOL
inline static std::string entity_aspect(int r0, int r1)
{
  return std::string( integerToString(r0+1) + "_" + integerToString(r1+1) );
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
  // By residue/molecule interaction pair data set setup
  activeSets_.clear();
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
      // Check if entity numbers overlap
      std::vector<bool> entity_selected;
      if (mode_ == BY_RES)
        entity_selected.assign( setup.Top().Nres(), false );
      else if (mode_ == BY_MOL)
        entity_selected.assign( setup.Top().Nmol(), false );
      for (Earray::const_iterator it1 = entities1_.begin(); it1 != entities1_.end(); ++it1)
        entity_selected[it1->num_] = true;
      for (Earray::const_iterator it2 = entities2_.begin(); it2 != entities2_.end(); ++it2) {
        if (entity_selected[it2->num_]) {
          mprinterr("Error: %s selected by mask2 '%s' must not overlap those selected by mask1 '%s'\n",
                    modeStr_[mode_], mask2_.MaskString(), mask1_.MaskString());
          return Action::ERR;
        }
      }
      // Set up DataSets for each entity pair
      int tgtOffset = 0;
      if (mode_ == BY_RES)
        tgtOffset = resOffset_;
      for (Earray::const_iterator it1 = entities1_.begin(); it1 != entities1_.end(); ++it1)
      {
        for (Earray::const_iterator it2 = entities2_.begin(); it2 != entities2_.end(); ++it2)
        {
          int offset = it1->num_ - it2->num_;
          if (offset < 0) offset = -offset;
          if (offset > tgtOffset) {
            if (debug_ > 1)
              mprintf("DEBUG: Pair %i - %i\n", it1->num_ + 1, it2->num_ + 1);
            MetaData meta(dsname_, entity_aspect(it1->num_, it2->num_));
            DataSet* ds = interactionSets_.AddInteractionSet(*masterDSL_, DataSet::FLOAT, meta, it2->num_, it1->num_, outfile_);
            if (ds == 0) {
              mprinterr("Error: Could not allocate data set %s[%s]\n",
                        meta.Name().c_str(), meta.Aspect().c_str());
              return Action::ERR;
            }
            activeSets_.push_back( ActiveSet(ds, it1, it2) );
          }
        }
      }
    } else {
      // Set up DataSets for each entity pair
      int tgtOffset = 0;
      if (mode_ == BY_RES)
        tgtOffset = resOffset_;
      for (Earray::const_iterator it1 = entities1_.begin(); it1 != entities1_.end(); ++it1)
      {
        for (Earray::const_iterator it2 = it1 + 1; it2 != entities1_.end(); ++it2)
        {
          int offset = it1->num_ - it2->num_;
          if (offset < 0) offset = -offset;
          if (offset > tgtOffset) {
            if (debug_ > 1)
              mprintf("DEBUG: Pair %i - %i\n", it1->num_ + 1, it2->num_ + 1);
            MetaData meta(dsname_, entity_aspect(it1->num_, it2->num_));
            DataSet* ds = interactionSets_.AddInteractionSet(*masterDSL_, DataSet::FLOAT, meta, it2->num_, it1->num_, outfile_);
            if (ds == 0) {
              mprinterr("Error: Could not allocate data set %s[%s]\n",
                        meta.Name().c_str(), meta.Aspect().c_str());
              return Action::ERR;
            }
            activeSets_.push_back( ActiveSet(ds, it1, it2) );
          }
        }
      }
    }
    mprintf("\t%zu active pairs of %s\n", activeSets_.size(), modeStr_[mode_]);
    if (debug_ > 0) {
      mprintf("DEBUG: Active sets:\n");
      for (DSarray::const_iterator it = activeSets_.begin(); it != activeSets_.end(); ++it)
        mprintf("\t%s\n", it->ds_->legend());
    }
    if (activeSets_.empty()) {
      mprintf("Warning: No active interaction pairs. Skipping.\n");
      return Action::SKIP;
    }
  } // END BY_RES, BY_MOL

  return Action::OK;
}

/** Get min distance between all atom pairs */
double Action_MinMaxDist::get_min_dist(AtomMask const& mask1, AtomMask const& mask2, Frame const& frameIn)
const
{
  double min_dist2 = std::numeric_limits<double>::max(); 
  for (AtomMask::const_iterator at1 = mask1.begin(); at1 != mask1.end(); ++at1) {
    const double* XYZ1 = frameIn.XYZ(*at1);
    for (AtomMask::const_iterator at2 = mask2.begin(); at2 != mask2.end(); ++at2) {
      if (*at1 != *at2) {
        const double* XYZ2 = frameIn.XYZ(*at2);
        double dist2 = DIST2(imageOpt_.ImagingType(), XYZ1, XYZ2, frameIn.BoxCrd());
        //mprintf("DEBUG: get_min_dist(m1,m2) %i %i = %f\n", *at1 + 1, *at2 + 1, sqrt(dist2));
        if (dist2 < min_dist2)
          min_dist2 = dist2;
      }
    }
  }
  return sqrt(min_dist2);
}

/** Get min distance between all atoms pairs */
double Action_MinMaxDist::get_min_dist(AtomMask const& mask1, Frame const& frameIn)
const
{
  double min_dist2 = std::numeric_limits<double>::max(); 
  for (AtomMask::const_iterator at1 = mask1.begin(); at1 != mask1.end(); ++at1) {
    const double* XYZ1 = frameIn.XYZ(*at1);
    for (AtomMask::const_iterator at2 = at1 + 1; at2 != mask1.end(); ++at2) {
      const double* XYZ2 = frameIn.XYZ(*at2);
      double dist2 = DIST2(imageOpt_.ImagingType(), XYZ1, XYZ2, frameIn.BoxCrd());
      if (dist2 < min_dist2)
        min_dist2 = dist2;
    }
  }
  return sqrt(min_dist2);
}

/** Get max distance between all atom pairs */
double Action_MinMaxDist::get_max_dist(AtomMask const& mask1, AtomMask const& mask2, Frame const& frameIn)
const
{
  double max_dist2 = 0;
  for (AtomMask::const_iterator at1 = mask1.begin(); at1 != mask1.end(); ++at1) {
    const double* XYZ1 = frameIn.XYZ(*at1);
    for (AtomMask::const_iterator at2 = mask2.begin(); at2 != mask2.end(); ++at2) {
      if (*at1 != *at2) {
        const double* XYZ2 = frameIn.XYZ(*at2);
        double dist2 = DIST2(imageOpt_.ImagingType(), XYZ1, XYZ2, frameIn.BoxCrd());
        if (dist2 > max_dist2)
          max_dist2 = dist2;
      }
    }
  }
  return sqrt(max_dist2);
}

/** Get max distance between all atoms pairs */
double Action_MinMaxDist::get_max_dist(AtomMask const& mask1, Frame const& frameIn)
const
{
  double max_dist2 = 0; 
  for (AtomMask::const_iterator at1 = mask1.begin(); at1 != mask1.end(); ++at1) {
    const double* XYZ1 = frameIn.XYZ(*at1);
    for (AtomMask::const_iterator at2 = at1 + 1; at2 != mask1.end(); ++at2) {
      const double* XYZ2 = frameIn.XYZ(*at2);
      double dist2 = DIST2(imageOpt_.ImagingType(), XYZ1, XYZ2, frameIn.BoxCrd());
      if (dist2 > max_dist2)
        max_dist2 = dist2;
    }
  }
  return sqrt(max_dist2);
}

// Action_MinMaxDist::DoAction()
Action::RetType Action_MinMaxDist::DoAction(int frameNum, ActionFrame& frm)
{
  if (imageOpt_.ImagingEnabled())
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );

  if (mode_ == BY_ATOM) {
    float dist;
    if (mask2_.MaskStringSet()) {
      if (distType_ == MIN_DIST)
        dist = get_min_dist(mask1_, mask2_, frm.Frm());
      else if (distType_ == MAX_DIST)
        dist = get_max_dist(mask1_, mask2_, frm.Frm());
    } else {
       if (distType_ == MIN_DIST)
        dist = get_min_dist(mask1_, frm.Frm());
      else if (distType_ == MAX_DIST)
        dist = get_max_dist(mask1_, frm.Frm());
    }
    byAtomSet_->Add( frameNum, &dist );
  } else {
    // BY_RES / BY_MOL
    for (DSarray::const_iterator set = activeSets_.begin(); set != activeSets_.end(); ++set) {
      float dist;
      if (distType_ == MIN_DIST)
        dist = get_min_dist(set->it1_->emask_, set->it2_->emask_, frm.Frm());
      else if (distType_ == MAX_DIST)
        dist = get_max_dist(set->it1_->emask_, set->it2_->emask_, frm.Frm());
    
      set->ds_->Add( frameNum, &dist );
    }
  }
   
  return Action::OK;
}
