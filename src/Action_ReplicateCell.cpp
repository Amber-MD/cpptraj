#include "Action_ReplicateCell.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"

// CONSTRUCTOR
Action_ReplicateCell::Action_ReplicateCell() : coords_(0), ncopies_(0), writeTraj_(false) {} 

void Action_ReplicateCell::Help() const {
  mprintf("\t[out <traj filename>] [name <dsname>]\n"
          "\t{ all | dir <XYZ> [dir <XYZ> ...] } [<mask>]\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Replicate unit cell in specified (or all) directions for atoms in <mask>.\n"
          "    <XYZ>: X= 1, 0, -1, replicate in specified direction (e.g. 100 is +X only)\n");
  mprintf("%s", ActionTopWriter::Options());
}

static inline int toDigit(char c) {
  switch (c) {
    case '0' : return 0;
    case '1' : return 1;
    case '2' : return 2;
    case '3' : return 3;
    case '4' : return 4;
    case '5' : return 5;
    case '6' : return 6;
    case '7' : return 7;
    case '8' : return 8;
    case '9' : return 9;
  }
  return 0; // SANITY CHECK
}

// Action_ReplicateCell::Init()
Action::RetType Action_ReplicateCell::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Set up output traj
  std::string trajfilename = actionArgs.GetStringKey("out");
  topWriter_.InitTopWriter( actionArgs, "replicated cell", debugIn );
  bool setAll = actionArgs.hasKey("all");
  std::string dsname = actionArgs.GetStringKey("name");
  if (!dsname.empty()) {
    coords_ = (DataSet_Coords*)init.DSL().AddSet(DataSet::COORDS, dsname, "RCELL");
    if (coords_ == 0) return Action::ERR;
  }
  if (trajfilename.empty() && coords_ == 0) {
    mprinterr("Error: Either 'out <traj filename> or 'name <dsname>' must be specified.\n");
    return Action::ERR;
  }
  // Get Mask
  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  // Determine which directions to set
  if (setAll) {
    for (int ix = -1; ix < 2; ix++)
      for (int iy = -1; iy < 2; iy++)
        for (int iz = -1; iz < 2; iz++) {
          directionArray_.push_back( ix );
          directionArray_.push_back( iy );
          directionArray_.push_back( iz );
        }
  } else {
    std::string dirstring = actionArgs.GetStringKey("dir");
    while (!dirstring.empty()) {
      std::vector<int> ixyz(3, -2);
      std::vector<int>::iterator iptr = ixyz.begin();
      for (std::string::const_iterator c = dirstring.begin();
                                       c != dirstring.end(); ++c)
      {
        if (iptr == ixyz.end()) {
          mprinterr("Error: 'dir' string has too many characters.\n");
          return Action::ERR;
        }
        int sign = 1;
        if      (*c == '+') ++c;
        else if (*c == '-') { sign = -1; ++c; }
        
        if (isdigit( *c ))
          *iptr = toDigit( *c ) * sign;
        else {
          mprinterr("Error: illegal character '%c' in 'dir' string '%s'; only numbers allowed.\n",
                    *c, dirstring.c_str());
          return Action::ERR;
        }
        ++iptr;
      }
      //mprintf("DEBUG: %s = %i %i %i\n", dirstring.c_str(), ixyz[0], ixyz[1], ixyz[2]);
      directionArray_.push_back( ixyz[0] );
      directionArray_.push_back( ixyz[1] );
      directionArray_.push_back( ixyz[2] );
      dirstring = actionArgs.GetStringKey("dir");
    }
  }
  ncopies_ = (int)(directionArray_.size() / 3);
  if (ncopies_ < 1) {
    mprinterr("Error: No directions (or 'all') specified.\n");
    return Action::ERR;
  }
  // Initialize output trajectory with remaining arguments
  if (!trajfilename.empty()) {
    outtraj_.SetDebug( debugIn );
    if ( outtraj_.InitEnsembleTrajWrite(trajfilename, actionArgs.RemainingArgs(), init.DSL(),
                                        TrajectoryFile::UNKNOWN_TRAJ, init.DSL().EnsembleNum()) )
      return Action::ERR;
    writeTraj_ = true;
#   ifdef MPI
    outtraj_.SetTrajComm( init.TrajComm() );
#   endif
  } else
    writeTraj_ = false;

  mprintf("    REPLICATE CELL: Replicating cell in %i directions:\n", ncopies_);
  mprintf("\t\t X  Y  Z\n");
  for (unsigned int i = 0; i != directionArray_.size(); i += 3)
    mprintf("\t\t%2i %2i %2i\n", directionArray_[i], 
            directionArray_[i+1], directionArray_[i+2]);
  mprintf("\tUsing atoms in mask '%s'\n", Mask1_.MaskString());
  if (writeTraj_)
    mprintf("\tWriting to trajectory %s\n", outtraj_.Traj().Filename().full());
  topWriter_.PrintOptions();
  if (coords_ != 0)
    mprintf("\tSaving coords to data set %s\n", coords_->legend());

  return Action::OK;
}

// Action_ReplicateCell::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_ReplicateCell::Setup(ActionSetup& setup) {
  if (setup.Top().SetupIntegerMask( Mask1_ )) return Action::ERR;
  mprintf("\t%s (%i atoms)\n",Mask1_.MaskString(), Mask1_.Nselected());
  if (Mask1_.None()) {
    mprintf("Warning: One or both masks have no atoms.\n");
    return Action::SKIP;
  }
  // Check unit cell info for this parm
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Warning: No box info, cannot replica cell for topology %s\n", setup.Top().c_str());
    return Action::SKIP;
  }
  // Create combined topology.
  if (combinedTop_.Natom() > 0) {
    // Topology already set up. Check that # atoms matches.
    if (Mask1_.Nselected() * ncopies_ != combinedTop_.Natom()) {
      mprintf("Warning: Unit cell can currently only be replicated for"
              " topologies with same # atoms.\n");
      return Action::SKIP;
    }
    // Otherwise assume top does not change.
  } else {
    // Set up topology and frame.
    Topology* stripParm = setup.Top().modifyStateByMask( Mask1_ );
    if (stripParm == 0) return Action::ERR;
    for (int cell = 0; cell != ncopies_; cell++)
      combinedTop_.AppendTop( *stripParm );
    combinedTop_.Brief("Combined parm:");
    delete stripParm;
    topWriter_.WriteTops( combinedTop_ );
    // Only coordinates for now. FIXME
    combinedFrame_.SetupFrameM(combinedTop_.Atoms());
    // Set up COORDS / output traj if necessary.
    if (coords_ != 0)
      coords_->CoordsSetup( combinedTop_, CoordinateInfo() );
    if (writeTraj_) {
      if ( outtraj_.SetupTrajWrite( &combinedTop_, CoordinateInfo(), setup.Nframes() ) )
        return Action::ERR;
    }
  }

  return Action::OK;
}

// Action_ReplicateCell::DoAction()
Action::RetType Action_ReplicateCell::DoAction(int frameNum, ActionFrame& frm) {
  int idx, newFrameIdx;
  unsigned int id;
  Vec3 frac, t2;

  int shift = Mask1_.Nselected() * 3;
# ifdef _OPENMP
# pragma omp parallel private(idx, newFrameIdx, id) firstprivate(frac, t2)
  {
# pragma omp for
# endif
  for (idx = 0; idx < Mask1_.Nselected(); idx++) {
    // Convert to fractional coords
    frac = frm.Frm().BoxCrd().FracCell() * Vec3(frm.Frm().XYZ( Mask1_[idx] ));
    //mprintf("DEBUG: Atom %i frac={ %g %g %g }\n", Mask1_[idx]+1, frac[0], frac[1], frac[2]);
    // replicate in each direction
    newFrameIdx = idx * 3;
    for (id = 0; id != directionArray_.size(); id+=3, newFrameIdx += shift)
    {
       // Convert back to Cartesian coords.
       t2 = frm.Frm().BoxCrd().UnitCell().TransposeMult(frac + Vec3(directionArray_[id  ],
                                             directionArray_[id+1],
                                             directionArray_[id+2]));
       combinedFrame_[newFrameIdx  ] = t2[0];
       combinedFrame_[newFrameIdx+1] = t2[1];
       combinedFrame_[newFrameIdx+2] = t2[2];
    }
  }
# ifdef _OPENMP
  }
# endif
  if (writeTraj_) {
    if (outtraj_.WriteSingle(frm.TrajoutNum(), combinedFrame_) !=0 )
      return Action::ERR;
  }
  if (coords_ != 0)
    coords_->AddFrame( combinedFrame_ );
  
  return Action::OK;
}
