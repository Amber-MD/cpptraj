#include <cmath>
#include <algorithm> // std::min, std::max
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DistRoutines.h"
#include "Constants.h" // RADDEG
#include "Trajout_Single.h"
#include "ParmFile.h"
#ifdef NASTRUCTDEBUG
#include "PDBfile.h" //FIXME remove
#endif
#ifdef MPI
# include "DataSet_float.h" // internal pointer needed for sync
#endif

/** CONSTRUCTOR */
Action_NAstruct::Action_NAstruct() :
  puckerMethod_(NA_Base::ALTONA),
  HBdistCut2_(12.25),     // Hydrogen Bond distance cutoff^2: 3.5^2
  // NOTE: should this be bigger? 3dna cutoff is 15 ang
  originCut2_(25.0),      // Origin cutoff^2 for base-pairing: 5.0^2
  // NOTE: should this be smaller? 3dna cutoff is 1.5 Ang
  staggerCut_(2.0),       // Vertical separation cutoff
  // NOTE: should this be smaller? 3dna cutoff is 30 deg.
  z_angle_cut_(1.134464), // Z angle cutoff in radians (65 deg)
  maxResSize_(0),
  debug_(0),
  nframes_(0),
  findBPmode_(FIRST),
  grooveCalcType_(PP_OO),
  bpConvention_(BP_3DNA),
  printheader_(true),
  seriesUpdated_(false),
  skipIfNoHB_(true),
  spaceBetweenFrames_(true),
  sscalc_(false),
  wc_hb_only_(false),
  bpout_(0),
  ssout_(0),
  stepout_(0),
  helixout_(0),
  masterDSL_(0),
# ifdef NASTRUCTDEBUG
  calcparam_(true),
# endif
  axesOut_(0),
  axesParm_(0),
  bpAxesOut_(0),
  bpAxesParm_(0),
  stepAxesOut_(0),
  stepAxesParm_(0),
  setupNframes_(0),
  setupTop_(0)
{}

/** DESTRUCTOR */
Action_NAstruct::~Action_NAstruct() {
  if (axesOut_ != 0) {
    axesOut_->EndTraj();
    delete axesOut_;
  }
  if (axesParm_ != 0) {
    delete axesParm_;
  }
  if (bpAxesOut_ != 0) {
    bpAxesOut_->EndTraj();
    delete bpAxesOut_;
  }
  if (bpAxesParm_ != 0) {
    delete bpAxesParm_;
  }
  if (stepAxesOut_ != 0) {
    stepAxesOut_->EndTraj();
    delete stepAxesOut_;
  }
  if (stepAxesParm_ != 0) {
    delete stepAxesParm_;
  }
}

/** Set up axes pseudo trajectory. */
int Action_NAstruct::init_axes_pseudoTraj(const char* description,
                                          const char* trajKeyword,
                                          const char* argKeyword,
                                          const char* parmKeyword,
                                          const char* topName,
                                          DataSetList const& DSL,
                                          ArgList& actionArgs,
                                          Trajout_Single** outputTraj,
                                          Topology** outputParm)
const
{
  std::string axesout = actionArgs.GetStringKey(trajKeyword);
  if (!axesout.empty()) {
    *outputTraj = new Trajout_Single();
    (*outputTraj)->SetDebug( debug_ );
#   ifdef MPI
    (*outputTraj)->SetTrajComm( trajComm_ );
#   endif
    std::string axesoutargStr;
    std::string axesoutarg = actionArgs.GetStringKey(argKeyword);
    while (!axesoutarg.empty()) {
      axesoutargStr.append(" " + axesoutarg);
      axesoutarg = actionArgs.GetStringKey(argKeyword);
    }
    ArgList axesOutArglist( axesoutargStr );
    if ((*outputTraj)->InitEnsembleTrajWrite( axesout, axesOutArglist, DSL,
                                         TrajectoryFile::UNKNOWN_TRAJ, DSL.EnsembleNum()) )
    {
      mprinterr("Error: Could not init %s trajectory '%s'\n", description, axesout.c_str());
      return 1;
    }
    std::string axesparmout = actionArgs.GetStringKey(parmKeyword);
    *outputParm = new Topology();
    (*outputParm)->SetDebug( debug_ );
    (*outputParm)->SetParmName( topName, axesparmout );
  }
  return 0;
}

/** Print help text. */
void Action_NAstruct::Help() const {
  mprintf("\t[<dataset name>] [resrange <range>] [sscalc] [naout <suffix>]\n"
          "\t[noheader] [resmap <ResName>:{A,C,G,T,U} ...] [calcnohb]\n"
          "\t[noframespaces] [baseref <file>] ...\n"
          "\t[bpmode {3dna|babcock}] [allhb]\n"
          "\t[hbcut <hbcut>] [origincut <origincut>] [altona | cremer]\n"
          "\t[zcut <zcut>] [zanglecut <zanglecut>] [groovecalc {simple | 3dna}]\n"
          "\t[axesout <file> [axesoutarg <arg> ...] [axesparmout <file>]]\n"
          "\t[bpaxesout <file> [bpaxesoutarg <arg> ...] [bpaxesparmout <file>]]\n"
          "\t[stepaxesout <file> [stepaxesoutarg <arg> ...] [stepaxesparmout <file>]]\n"
          "\t[axisnameo <name>] [axisnamex <name>] [axisnamey <name>] [axisnamez <name>]\n"
          "\t[{ %s |\n"
          "\t   allframes |\n"
          "\t   specifiedbp pairs <b1>-<b2>,... }]\n", DataSetList::RefArgs);
  mprintf("  Perform nucleic acid structure analysis. Base pairing can be determined\n"
          "  in multiple ways:\n"
          "    - If 'first' (default) or a reference is specified, determine base\n"
          "      pairing using geometric criteria in a manner similar to 3DNA.\n"
          "    - If 'allframes' is specified, base pairing will be determined\n"
          "      using geometric criteria for every single frame.\n"
          "    - If 'specifiedbp' is specified, base pairing is given by subsequent\n"
          "      'pairs <b1>-<b2>,...' arguments, where <b1> and <b2> are the residue\n"
          "      numbers of bases in the base pair, e.g. 'pairs 1-16,2-15,3-14,4-13'.\n"
          "    - If 'reference', 'ref', or 'refindex' is specified, use a reference\n"
          "      structure to determine base pairing.\n"
          "  If 'calcnohb' is specified NA parameters will be calculated even if no\n"
          "  hydrogen bonds present between base pairs.\n"
          "  If 'allhb' is specified report the total number of hydrogen bonds detected\n"
          "  instead of just the number of Watson-Crick-Franklin hydrogen bonds.\n"
          "  Base pair parameters are written to 'BP.<suffix>', base pair step parameters\n"
          "  are written to 'BPstep.<suffix>', and helix parameters are written to\n"
          "  Helix.<suffix>'.\n"
          "  If 'noframespaces' is specified there will be no spaces between frames\n"
          "  in the 'naout' files.\n");
}

// Action_NAstruct::Init()
Action::RetType Action_NAstruct::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  masterDSL_ = init.DslPtr();
  // Get keywords
  sscalc_ = actionArgs.hasKey("sscalc");
  std::string outputsuffix = actionArgs.GetStringKey("naout");
  if (!outputsuffix.empty()) {
    // Set up output files.
    FileName FName( outputsuffix );
    bpout_ = init.DFL().AddCpptrajFile(FName.PrependFileName("BP."), "Base Pair");
    stepout_ = init.DFL().AddCpptrajFile(FName.PrependFileName("BPstep."), "Base Pair Step");
    helixout_ = init.DFL().AddCpptrajFile(FName.PrependFileName("Helix."), "Helix");
    if (bpout_ == 0 || stepout_ == 0 || helixout_ == 0) return Action::ERR;
    if (sscalc_) {
      ssout_ = init.DFL().AddCpptrajFile(FName.PrependFileName("SS."), "Single Strand");
      if (ssout_ == 0) return Action::ERR;
    }
  }
  std::string bpmode = actionArgs.GetStringKey("bpmode");
  if (bpmode.empty())
    bpConvention_ = BP_3DNA;
  else if (bpmode == "3dna")
    bpConvention_ = BP_3DNA;
  else if (bpmode == "babcock")
    bpConvention_ = BP_BABCOCK;
  else {
    mprinterr("Error: Unrecognized keyword '%s' for 'bpmode'.\n", bpmode.c_str());
    return Action::ERR;
  }
  wc_hb_only_ = !actionArgs.hasKey("allhb");
  double hbcut = actionArgs.getKeyDouble("hbcut", -1);
  if (hbcut > 0) 
    HBdistCut2_ = hbcut * hbcut;
  double origincut = actionArgs.getKeyDouble("origincut", -1);
  if (origincut > 0)
    originCut2_ = origincut * origincut;
  double zcut = actionArgs.getKeyDouble("zcut", -1);
  if (zcut > 0)
    staggerCut_ = zcut;
  double zanglecut_deg = actionArgs.getKeyDouble("zanglecut", -1);
  if (zanglecut_deg > 0)
    z_angle_cut_ = zanglecut_deg * Constants::DEGRAD;
  std::string groovecalc = actionArgs.GetStringKey("groovecalc");
  if (!groovecalc.empty()) {
    if (groovecalc == "simple") grooveCalcType_ = PP_OO;
    else if (groovecalc == "3dna") grooveCalcType_ = HASSAN_CALLADINE;
    else {
      mprinterr("Error: Invalid value for 'groovecalc' %s; expected simple or 3dna.\n",
                groovecalc.c_str());
      return Action::ERR;
    }
  } else
    grooveCalcType_ = PP_OO;
  if      (actionArgs.hasKey("altona")) puckerMethod_=NA_Base::ALTONA;
  else if (actionArgs.hasKey("cremer")) puckerMethod_=NA_Base::CREMER;
  // See if we want axes pseudo-trajectories
  if (init_axes_pseudoTraj("base axes", "axesout", "axesoutarg", "axesparmout",
                           "BaseAxes", init.DSL(), actionArgs,
                           &axesOut_, &axesParm_))
    return Action::ERR;
  if (init_axes_pseudoTraj("basepair axes", "bpaxesout", "bpaxesoutarg", "bpaxesparmout",
                           "BasePairAxes", init.DSL(), actionArgs,
                           &bpAxesOut_, &bpAxesParm_))
    return Action::ERR;
  if (init_axes_pseudoTraj("step axes", "stepaxesout", "stepaxesoutarg", "stepaxesparmout",
                           "StepAxes", init.DSL(), actionArgs,
                           &stepAxesOut_, &stepAxesParm_))
    return Action::ERR;
  axisNameO_ = actionArgs.GetStringKey("axisnameo", "Orig");
  axisNameX_ = actionArgs.GetStringKey("axisnamex", "X");
  axisNameY_ = actionArgs.GetStringKey("axisnamey", "Y");
  axisNameZ_ = actionArgs.GetStringKey("axisnamez", "Z");
  // Get residue range
  resRange_.SetRange(actionArgs.GetStringKey("resrange"));
  if (!resRange_.Empty())
    resRange_.ShiftBy(-1); // User res args start from 1
  printheader_ = !actionArgs.hasKey("noheader");
  skipIfNoHB_ = !actionArgs.hasKey("calcnohb");
  spaceBetweenFrames_ = !actionArgs.hasKey("noframespaces");
  // Determine how base pairs will be found.
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty())
    findBPmode_ = REFERENCE;
  else if (actionArgs.hasKey("allframes"))
    findBPmode_ = ALL;
  else if (actionArgs.hasKey("specifiedbp"))
    findBPmode_ = SPECIFIED;
  else if (actionArgs.hasKey("guessbp")) {
    mprinterr("Error: 'guessbp' is deprecated. Consider using 'specifiedbp' instead.\n");
    return Action::ERR;
  } else if (actionArgs.hasKey("first"))
    findBPmode_ = FIRST;
  else 
    findBPmode_ = FIRST;
# ifdef MPI
  if (findBPmode_ == ALL && trajComm_.Size() > 1) {
    mprinterr("Error: Currently 'allframes' does not work with > 1 process per trajectory"
              " (currently %i)\n", trajComm_.Size());
    return Action::ERR;
  }
# endif
  // Check for user-specified base pairs
  if (findBPmode_ == SPECIFIED) {
    std::string pairsarg = actionArgs.GetStringKey("pairs");
    // Format is b1-b2,...
    while (!pairsarg.empty()) {
      ArgList pairslist( pairsarg, "," );
      for (int iarg = 0; iarg < pairslist.Nargs(); iarg++) {
        ArgList bpair( pairslist[iarg], "-" );
        if (bpair.Nargs() != 2) {
          mprinterr("Error: Malformed base pair argument: %s\n", pairslist[iarg].c_str());
          return Action::ERR;
        }
        if (!validInteger(bpair[0])) {
          mprinterr("Error: Expected an integer, got '%s'\n", bpair[0].c_str());
          return Action::ERR;
        }
        if (!validInteger(bpair[1])) {
          mprinterr("Error: Expected an integer, got '%s'\n", bpair[1].c_str());
          return Action::ERR;
        }
        int b1idx = convertToInteger(bpair[0]);
        int b2idx = convertToInteger(bpair[1]);
        if (b1idx < 1 || b2idx < 1) {
          mprinterr("Error: Base pair arg '%s', base #s must be > 0.\n", pairslist[iarg].c_str());
          return Action::ERR;
        }
        specifiedPairs_.push_back( std::pair<unsigned int,unsigned int>( b1idx, b2idx ) );
      } // END loop over specified pairs
      pairsarg = actionArgs.GetStringKey("pairs");
    } // END checking for 'pairs' keywords
    if (specifiedPairs_.empty()) {
      mprinterr("Error: No 'pairs' arguments for 'specifiedbp'\n");
      return Action::ERR;
    }
  }
  // Check for base pair mode incompatibilities
  if (findBPmode_ == ALL) {
    if (bpAxesOut_ != 0) {
      mprinterr("Error: Cannot use 'allframes' mode with 'bpaxesout' since the # of\n"
                "Error:  base pairs can change each frame.\n");
      return Action::ERR;
    }
    if (stepAxesOut_ != 0) {
      mprinterr("Error: Cannot use 'allframes' mode with 'stepaxesout' since the # of\n"
                "Error:  base pair steps can change each frame.\n");
      return Action::ERR;
    }
  }
  // For guess/specify modes, get base pairing type
  std::string bptype = actionArgs.GetStringKey("bptype");
  while (!bptype.empty()) {
    mprintf("Warning: 'bptype' is deprecated. Ignoring.\n");
    bptype = actionArgs.GetStringKey("bptype");
  }
  // Get custom residue maps
  ArgList maplist;
  NA_Base::NAType mapbase;
  while ( actionArgs.Contains("resmap") ) {
    // Split maparg at ':'
    maplist.SetList( actionArgs.GetStringKey("resmap"), ":" );
    // Expect only 2 args
    if (maplist.Nargs()!=2) {
      mprinterr("Error: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",
                maplist.ArgLine());
      return Action::ERR;
    }
    // Check that second arg is A,C,G,T,or U
    if      (maplist[1] == "A") mapbase = NA_Base::ADE;
    else if (maplist[1] == "C") mapbase = NA_Base::CYT;
    else if (maplist[1] == "G") mapbase = NA_Base::GUA;
    else if (maplist[1] == "T") mapbase = NA_Base::THY;
    else if (maplist[1] == "U") mapbase = NA_Base::URA;
    else {
      mprinterr("Error: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",
                maplist.ArgLine());
      return Action::ERR;
    }
    // Check that residue name is <= 4 chars
    if (maplist[0].size() > 4) {
      mprinterr("Error: resmap resname > 4 chars (%s)\n",maplist.ArgLine());
      return Action::ERR;
    }
    std::pair<RefMapType::iterator, bool> ret = 
      nameToRef_.insert( std::pair<std::string,NA_Base::NAType>(maplist[0], mapbase ) );
    if (!ret.second) {
      mprinterr("Error: Already tried to map residue '%s'\n", maplist[0].c_str());
      return Action::ERR;
    }
  }
  // Get custom base references
  while ( actionArgs.Contains("baseref") ) {
    std::string brefname = actionArgs.GetStringKey("baseref");
    if ( refBases_.LoadFromFile( brefname ) ) return Action::ERR;
  }
  // Get Masks
  // DataSet name
  dataname_ = actionArgs.GetStringNext();

  mprintf("    NAstruct: ");
  if (resRange_.Empty())
    mprintf("Scanning all NA residues\n");
  else
    mprintf("Scanning residues %s\n",resRange_.RangeArg());
  if (bpout_ != 0) {
    mprintf("\tBase pair parameters written to %s\n", bpout_->Filename().full());
    mprintf("\tBase pair step parameters written to %s\n", stepout_->Filename().full());
    mprintf("\tHelical parameters written to %s\n", helixout_->Filename().full());
    if (!printheader_) mprintf("\tHeader line will not be written.\n");
    if (!spaceBetweenFrames_) mprintf("\tNo spaces will be written between frames.\n");
  }
  if (sscalc_) {
    mprintf("\tWill determine parameters for consecutive bases in strands.\n");
    if (ssout_ != 0)
      mprintf("\tSingle strand parameters written to %s\n", ssout_->Filename().full());
  }
  switch (bpConvention_) {
    case BP_3DNA : mprintf("\tUsing 3DNA conventions for base pairing (no XY flip for parallel strands).\n"); break;
    case BP_BABCOCK : mprintf("\tUsing Babcock et al. conventions for base pairing (XY flip for parallel strands).\n"); break;
  }
  mprintf("\tHydrogen bond cutoff for determining base pairs is %.2f Angstroms.\n",
          sqrt( HBdistCut2_ ) );
  if (wc_hb_only_)
    mprintf("\tOnly reporting total # of Watson-Crick-Franklin hydrogen bonds.\n");
  else
    mprintf("\tReporting total # of all hydrogen bonds.\n");
  mprintf("\tBase reference axes origin cutoff for determining base pairs is %.2f Angstroms.\n",
          sqrt( originCut2_ ) );
  mprintf("\tBase Z height cutoff (stagger) for determining base pairs is %.2f Angstroms.\n",
          staggerCut_);
  mprintf("\tBase Z angle cutoff for determining base pairs is %.2f degrees.\n",
          z_angle_cut_ * Constants::RADDEG);
  if (!nameToRef_.empty()) {
    static const char natypeNames[] = { '?', 'A', 'C', 'G', 'T', 'U' };
    mprintf("\tWill attempt to map the following residues to existing references:\n");
    for (RefMapType::const_iterator it = nameToRef_.begin();
                                    it != nameToRef_.end(); ++it)
      mprintf("\t  %-8s : %c\n", it->first.c_str(), natypeNames[it->second]);
  }
  if (findBPmode_ == REFERENCE) {
    // Use reference to determine base pairing
    mprintf("\tUsing reference %s to determine base-pairing.\n", REF.refName());
    ActionSetup ref_setup(REF.ParmPtr(), REF.CoordsInfo(), 1);
    if (Setup( ref_setup )) return Action::ERR;
    // Set up base axes
    if ( SetupBaseAxes(REF.Coord()) ) return Action::ERR;
    // Determine Base Pairing
    if ( DetermineBasePairing() ) return Action::ERR;
    mprintf("\tSet up %zu base pairs.\n", BasePairs_.size() );
  } else if (findBPmode_ == SPECIFIED) {
    mprintf("\tUser specified base pairs:");
    for (PairArray::const_iterator it = specifiedPairs_.begin(); it != specifiedPairs_.end(); ++it)
      mprintf(" %u-%u", it->first, it->second);
    mprintf("\n");
  } else if (findBPmode_ == ALL)
    mprintf("\tBase pairs will be determined for each frame.\n");
  else if (findBPmode_ == FIRST)
    mprintf("\tBase pairs will be determined from first frame.\n");
  if (skipIfNoHB_)
    mprintf("\tParameters will not be calculated when no hbonds present between base pairs.\n");
  else
    mprintf("\tParameters will be calculated between base pairs even when no hbonds present.\n");
  if (puckerMethod_==NA_Base::ALTONA)
    mprintf("\tCalculating sugar pucker using Altona & Sundaralingam method.\n");
  else if (puckerMethod_==NA_Base::CREMER)
    mprintf("\tCalculating sugar pucker using Cremer & Pople method.\n");
  if (grooveCalcType_ == PP_OO)
    mprintf("\tUsing simple groove width calculation (P-P and O-O base pair distances).\n");
  else if (grooveCalcType_ == HASSAN_CALLADINE)
    mprintf("\tUsing groove width calculation of El Hassan & Calladine.\n");
  if (axesOut_ != 0) {
    mprintf("\tWriting base axes pseudo trajectory to '%s'\n", axesOut_->Traj().Filename().full());
    if (!axesParm_->OriginalFilename().empty())
      mprintf("\tWriting base axes pseudo topology to '%s'\n", axesParm_->OriginalFilename().full());
  }
  if (bpAxesOut_ != 0) {
    mprintf("\tWriting base pair axes pseudo trajectory to '%s'\n", bpAxesOut_->Traj().Filename().full());
    if (!bpAxesParm_->OriginalFilename().empty())
      mprintf("\tWriting base pair axes pseudo topology to '%s'\n", bpAxesParm_->OriginalFilename().full());
  }
  if (stepAxesOut_ != 0) {
    mprintf("\tWriting base pair step axes pseudo trajectory to '%s'\n", stepAxesOut_->Traj().Filename().full());
    if (!stepAxesParm_->OriginalFilename().empty())
      mprintf("\tWriting base pair step axes pseudo topology to '%s'\n", stepAxesParm_->OriginalFilename().full());
  }
  if (axesOut_ != 0 || bpAxesOut_ != 0 || stepAxesOut_ != 0) {
    mprintf("\tAxes pseudo atom names: origin='%s' X='%s' Y='%s' Z='%s'\n",
            *axisNameO_, *axisNameX_, *axisNameY_, *axisNameZ_);
  }
  mprintf("# Citations: Babcock MS; Pednault EPD; Olson WK; \"Nucleic Acid Structure\n"
          "#             Analysis: Mathematics for Local Cartesian and Helical Structure\n"
          "#             Parameters That Are Truly Comparable Between Structures\",\n"
          "#             J. Mol. Biol. (1994) 237, 125-156.\n"
          "#            Olson WK; Bansal M; Burley SK; Dickerson RE; Gerstein M;\n"
          "#             Harvey SC; Heinemann U; Lu XJ; Neidle S; Shekked Z; Sklenar H;\n"
          "#             Suzuki M; Tung CS; Westhof E; Wolberger C; Berman H; \"A Standard\n"
          "#             Reference Frame for the Description of Nucleic Acid Base-pair\n"
          "#             Geometry\", J. Mol. Biol. (2001) 313, 229-237.\n");
  if (grooveCalcType_ == HASSAN_CALLADINE)
    mprintf("#            El Hassan MA; Calladine CR; \"Two Distinct Modes of\n"
            "#             Protein-induced Bending in DNA.\"\n"
            "#             J. Mol. Biol. (1998) 282, 331-343.\n");
  init.DSL().SetDataSetsPending(true);
  return Action::OK;
}

// -----------------------------------------------------------------------------
#ifdef NASTRUCTDEBUG
/// Write given NA_Axis to a PDB file.
static void WriteAxes(PDBfile& outfile, int resnum, const char* resname, NA_Axis const& axis)
{
  // Origin
  Vec3 oxyz = axis.Oxyz();
  outfile.WriteATOM("Orig", resnum, oxyz[0], oxyz[1], oxyz[2], resname, 0.0);
  // X vector
  Vec3 vec = axis.Rx() + oxyz;
  outfile.WriteATOM("X", resnum, vec[0], vec[1], vec[2], resname, 0.0);
  // Y vector
  vec = axis.Ry() + oxyz;
  outfile.WriteATOM("Y", resnum, vec[0], vec[1], vec[2], resname, 0.0);
  // Z vector
  vec = axis.Rz() + oxyz;
  outfile.WriteATOM("Z", resnum, vec[0], vec[1], vec[2], resname, 0.0);
}
// -----------------------------------------------------------------------------
#endif

/// Add NA_Axis to given Frame
static void axesToFrame(Frame& frame, NA_Axis const& axis) {
  // Origin
  Vec3 oxyz = axis.Oxyz();
  frame.AddVec3( oxyz );
  // X vector
  frame.AddVec3( axis.Rx() + oxyz );
  // Y vector
  frame.AddVec3( axis.Ry() + oxyz );
  // Z vector
  frame.AddVec3( axis.Rz() + oxyz );
}

// Action_NAstruct::SetupBaseAxes()
/** For each residue in Bases (set up in Setup()), get the corresponding input
  * coords and fit the reference coords on top of input coords. This sets up 
  * the reference axes for each base.
  */
int Action_NAstruct::SetupBaseAxes(Frame const& InputFrame) {
  Frame refFrame(maxResSize_); // Hold copy of base reference coords for RMS fit
  Frame inpFrame(maxResSize_); // Hold copy of input base coords for RMS fit
# ifdef NASTRUCTDEBUG
  PDBfile basesfile;
  basesfile.OpenWrite("bases.pdb");
  mprintf("\n=================== Setup Base Axes ===================\n");
# endif
  for (std::vector<NA_Base>::iterator base = Bases_.begin(); 
                                      base != Bases_.end(); ++base)
  {
    // Set input coords for entire NA residue. 
    base->SetInputFrame( InputFrame );
    // Set input coords for RMS fit.
    inpFrame.SetCoordinates( base->Input(), base->InputFitMask() );
    // Set ref coords for RMS fit. 
    refFrame.SetCoordinates( base->Ref(), base->RefFitMask() );
#   ifdef NASTRUCTDEBUG
    mprintf("Base %i:%4s\n", base->ResNum()+1, base->ResName()); 
    base->InputFitMask().PrintMaskAtoms("InpMask");
    base->RefFitMask().PrintMaskAtoms("RefMask");
    mprintf("%-2s %4s %8s %8s %8s %2s %8s %8s %8s\n","#","Atom","Ex","Ey","Ez","#","Rx","Ry","Rz");
    AtomMask::const_iterator refatom = base->RefFitMask().begin();
    for (AtomMask::const_iterator inpatom = base->InputFitMask().begin();
                                  inpatom != base->InputFitMask().end(); ++inpatom)
    {
      const double* XYZ = base->Input().XYZ(*inpatom);
      mprintf("%-2i %4s %8.3f %8.3f %8.3f", *inpatom+1, base->atomName(*inpatom),
              XYZ[0], XYZ[1], XYZ[2]);
      XYZ = base->Ref().XYZ(*refatom);
      mprintf(" %2i %8.3f %8.3f %8.3f\n", *refatom+1, XYZ[0], XYZ[1], XYZ[2]);
      ++refatom;
    }
#   endif 
    /* Now that we have a set of reference coords and the corresponding input
     * coords, RMS fit the reference coords to the input coords to obtain the
     * appropriate rotation and translations that will put the reference coords 
     * on top of input (experimental) coords. Per 3DNA procedure, not all 
     * reference atoms are used in the RMS fit; only ring atoms are used. 
     */ 
    Matrix_3x3 RotMatrix;
    Vec3 TransVec, refTrans;
    double rmsd = refFrame.RMSD( inpFrame, RotMatrix, TransVec, refTrans, false);
    /* RotMatrix and TransVec now contain rotation and translation
     * that will orient refcoord to expframe. The first translation is that of
     * the reference frame to the absolute origin, the second translation is
     * that of the reference frame to the exp. coords after rotation.
     * The rotation matrix contains the coordinates of the X, Y, and Z unit 
     * vectors of the base axes.
     */
    // Store the Rotation matrix and the rotated and translated origin.
    base->Axis().StoreRotMatrix( RotMatrix, (RotMatrix*TransVec)+refTrans );
    if (debug_>0) { 
      mprintf("Base %i: RMS of RefCoords from ExpCoords is %f\n",base->ResNum(), rmsd);
      base->Axis().PrintAxisInfo("BaseAxes");
    }

#   ifdef NASTRUCTDEBUG
    // Overlap ref coords onto input coords.
    Frame reftemp = base->Ref(); 
    reftemp.Trans_Rot_Trans(TransVec, RotMatrix, refTrans);
    // DEBUG - Write reference frame to file
    for (int i = 0; i < reftemp.Natom(); i++) {
      const double* XYZ = reftemp.XYZ(i);
      basesfile.WriteATOM( base->RefName(i), base->ResNum()+1, XYZ[0], XYZ[1], XYZ[2], 
                           base->ResName(), 0.0 );
    }
#   endif
  } // END loop over bases

  return 0;
}

// -----------------------------------------------------------------------------
Action_NAstruct::HbondType Action_NAstruct::GCpair(NA_Base const& bG, int ig,
                                                   NA_Base const& bC, int ic)
{
  if (bG.AtomName(ig) == "O6" && bC.AtomName(ic) == "N4") return WC;
  if (bG.AtomName(ig) == "N1" && bC.AtomName(ic) == "N3") return WC;
  if (bG.AtomName(ig) == "N2" && bC.AtomName(ic) == "O2") return WC;
  return OTHER;
}

Action_NAstruct::HbondType Action_NAstruct::ATpair(NA_Base const& bA, int ia,
                                                   NA_Base const& bT, int it)
{
  if (bA.AtomName(ia) == "N6" && bT.AtomName(it) == "O4") return WC;
  if (bA.AtomName(ia) == "N1" && bT.AtomName(it) == "N3") return WC;
  return OTHER;
}

Action_NAstruct::HbondType Action_NAstruct::ID_HBtype(NA_Base const& base1, int b1,
                                                      NA_Base const& base2, int b2)
{
  if      ( base1.Type() == NA_Base::GUA && base2.Type() == NA_Base::CYT )
    return (GCpair(base1, b1, base2, b2));
  else if ( base1.Type() == NA_Base::CYT && base2.Type() == NA_Base::GUA )
    return (GCpair(base2, b2, base1, b1));
  else if ( base1.Type() == NA_Base::ADE && base2.Type() == NA_Base::THY )
    return (ATpair(base1, b1, base2, b2));
  else if ( base1.Type() == NA_Base::THY && base2.Type() == NA_Base::ADE )
    return (ATpair(base2, b2, base1, b1));
  else if ( base1.Type() == NA_Base::ADE && base2.Type() == NA_Base::URA ) // A-U has same WC pattern as A-T
    return (ATpair(base1, b1, base2, b2));
  else if ( base1.Type() == NA_Base::URA && base2.Type() == NA_Base::ADE ) // A-U has same WC pattern as A-T
    return (ATpair(base2, b2, base1, b1));
  return OTHER;
}
    
/** Given two NA_Bases for which IDs have been given and input coords set,
  * calculate the number of hydrogen bonds between them.
  */
// TODO Identify type of base pairing (WCF, Hoog., etc)
int Action_NAstruct::CalcNumHB(NA_Base const& base1, NA_Base const& base2, int& n_WC) {
  int Nhbonds = 0;
  n_WC = 0;

  for (int b1 = 0; b1 != base1.Natom(); b1++) {
    if ( base1.HbondType(b1) != NA_Base::NONE ) {
      const double* xyz1 = base1.HBxyz(b1);
      for (int b2 = 0; b2 != base2.Natom(); b2++) {
        if ( base2.HbondType(b2) != NA_Base::NONE &&
             base2.HbondType(b2) != base1.HbondType(b1) )
        {
          const double* xyz2 = base2.HBxyz(b2);
          double dist2 = DIST2_NoImage(xyz1, xyz2);
#         ifdef NASTRUCTDEBUG
          mprintf("\t\t%s:%s -- %s:%s = %f",
                    base1.ResName(), base1.atomName(b1),
                    base2.ResName(), base2.atomName(b2), sqrt(dist2));
#         endif
          if (dist2 < HBdistCut2_) {
            ++Nhbonds;
            HbondType hbtype = ID_HBtype(base1, b1, base2, b2);
            if (hbtype == WC) n_WC++;
#           ifdef NASTRUCTDEBUG
            mprintf(" (%i)", (int)hbtype);
#           endif
          }
#         ifdef NASTRUCTDEBUG
          mprintf("\n");
#         endif
        }
      }
    }
  }
  return Nhbonds;
}

// -----------------------------------------------------------------------------
Action_NAstruct::BPmap::iterator
  Action_NAstruct::AddBasePair(int base1idx, NA_Base const& base1,
                               int base2idx, NA_Base const& base2)
{
  Rpair respair(base1.ResNum(), base2.ResNum());
  // Bases are paired. Try to find existing base pair.
  BPmap::iterator entry = BasePairs_.lower_bound( respair );
  if (entry == BasePairs_.end() || entry->first != respair) {
    // New base pair
#   ifdef NASTRUCTDEBUG
    mprintf("      New base pair: %i to %i", base1.ResNum()+1, base2.ResNum()+1);
#   endif
    MetaData md(dataname_, BasePairs_.size() + 1); // Name, index
    md.SetLegend( base1.BaseName() + base2.BaseName() );
    BPtype BP;
    md.SetAspect("shear");
    BP.shear_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    md.SetAspect("stretch");
    BP.stretch_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    md.SetAspect("stagger");
    BP.stagger_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    md.SetAspect("buckle");
    BP.buckle_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    md.SetAspect("prop");
    BP.prop_    = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    md.SetAspect("open");
    BP.opening_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    md.SetAspect("hb");
    BP.hbonds_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::INTEGER, md);
    md.SetAspect("bp");
    BP.isBP_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::INTEGER, md);
    if (grooveCalcType_ == PP_OO) {
      md.SetAspect("major");
      BP.major_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
      md.SetAspect("minor");
      BP.minor_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
    } else {
      BP.major_ = 0;
      BP.minor_ = 0;
    }
#   ifdef NASTRUCTDEBUG
    //md.SetAspect("oxyz");
    //BP.axes_oxyz_ = (DataSet*)masterDSL_->AddSet(DataSet::VECTOR, md);
    md.SetAspect("nxyz");
    BP.axes_nxyz_ = (DataSet*)masterDSL_->AddSet(DataSet::VECTOR, md);
#   endif
    BP.bpidx_ = BasePairs_.size();
    BP.base1idx_ = base1idx;
    BP.base2idx_ = base2idx;
    entry = BasePairs_.insert( entry, std::pair<Rpair, BPtype>(respair, BP) );
  }
# ifdef NASTRUCTDEBUG
  else
    mprintf("      Existing base pair: %i to %i", base1.ResNum()+1, base2.ResNum()+1);
# endif
  return entry;
}

/** Search Bases_ for the specified residue number, return index into Bases_. */
int Action_NAstruct::find_index_in_bases(int resnum) const {
  for (Barray::const_iterator it = Bases_.begin(); it != Bases_.end(); ++it) {
    if (it->ResNum() == resnum) {
      return (int)(it - Bases_.begin());
    }
  }
  return -1;
}

/** User-specified base pairing. */
int Action_NAstruct::SpecifiedBasePairing() {
  int n_wc_hb;
# ifdef NASTRUCTDEBUG  
  mprintf("\n=================== Specified Base Pairing ===================\n");
# endif
  for (PairArray::const_iterator it = specifiedPairs_.begin();
                                 it != specifiedPairs_.end();
                               ++it)
  {
    // User-specified base pair #s start from 1
    int b1idx = find_index_in_bases(it->first - 1);
    int b2idx = find_index_in_bases(it->second - 1);
    if (b1idx < 0) {
      mprinterr("Error: Specified base residue # %u not found in set up bases.\n", it->first);
      return 1;
    }
    NA_Base& base1 = Bases_[b1idx];
    if (b2idx < 0) {
      mprinterr("Error: Specified base residue # %u not found in set up bases.\n", it->second);
      return 1;
    }
    NA_Base& base2 = Bases_[b2idx];
/*#   ifdef NASTRUCTDEBUG
    // Glycosidic N-N distance
    if (base1.HasNXatom() && base2.HasNXatom()) {
      double n_n_dist2 = DIST2_NoImage(base1.NXxyz(), base2.NXxyz());
      mprintf("DEBUG: NX-NX distance= %f\n", sqrt(n_n_dist2));
    }
#   endif*/
    NA_Axis b1Axis = base1.Axis();
    NA_Axis b2Axis = base2.Axis();
    // Determine if base Z axis vectors are aligned with strand direction
    bool is_z = false;
    int b1_5to3 = axis_points_5p_to_3p( base1 );
    int b2_5to3 = axis_points_5p_to_3p( base2 );
    // TODO trap errors here
    // If antiparallel and both bases are aligned 3' to 5', may be ZDNA
    if (b1_5to3 == 0 && b2_5to3 == 0) {
#     ifdef NASTRUCTDEBUG
      mprintf("Both bases aligned 3' to 5', ZDNA\n");
#     endif
      b1Axis.FlipXZ();
      b2Axis.FlipXZ();
      is_z = true;
    }
    // Determine if base Z axis vectors point in same (theta <= 90) or
    // opposite (theta > 90) directions.
    bool is_antiparallel;
    double z_theta = b1Axis.Rz().Angle( b2Axis.Rz() );
#   ifdef NASTRUCTDEBUG
    double z_deviation_from_linear;
#   endif
    if (z_theta > Constants::PIOVER2) { // If theta(Z) > 90 deg.
#     ifdef NASTRUCTDEBUG
      mprintf("\t%s is anti-parallel to %s (%g deg)\n", base1.ResName(), base2.ResName(),
              z_theta * Constants::RADDEG);
      z_deviation_from_linear = Constants::PI - z_theta;
#     endif
      is_antiparallel = true;
      // Antiparallel - flip Y and Z axes of complimentary base
      b2Axis.FlipYZ();
    } else {
#     ifdef NASTRUCTDEBUG
      mprintf("\t%s is parallel to %s (%g deg)\n", base1.ResName(), base2.ResName(),
              z_theta * Constants::RADDEG);
      z_deviation_from_linear = z_theta;
#     endif
      is_antiparallel = false;
      // Parallel - no flip needed if 3dna.
      // If using Babcock convention, flip X and Y axes.
      if (bpConvention_ == BP_BABCOCK)
        b2Axis.FlipXY();
    }
#   ifdef NASTRUCTDEBUG
    mprintf("\tDeviation from linear: %g deg.\n", z_deviation_from_linear * Constants::RADDEG);
    // Calculate parameters between axes.
    double Param[6];
    calculateParameters(b2Axis, b1Axis, 0, Param);
    mprintf("    Shear= %6.2f  Stretch= %6.2f  Stagger= %6.2f  Buck= %6.2f  Prop= %6.2f  Open= %6.2f\n",
            Param[0], Param[1], Param[2], Param[5]*Constants::RADDEG, Param[4]*Constants::RADDEG, Param[3]*Constants::RADDEG);
#   endif
    int NHB = CalcNumHB(base1, base2, n_wc_hb);
    BPmap::iterator entry = AddBasePair(b1idx, base1, b2idx, base2);
#   ifdef NASTRUCTDEBUG
    mprintf(", %i hbonds\n", NHB);
#   endif
    entry->second.nhb_ = NHB;
    entry->second.n_wc_hb_ = n_wc_hb;
    entry->second.isAnti_ = is_antiparallel;
    entry->second.isZ_ = is_z;
  } // END loop over pairs 

  return 0;

}

// Action_NAstruct::DetermineBasePairing()
/** Determine which bases are paired from the individual base axes and set up
  * entry in BasePairs_ if one not already present.
  */
int Action_NAstruct::DetermineBasePairing() {
  int n_wc_hb;
# ifdef NASTRUCTDEBUG  
  mprintf("\n=================== Setup Base Pairing ===================\n");
# endif
  // Loop over all possible pairs of bases 
  for (Barray::const_iterator base1 = Bases_.begin(); base1 != Bases_.end(); ++base1)
  {
    for (Barray::const_iterator base2 = base1 + 1; base2 != Bases_.end(); ++base2)
    {
      double dist2 = DIST2_NoImage(base1->Axis().Oxyz(), base2->Axis().Oxyz());
#     ifdef NASTRUCTDEBUG
      double axes_distance = sqrt(dist2);
      mprintf("\n  ----- Axes distance for %i:%s -- %i:%s is %f -----\n",
              base1->ResNum()+1, base1->ResName(), 
              base2->ResNum()+1, base2->ResName(), axes_distance);
#     endif
      if (dist2 < originCut2_) {
//#       ifdef NASTRUCTDEBUG
//        mprintf("  Axes distance for %i:%s -- %i:%s is %f\n",
//                base1->ResNum()+1, base1->ResName(), 
//                base2->ResNum()+1, base2->ResName(), sqrt(dist2));
//#       endif
#       ifdef NASTRUCTDEBUG
        // Glycosidic N-N distance
        if (base1->HasNXatom() && base2->HasNXatom()) {
          double n_n_dist2 = DIST2_NoImage(base1->NXxyz(), base2->NXxyz());
          mprintf("DEBUG: NX-NX distance= %f\n", sqrt(n_n_dist2));
        }
#       endif
        NA_Axis b1Axis = base1->Axis();
        NA_Axis b2Axis = base2->Axis();
        // Determine if base Z axis vectors are aligned with strand direction
        bool is_z = false;
        int b1_5to3 = axis_points_5p_to_3p( *base1 );
        int b2_5to3 = axis_points_5p_to_3p( *base2 );
        // TODO trap errors here
        // If antiparallel and both bases are aligned 3' to 5', may be ZDNA
        if (b1_5to3 == 0 && b2_5to3 == 0) {
#         ifdef NASTRUCTDEBUG
          mprintf("Both bases aligned 3' to 5', ZDNA\n");
#         endif
          b1Axis.FlipXZ();
          b2Axis.FlipXZ();
          is_z = true;
        }
        // Determine if base Z axis vectors point in same (theta <= 90) or
        // opposite (theta > 90) directions.
        bool is_antiparallel;
        double z_theta = b1Axis.Rz().Angle( b2Axis.Rz() );
        double z_deviation_from_linear;
        if (z_theta > Constants::PIOVER2) { // If theta(Z) > 90 deg.
#         ifdef NASTRUCTDEBUG
          mprintf("\t%s is anti-parallel to %s (%g deg)\n", base1->ResName(), base2->ResName(),
                  z_theta * Constants::RADDEG);
#         endif
          is_antiparallel = true;
          z_deviation_from_linear = Constants::PI - z_theta;
          // Antiparallel - flip Y and Z axes of complimentary base
          b2Axis.FlipYZ();
          
        } else {
#         ifdef NASTRUCTDEBUG
          mprintf("\t%s is parallel to %s (%g deg)\n", base1->ResName(), base2->ResName(),
                  z_theta * Constants::RADDEG);
#         endif
          is_antiparallel = false;
          z_deviation_from_linear = z_theta;
          // Parallel - no flip needed if 3dna.
          // If using Babcock convention, flip X and Y axes.
          if (bpConvention_ == BP_BABCOCK)
            b2Axis.FlipXY();
        }
#       ifdef NASTRUCTDEBUG
        mprintf("\tDeviation from linear: %g deg.\n", z_deviation_from_linear * Constants::RADDEG);
#       endif       
        // Calculate parameters between axes.
        double Param[6];
        //calculateParameters(base1->Axis(), base2->Axis(), 0, Param);
        calculateParameters(b2Axis, b1Axis, 0, Param);
#       ifdef NASTRUCTDEBUG
        mprintf("    Shear= %6.2f  Stretch= %6.2f  Stagger= %6.2f  Buck= %6.2f  Prop= %6.2f  Open= %6.2f\n",
                Param[0], Param[1], Param[2], Param[5]*Constants::RADDEG, Param[4]*Constants::RADDEG, Param[3]*Constants::RADDEG);
#       endif
        // Stagger (vertical separation) must be less than a cutoff.
        if ( fabs(Param[2]) < staggerCut_ ) {
          // Deviation from linear must be less than cutoff
          if (z_deviation_from_linear < z_angle_cut_) {
            int NHB = CalcNumHB(*base1, *base2, n_wc_hb);
            if (NHB > 0) {
              BPmap::iterator entry = AddBasePair(base1-Bases_.begin(), *base1,
                                                  base2-Bases_.begin(), *base2);
#             ifdef NASTRUCTDEBUG
              mprintf(", %i hbonds, axes distance= %g\n", NHB, axes_distance);
#             endif
              entry->second.nhb_ = NHB;
              entry->second.n_wc_hb_ = n_wc_hb;
              entry->second.isAnti_ = is_antiparallel;
              entry->second.isZ_ = is_z;
            } // END if # hydrogen bonds > 0
          } // END if Z angle < cut
        } // END if stagger < stagger cut
      } // END if base to base origin distance < cut
    } // END base2 loop
  } // END base1 loop

  return 0;
}

// -----------------------------------------------------------------------------
// AverageMatrices()
static Matrix_3x3 AverageMatrices(Matrix_3x3 const& RotatedR1, Matrix_3x3 const& RotatedR2) {
  Matrix_3x3 R;
  // Average R1 and R2 to get the middle frame
  for (int i = 0; i < 9; i++)
    R[i] = (RotatedR1[i] + RotatedR2[i]) / 2;
  // Normalize X, Y and Z vectors
  double r2 = sqrt( R[0]*R[0] + R[3]*R[3] + R[6]*R[6] );
  R[0] /= r2;
  R[3] /= r2;
  R[6] /= r2;
  r2 = sqrt( R[1]*R[1] + R[4]*R[4] + R[7]*R[7] );
  R[1] /= r2;
  R[4] /= r2;
  R[7] /= r2;
  r2 = sqrt( R[2]*R[2] + R[5]*R[5] + R[8]*R[8] );
  R[2] /= r2;
  R[5] /= r2;
  R[8] /= r2;
  return R;
}

// Action_NAstruct::calculateParameters()
/** Given two base axes, calculate translational and rotational parameters
  * between them. Store base pair axis info in BPaxis if given.
  */
int Action_NAstruct::calculateParameters(NA_Axis const& Axis1, NA_Axis const& Axis2, 
                                         NA_Axis* BPaxis, double *Param) 
{
# ifdef NASTRUCTDEBUG
  NA_Axis tempAxis;
  PDBfile paramfile;
  if (calcparam_)
    paramfile.OpenWrite("Param.pdb");
  Axis1.Oxyz().Print("O1");
  Axis1.Rot().Print("R1");
  Axis2.Oxyz().Print("O2");
  Axis2.Rot().Print("R2");
# endif
  // Hinge axis is cross product between Z1 and Z2
  Vec3 hingeAxis = Axis1.Rz().Cross( Axis2.Rz() );
# ifdef NASTRUCTDEBUG
  hingeAxis.Print("hinge");
# endif
  // Normalize hinge axis
  hingeAxis.Normalize();
# ifdef NASTRUCTDEBUG
  hingeAxis.Print("norm(hinge)");
# endif
  // Roll/Tilt is Angle between Z1 and Z2
  double rolltilt = Axis1.Rz().Angle( Axis2.Rz() );
# ifdef NASTRUCTDEBUG
  mprintf("\tAngle between Z1 and Z2= %f\n",rolltilt*Constants::RADDEG);
# endif
  // Calculate forward and backwards half rolltilt rotation around
  // hinge axis.
  Matrix_3x3 R;
  R.CalcRotationMatrix(hingeAxis, -0.5*rolltilt);
# ifdef NASTRUCTDEBUG
  R.Print("Rhalf");
# endif
  // Rotate R2 by -0.5 * rolltilt degrees around the hinge
  Matrix_3x3 RotatedR2 = R * Axis2.Rot();
  // Rotate R1 by 0.5 * rolltilt degrees around the hinge (inverse rotation)
  R.Transpose();
  Matrix_3x3 RotatedR1 = R * Axis1.Rot();
# ifdef NASTRUCTDEBUG
  // Print rotated R1 and R2
  RotatedR1.Print("Rotated R1");
  RotatedR2.Print("Rotated R2");
  if (calcparam_) {
    tempAxis.StoreRotMatrix(RotatedR1, Axis1.Oxyz()); 
    WriteAxes(paramfile, 1, "R1'", tempAxis);
    tempAxis.StoreRotMatrix(RotatedR2, Axis2.Oxyz());
    WriteAxes(paramfile, 2, "R2'", tempAxis);
  }
# endif
  // Average rotated R1 and R2 to get the middle frame
  Matrix_3x3 Rm = AverageMatrices(RotatedR1, RotatedR2);
  // Take average of origins
  Vec3 OM = (Axis1.Oxyz() + Axis2.Oxyz()) / 2.0;
# ifdef NASTRUCTDEBUG
  OM.Print("Origin Mean");
  // Print Rm and hinge axis
  Rm.Print("Rm");
  if (calcparam_) {
    // Use R to store hinge axis in Z
    R.Zero();
    R[2] = hingeAxis[0]; 
    R[5] = hingeAxis[1]; 
    R[8] = hingeAxis[2];
    tempAxis.StoreRotMatrix(R, OM);
    WriteAxes(paramfile, 3, "Hng", tempAxis);
    // Store middle frame
    tempAxis.StoreRotMatrix(Rm, OM);
    WriteAxes(paramfile, 4, "Rm", tempAxis);
  }
# endif

  // If BPaxis is not null, store Rm and OM as BP axis.
  if (BPaxis != 0)  
    BPaxis->StoreRotMatrix(Rm, OM);

  // Shift Slide Rise / Shear Stretch Stagger
  OM = Axis2.Oxyz() - Axis1.Oxyz();
  // Since this is really vector times matrix, use matrix transpose times vec
  Vec3 Vec = Rm.TransposeMult( OM );
# ifdef NASTRUCTDEBUG
  OM.Print("O21");
  Vec.Print("Vec");
# endif
  Param[0] = Vec[0];
  Param[1] = Vec[1];
  Param[2] = Vec[2];

  // Set Z1 to Z from middle frame
  Vec3 Z1 = Rm.Col3();
  // Twist / Opening
  // Angle between rotated Y1 and rotated Y2
  // Sign of twistopen related to (Y1'xY2') dot Z of middle frame
  Vec3 Y1 = RotatedR1.Col2();
  Vec3 Y2 = RotatedR2.Col2();
  double twistopen = Y1.SignedAngle(Y2, Z1);
# ifdef NASTRUCTDEBUG
  mprintf("\tFinal Twist/Opening is %10.4f\n",twistopen*Constants::RADDEG);
# endif
  Param[3] = twistopen;

  // Phase angle
  // Angle between hinge axis and middle frame Y axis
  // Sign of phi related to (hingeAxis x Ym) dot Z of middle frame
  Y1 = Rm.Col2();
  double phi = hingeAxis.SignedAngle(Y1, Z1);
  double sinphi = sin( phi );
  double cosphi = cos( phi );
# ifdef NASTRUCTDEBUG
  mprintf("\tPhase angle is %f, sinphi is %f, cosphi is %f\n",phi*Constants::RADDEG,sinphi,cosphi);
# endif

  // Roll / Propeller
  double rollprop = rolltilt * cosphi;
  Param[4] = rollprop;

  // Tilt / Buckle
  double tiltbuck = rolltilt * sinphi;
  Param[5] = tiltbuck;

# ifdef NASTRUCTDEBUG
  mprintf("\tRoll/Propeller %10.4f\n",rollprop*Constants::RADDEG);
  mprintf("\tTilt/Buckle %10.4f\n",tiltbuck*Constants::RADDEG);
  if (calcparam_) calcparam_=false;
# endif
  return 0;
}

// Action_NAstruct::helicalParameters()
int Action_NAstruct::helicalParameters(NA_Axis const& Axis1, NA_Axis const& Axis2, double *Param) 
{
  // O1 = X2 - X1
  Vec3 O1 = Axis2.Rx() - Axis1.Rx();
  // O2 = Y2 - Y1
  Vec3 O2 = Axis2.Ry() - Axis1.Ry();
  // Local helical axis: (X2-X1) x (Y2-Y1)
  Vec3 helicalAxis = O1.Cross( O2 );
# ifdef NASTRUCTDEBUG
  O1.Print("X2 - X1");
  O2.Print("Y2 - Y1");
  helicalAxis.Print("(X2-X1) x (Y2-Y1)");
# endif
  helicalAxis.Normalize( );
# ifdef NASTRUCTDEBUG
  helicalAxis.Print("NORM[(X2-X1)x(Y2-Y1)]");
# endif

  // Tip/inclination is angle between helical axis and z1
  double tipinc = helicalAxis.Angle( Axis1.Rz() );
  // Hinge axis is normalized cross product of helical axis to z1
  Vec3 hingeAxis = helicalAxis.Cross( Axis1.Rz() );
  hingeAxis.Normalize();
  // Rotate R1 around hinge axis by -tipinc
  Matrix_3x3 R;
  R.CalcRotationMatrix(hingeAxis, -tipinc);
  Matrix_3x3 RotatedR1 = R * Axis1.Rot();
# ifdef NASTRUCTDEBUG
  mprintf("\tTip/Inclination: %f\n",tipinc*Constants::RADDEG);
  hingeAxis.Print("Hinge axis 1");
  RotatedR1.Print("Rotated R1");
# endif

  // Tip/inclination should be same for z2
  //mprintf("\tTipCheck= %f\n",dot_product_angle(helicalAxis, Z2)*Constants::RADDEG);
  // Hinge axis (Vec) is normalized cross product from h to z2
  Vec3 Vec = helicalAxis.Cross( Axis2.Rz() );
  Vec.Normalize();
  // Rotate R2 around hinge axis by -tipinc
  R.CalcRotationMatrix(Vec, -tipinc); 
  Matrix_3x3 RotatedR2 = R * Axis2.Rot();
# ifdef NASTRUCTDEBUG
  Vec.Print("Hinge axis 2");
  RotatedR2.Print("Rotated R2");
# endif

  // Average Rotated R1 and R2 to get middle helical frame
  R = AverageMatrices(RotatedR1, RotatedR2);

  // Helical twist is angle from Rotated Y1 to Rotated Y2
  // Sign is given by (Y1'xY2' dot helicalAxis)
  Vec3 Y1 = RotatedR1.Col2();
  Vec3 Y2 = RotatedR2.Col2();
  double Twist = Y1.SignedAngle(Y2, helicalAxis);
  Param[5] = Twist;

  // Calc Vec = O2 - O1
  Vec = Axis2.Oxyz() - Axis1.Oxyz();
  // Project (O2-O1) onto helical axis
  double Rise = Vec * helicalAxis;
  Param[2] = Rise;
# ifdef NASTRUCTDEBUG
  R.Print("Hm");
  mprintf("\tTwist is %f\n",Twist*Constants::RADDEG);
  mprintf("\tRise is %f\n",Rise);
# endif

  // Phase angle is angle from hinge Axis 1 to RotatedR1 Y
  // Sign is given by (hingeAxis x Y1') dot helicalAxis
  double phase = hingeAxis.SignedAngle(Y1, helicalAxis);

  // Tip is tipinc * cos( phase )
  double Tip = tipinc * cos( phase );
  Param[4] = Tip;
  // Inclination is tipinc * sin( phase )
  double Inc = tipinc * sin( phase );
  Param[3] = Inc;
# ifdef NASTRUCTDEBUG
  mprintf("\tPhase angle is %f\n",phase*Constants::RADDEG);
  mprintf("\tTip is %f\n",Tip*Constants::RADDEG);
  mprintf("\tInclination is %f\n",Inc*Constants::RADDEG);
# endif

  Vec3 Z1 = helicalAxis * Rise;
  // Calc vector AB (store in O1) 
  // Vec contains O2-O1
  O1 = Vec - Z1; 

  // Calc vector AD
  double AD_angle = Constants::PIOVER2 - (0.5 * Twist);
  // rotation of AD_angle around helicalAxis
  // NOTE: Assuming we dont need RotatedR2 anymore
  RotatedR2.CalcRotationMatrix(helicalAxis, AD_angle);
  // rotate AB, = AD (store in O2)
  O2 = RotatedR2 * O1;
  O2.Normalize();
# ifdef NASTRUCTDEBUG
  O1.Print("AB");
  mprintf("\tAD_angle is %f\n",AD_angle*Constants::RADDEG);
  O2.Print("AD");
# endif

  // Calc magnitude of AD; 0.5 * |AB| / sin( 0.5 * Twist )
  double AD_mag = (0.5 * sqrt(O1.Magnitude2())) / sin( 0.5 * Twist );

  // Calc origin of local helical frame for BP 1
  // O1 = Origin1 + (AD_mag * AD)
  O1 = Axis1.Oxyz() + (O2 * AD_mag); 

  // Calc origin of local helical frame for BP 2
  // O2 = O1 + (Rise * helicalAxis)
  // Z1 contains helicalAxis * Rise
  O2 = O1 + Z1;

  // Calculate origin of middle helical frame, store in Vec 
  Vec = (O2 + O1) / 2.0;
# ifdef NASTRUCTDEBUG
  mprintf("\t|AD| = %f\n",AD_mag);
  O1.Print("o1_h");
  O2.Print("o2_h");
  Vec.Print("Om_h");
# endif

  // Calc vector from O1 to Origin1
  Vec = Axis1.Oxyz() - O1;

  // X-disp is projection of vector from O1 to Origin1 onto 
  // X axis of RotatedR1.
  // TODO: X1 could just be O1
  Vec3 X1 = RotatedR1.Col1();
  double X_disp = Vec * X1;
  Param[0] = X_disp;

  // Y-disp is projection of vector from O1 to Origin1 onto 
  // Y axis of RotatedR1.
  X1 = RotatedR1.Col2();
  double Y_disp = Vec * X1;
  Param[1] = Y_disp;
# ifdef NASTRUCTDEBUG
  mprintf("\tX-displacement= %f\n",X_disp);
  mprintf("\tY-displacement= %f\n",Y_disp);
# endif

  return 0;
}

/** Get index of base in Bases_ N steps away from base specified by idx.
  * Positive value is towards 3', negative value is towards 5'.
  */
int Action_NAstruct::GetBaseIdxStep(int idx, int Nsteps) const {
  if (Nsteps == 0 || idx == -1) return idx;
  NA_Base const& base = Bases_[idx];
  if (Nsteps > 0)
    return GetBaseIdxStep( base.C3resIdx(), Nsteps - 1);
  else // Nsteps < 0
    return GetBaseIdxStep( base.C5resIdx(), Nsteps + 1);
}

/** For bases in a single strand, get the values of buckle, propeller twist,
  * opening, shear, stretch, and stagger.
  */
int Action_NAstruct::DetermineStrandParameters(int frameNum) {
  double Param[6];
  NA_Axis commonAxis;
  // Loop over strand pairs
  for (Smap::const_iterator it = StrandPairs_.begin(); it != StrandPairs_.end(); ++it)
  {
    int b1idx = it->first.first;
    int b2idx = it->first.second;
    Stype const& SP = it->second;
    //mprintf("DEBUG: Strand %u, bases %u to %u\n", SP.strandidx_, b1idx, b2idx);

    // Get bases
    NA_Base const& base1 = Bases_[b1idx];
    NA_Base const& base2 = Bases_[b2idx];
    // Calc parameters between bases in the strand
    calculateParameters(base2.Axis(), base1.Axis(), &commonAxis, Param);
    // Store data
    Param[3] *= Constants::RADDEG;
    Param[4] *= Constants::RADDEG;
    Param[5] *= Constants::RADDEG;
    //mprintf("DBG: BP %i # hbonds = %i\n", nbasepair+1, NumberOfHbonds_[nbasepair]);
    // Convert everything to float to save space
    float shear = (float)Param[0];
    float stretch = (float)Param[1];
    float stagger = (float)Param[2];
    float opening = (float)Param[3];
    float prop = (float)Param[4];
    float buckle = (float)Param[5];
    // Add to DataSets
    //mprintf("DEBUG:\tShear= %f  stretch= %f  stagger= %f\n", shear, stretch, stagger);
    //mprintf("DEBUG:\tOpeni= %f  propell= %f  buckle_= %f\n", opening, prop, buckle);
    SP.dx_->Add(frameNum, &shear);
    SP.dy_->Add(frameNum, &stretch);
    SP.dz_->Add(frameNum, &stagger);
    SP.rx_->Add(frameNum, &opening);
    SP.ry_->Add(frameNum, &prop);
    SP.rz_->Add(frameNum, &buckle);
  }
  return 0;
}

/** \return 1 if base Z axis is pointing in 5' to 3' direction.
  * \return 0 if base Z axis is pointing 3' to 5'
  * \return -1 if an error occurs.
  */
int Action_NAstruct::axis_points_5p_to_3p(NA_Base const& base1) const {
  if (base1.HasC1atom()) {
    // Ensure base Z vector points 5' to 3'
    int c3residx = base1.C3resIdx();
    int c5residx = base1.C5resIdx();
    if (c3residx > -1 || c5residx > -1) {
      const double *c5res_c1xyz = 0;
      const double *c3res_c1xyz = 0;
      if (c3residx == -1)
        c3res_c1xyz = base1.C1xyz();
      else if (Bases_[c3residx].HasC1atom())
        c3res_c1xyz = Bases_[c3residx].C1xyz();
      if (c5residx == -1)
        c5res_c1xyz = base1.C1xyz();
      else if (Bases_[c5residx].HasC1atom())
        c5res_c1xyz = Bases_[c5residx].C1xyz();
      if (c5res_c1xyz == 0 || c3res_c1xyz == 0) {
        mprinterr("Error: 5' res and/ord 3' res missing C1 atom coords.\n");
        return -1;
      }
#     ifdef NASTRUCTDEBUG
      mprintf("DEBUG: Res %i  c5res_c1xyz = %f %f %f  c3res_c1xyz = %f %f %f\n",
              base1.ResNum()+1,
              c5res_c1xyz[0], c5res_c1xyz[1], c5res_c1xyz[2],
              c3res_c1xyz[0], c3res_c1xyz[1], c3res_c1xyz[2]);
#     endif
      Vec3 strand_vec( c3res_c1xyz[0] - c5res_c1xyz[0],
                       c3res_c1xyz[1] - c5res_c1xyz[1],
                       c3res_c1xyz[2] - c5res_c1xyz[2] );
      strand_vec.Normalize();
      double s_angle = base1.Axis().Rz().Angle( strand_vec );
#     ifdef NASTRUCTDEBUG
      strand_vec.Print("strand vector");
      base1.Axis().Rz().Print("Axis Z");
      mprintf("DEBUG: Angle between strand and Axis Z = %f\n", s_angle * Constants::RADDEG);
#     endif
      if (s_angle > Constants::PIOVER2) {
        // Z has flipped, likely due to rotation around chi.
#       ifdef NASTRUCTDEBUG
        mprintf("DEBUG: Res %i Z axis is aligned 3' to 5'\n", base1.ResNum()+1);
#       endif
        return 0;
      } else {
        // Z is aligned 5' to 3'
#       ifdef NASTRUCTDEBUG
        mprintf("DEBUG: Res %i Z axis is aligned 5' to 3'\n", base1.ResNum()+1);
#       endif
        return 1;
      }
    }
  }
  // If we are here, either not enough coords or not enough res. Assume aligned. TODO OK?
  return 1;
}

// Action_NAstruct::DeterminePairParameters()
/** For each base pair, get the values of buckle, propeller twist,
  * opening, shear, stretch, and stagger. Also determine the origin and 
  * rotation matrix for each base pair reference frame.
  */
int Action_NAstruct::DeterminePairParameters(int frameNum) {
  double Param[6];
# ifdef NASTRUCTDEBUG
  mprintf("\n=================== Determine BP Parameters ===================\n");
# endif
  // NOTE: iterator cannot be const because bpaxis_ needs to be updated
  for (BPmap::iterator it = BasePairs_.begin(); it != BasePairs_.end(); ++it)
  {
    if (it->second.nhb_ < 1 && skipIfNoHB_) continue;
    BPtype& BP = it->second;
    int b1 = BP.base1idx_;
    int b2 = BP.base2idx_;
    NA_Base& base1 = Bases_[b1];
    NA_Base& base2 = Bases_[b2]; //TODO copy? 

#   ifdef NASTRUCTDEBUG
    mprintf("BasePair %i:%s to %i:%s", b1+1, base1.ResName(), b2+1, base2.ResName());
    if (BP.isAnti_)
      mprintf(" Anti-parallel,");
    else
      mprintf(" Parallel,");
    if (BP.isZ_)
      mprintf(" aligned 3' to 5' (Z).\n");
    else
      mprintf(" aligned 5' to 3' (A/B).\n");
    base2.Axis().Rot().Print("Original base2 axis");
#   endif
    // Check Z-DNA - need to flip XZ
    if (BP.isZ_) {
      base1.Axis().FlipXZ();
      base2.Axis().FlipXZ();
    }
    // Check Antiparallel / Parallel
    // Flip YZ (rotate around X) for antiparallel
    // Flip XY (rotate around Z) for parallel (Babcock convention only)
    if (BP.isAnti_)
      base2.Axis().FlipYZ();
    else if (bpConvention_ == BP_BABCOCK)
      base2.Axis().FlipXY();
    if (grooveCalcType_ == PP_OO) {
      // Calc direct P--P distance
      float dPtoP = 0.0;
      //mprintf("\tDEBUG: %i %i:", BaseAxes[base1].Pidx(), BaseAxes[base2].Pidx() );
      if ( base1.HasPatom() && base2.HasPatom() ) {
        double DP = DIST2_NoImage( base1.Pxyz(), base2.Pxyz() );
        //mprintf(" %i to %i P--P D= %f", BaseAxes[base1].ResNum()+1, BaseAxes[base2].ResNum()+1,
        //        sqrt(dPtoP) );
        DP = sqrt(DP);
        dPtoP = (float)DP;
      }
      //mprintf("\n");
      float dOtoO = 0.0;
      //mprintf("\tDEBUG: %i %i:", BaseAxes[base1].O4idx(), BaseAxes[base2].O4idx() );
      if ( base1.HasO4atom() && base2.HasO4atom() ) {
        double DO4 = DIST2_NoImage( base1.O4xyz(), base2.O4xyz() );
        //mprintf(" %i to %i O4'--O4' D= %f", BaseAxes[base1].ResNum()+1, BaseAxes[base2].ResNum()+1,
        //        sqrt(dOtoO) );
        DO4 = sqrt(DO4);
        dOtoO = (float)DO4;
      }
      BP.major_->Add(frameNum, &dPtoP);
      BP.minor_->Add(frameNum, &dOtoO);
    }
    //mprintf("\n");
    // Calc BP parameters, set up basepair axes
    //calculateParameters(BaseAxes[base1],BaseAxes[base2],&BasePairAxes[nbasepair],Param);
    calculateParameters(base2.Axis(), base1.Axis(), &(BP.bpaxis_), Param);
    // Store data
    Param[3] *= Constants::RADDEG;
    Param[4] *= Constants::RADDEG;
    Param[5] *= Constants::RADDEG;
    //mprintf("DBG: BP %i # hbonds = %i\n", nbasepair+1, NumberOfHbonds_[nbasepair]);
    // Convert everything to float to save space
    float shear = (float)Param[0];
    float stretch = (float)Param[1];
    float stagger = (float)Param[2];
    float opening = (float)Param[3];
    float prop = (float)Param[4];
    float buckle = (float)Param[5];
    // Add to DataSets
    BP.shear_->Add(frameNum, &shear);
    BP.stretch_->Add(frameNum, &stretch);
    BP.stagger_->Add(frameNum, &stagger);
    BP.opening_->Add(frameNum, &opening);
    BP.prop_->Add(frameNum, &prop);
    BP.buckle_->Add(frameNum, &buckle);
    if (wc_hb_only_)
      BP.hbonds_->Add(frameNum, &(BP.n_wc_hb_));
    else
      BP.hbonds_->Add(frameNum, &(BP.nhb_));
    static const int ONE = 1;
    if (BP.nhb_ > 0)
      BP.isBP_->Add(frameNum, &ONE);
#   ifdef NASTRUCTDEBUG
    double bp_axes_vec[6]; // Hold base pair axis Z vec and origin
    bp_axes_vec[0] = BP.bpaxis_.Rz()[0];
    bp_axes_vec[1] = BP.bpaxis_.Rz()[1];
    bp_axes_vec[2] = BP.bpaxis_.Rz()[2];
    bp_axes_vec[3] = BP.bpaxis_.Oxyz()[0];
    bp_axes_vec[4] = BP.bpaxis_.Oxyz()[1];
    bp_axes_vec[5] = BP.bpaxis_.Oxyz()[2];
    BP.axes_nxyz_->Add(frameNum, bp_axes_vec);
#   endif
  }
  // Calculate base parameters.
  for (Barray::iterator base = Bases_.begin(); base != Bases_.end(); ++base)
    base->CalcPucker( frameNum, puckerMethod_ );

  return 0;
}

// Action_NAstruct::DetermineStepParameters() 
/** Determine base pair steps and values of Tilt, Roll, Twist, Shift,
  * Slide, and Rise.
  */
int Action_NAstruct::DetermineStepParameters(int frameNum) {
  double Param[6];
# ifdef NASTRUCTDEBUG
  mprintf("\n=================== Determine BPstep Parameters ===================\n");
# endif
  if (BasePairs_.size() < 2) return 0;
  // Determine which base pairs are in step configuration.
  // Step will be defined as:
  //   base1 -- base2
  //     |        |
  //   base3 -- base4
  for (BPmap::const_iterator bp1 = BasePairs_.begin(); bp1 != BasePairs_.end(); ++bp1) {
    BPtype const& BP1 = bp1->second;
    if (BP1.nhb_ < 1 && skipIfNoHB_) continue; // Base pair not valid this frame.
    NA_Base const& base1 = Bases_[BP1.base1idx_];
    NA_Base const& base2 = Bases_[BP1.base2idx_];
    
    // Find step
    int idx1 = base1.C3resIdx();
    int idx2;
    if (BP1.isAnti_)
      idx2 = base2.C5resIdx();
    else
      idx2 = base2.C3resIdx();
    if (idx1 != -1 && idx2 != -1) {
      Rpair respair(Bases_[idx1].ResNum(), Bases_[idx2].ResNum());
      BPmap::const_iterator bp2 = BasePairs_.find( respair );
      if (bp2 != BasePairs_.end() && (bp2->second.nhb_ > 0 || !skipIfNoHB_)) {
        BPtype const& BP2 = bp2->second;
        NA_Base const& base3 = Bases_[BP2.base1idx_];
        NA_Base const& base4 = Bases_[BP2.base2idx_];
#       ifdef NASTRUCTDEBUG
        mprintf("  BP step (%s--%s)-(%s--%s)\n",
                base1.BaseName().c_str(), base2.BaseName().c_str(),
                base3.BaseName().c_str(), base4.BaseName().c_str());
#       endif
        // NOTE: Unlike base pairs which are indexed by residue numbers, base
        //       pair steps are indexed by base pair indices.
        Rpair steppair(BP1.bpidx_, BP2.bpidx_);
        // Base pair step. Try to find existing base pair step.
        StepMap::iterator entry = Steps_.lower_bound( steppair );
        if (entry == Steps_.end() || entry->first != steppair) {
          // New base pair step
          StepType BS;
          MetaData md = NewStepType(BS, BP1.base1idx_, BP1.base2idx_,
                                        BP2.base1idx_, BP2.base2idx_, Steps_.size()+1);
          // H-C groove width calc setup
          if (grooveCalcType_ == HASSAN_CALLADINE) {
            // Major groove
            BS.p_p2_ = -1;
            BS.P_m2_ = GetBaseIdxStep( BS.b3idx_, -2 );
            if (BP1.isAnti_)
              // For anti, want base3 - 2 (to 5'), base2 - 2 (to 5')
              BS.p_p2_ = GetBaseIdxStep( BS.b2idx_, -2 );
            else
              // For para, want base3 - 2 (to 5'), base4 + 2 (to 3')
              BS.p_p2_ = GetBaseIdxStep( BS.b4idx_, +2 );
            if (BS.P_m2_ != -1 && BS.p_p2_ != -1 && 
                Bases_[BS.P_m2_].HasPatom() && Bases_[BS.p_p2_].HasPatom())
            {
              md.SetAspect("major");
              BS.majGroove_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
              //mprintf("DEBUG: Groove '%s' P-2= %i, p+2= %i\n",
              //        BS.majGroove_->legend(), BS.P_m2_+1, BS.p_p2_+1);
            }
            // Minor groove
            BS.p_m1_ = -1;
            BS.p_m2_ = -1;
            BS.P_p1_ = GetBaseIdxStep( BS.b3idx_, +1 );
            BS.P_p2_ = GetBaseIdxStep( BS.b3idx_, +2 );
            if (BP1.isAnti_) {
              BS.p_m1_ = GetBaseIdxStep( BS.b2idx_, +1 );
              BS.p_m2_ = GetBaseIdxStep( BS.b2idx_, +2 );
            } else {
              BS.p_m1_ = GetBaseIdxStep( BS.b4idx_, -1 );
              BS.p_m2_ = GetBaseIdxStep( BS.b4idx_, -2 );
            }
            if (BS.P_p1_ != -1 && BS.P_p2_ != -1 && BS.p_m1_ != -1 && BS.p_m2_ != -1 &&
                Bases_[BS.P_p1_].HasPatom() && Bases_[BS.P_p2_].HasPatom() &&
                Bases_[BS.p_m1_].HasPatom() && Bases_[BS.p_m2_].HasPatom())
            {
              md.SetAspect("minor");
              BS.minGroove_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
              //mprintf("DEBUG: Groove '%s' P+1= %i, P+2= %i, p-1= %i, p-2= %i\n",
              //        BS.minGroove_->legend(), BS.P_p1_+1, BS.P_p2_+1, BS.p_m1_+1, BS.p_m2_+1);
            }
          }
          entry = Steps_.insert( entry, std::pair<Rpair, StepType>(steppair, BS) );
#         ifdef NASTRUCTDEBUG
          mprintf("  New base pair step: %s\n", md.Legend().c_str());
#         endif
        }
        StepType& currentStep = entry->second;
        // Calc step parameters
        calculateParameters(BP1.bpaxis_, BP2.bpaxis_, &(currentStep.stepaxis_), Param);
        NA_Axis& midFrame = currentStep.stepaxis_;
        // Calculate zP: difference in step phosphate atoms along the Z axis
        // of the step middle frame.
        float Zp = 0.0;
        NA_Base const* pbase1 = &base3;
        NA_Base const* pbase2 = 0;
        if (BP1.isAnti_)
          pbase2 = &base2;
        else
          pbase2 = &base4;
        if (pbase1->HasPatom() && pbase2->HasPatom()) { // TODO warn if atoms missing?
          // Rotate the coordinates of the base phosphate atoms so they are
          // in the step middle frame.
          Vec3 p1xyz = midFrame.Rot().TransposeMult( Vec3(pbase1->Pxyz()) - midFrame.Oxyz() );
          Vec3 p2xyz = midFrame.Rot().TransposeMult( Vec3(pbase2->Pxyz()) - midFrame.Oxyz() );
          // If anti-parallel, ensure axes on second strand are properly flipped
          if (BP1.isAnti_) {
            // FlipYZ
            p2xyz[1] = -p2xyz[1];
            p2xyz[2] = -p2xyz[2];
          } else if (bpConvention_ == BP_BABCOCK) {
            // FlipXY
            p2xyz[0] = -p2xyz[0];
            p2xyz[1] = -p2xyz[1];
          }
          Vec3 xyzP = (p1xyz + p2xyz) / 2.0;
#         ifdef NASTRUCTDEBUG
          mprintf("ZPCALC: SI '%3s' P %6.2f %6.2f %6.2f  |  SII '%3s' P %6.2f %6.2f %6.2f  |  %6.2f %6.2f %6.2f\n",
                  pbase1->BaseName().c_str(), p1xyz[0], p1xyz[1], p1xyz[2],
                  pbase2->BaseName().c_str(), p2xyz[0], p2xyz[1], p2xyz[2],
                  xyzP[0], xyzP[1], xyzP[2]);
#         endif
          Zp = (float)xyzP[2];
        }
        currentStep.Zp_->Add(frameNum, &Zp);
        // TEST: Calculate major groove ----------
        if (grooveCalcType_ == HASSAN_CALLADINE) {
          if (currentStep.majGroove_ != 0) {
            double MGW = DIST2_NoImage( Bases_[currentStep.P_m2_].Pxyz(),
                                        Bases_[currentStep.p_p2_].Pxyz() );
            //mprintf("DEBUG:\t\tMajorGroove= %4.1f\n", sqrt(MGW));
            float fval = (float)sqrt( MGW );
            currentStep.majGroove_->Add(frameNum, &fval);
          }
          if (currentStep.minGroove_ != 0) {
            double d1 = sqrt(DIST2_NoImage( Bases_[currentStep.P_p1_].Pxyz(),
                                            Bases_[currentStep.p_m2_].Pxyz() ));
            double d2 = sqrt(DIST2_NoImage( Bases_[currentStep.P_p2_].Pxyz(),
                                            Bases_[currentStep.p_m1_].Pxyz() ));
            double mGW = 0.5 * (d1 + d2);
            //mprintf("DEBUG:\t\tMinorGroove= %4.1f\n", mGW);
            float fval = (float)mGW;
            currentStep.minGroove_->Add(frameNum, &fval);
          }
        }
        // ---------------------------------------
        // Store data
        Param[3] *= Constants::RADDEG;
        Param[4] *= Constants::RADDEG;
        Param[5] *= Constants::RADDEG;
        // Convert everything to float to save space
        float shift = (float)Param[0];
        float slide = (float)Param[1];
        float rise = (float)Param[2];
        float twist = (float)Param[3];
        float roll = (float)Param[4];
        float tilt = (float)Param[5];
        currentStep.shift_->Add(frameNum, &shift);
        currentStep.slide_->Add(frameNum, &slide);
        currentStep.rise_->Add(frameNum, &rise);
        currentStep.twist_->Add(frameNum, &twist);
        currentStep.roll_->Add(frameNum, &roll);
        currentStep.tilt_->Add(frameNum, &tilt);
        // Calc helical parameters
        helicalParameters(BP1.bpaxis_, BP2.bpaxis_, Param);
        Param[3] *= Constants::RADDEG;
        Param[4] *= Constants::RADDEG;
        Param[5] *= Constants::RADDEG;
        // Convert to float
        float xdisp = (float)Param[0];
        float ydisp = (float)Param[1];
        float hrise = (float)Param[2];
        float incl = (float)Param[3];
        float tip = (float)Param[4];
        float htwist = (float)Param[5];
        currentStep.xdisp_->Add(frameNum, &xdisp);
        currentStep.ydisp_->Add(frameNum, &ydisp);
        currentStep.hrise_->Add(frameNum, &hrise);
        currentStep.incl_->Add(frameNum, &incl);
        currentStep.tip_->Add(frameNum, &tip);
        currentStep.htwist_->Add(frameNum, &htwist);
      } // END second base pair found
    } // END second base pair valid
  } // END loop over base pairs 
  return 0;
}
// ----------------------------------------------------------------------------

/** Starting at atom, try to travel phosphate backbone to atom in next res.
  * \return Atom # of atom in next residue or -1 if no next residue.
  */
int Action_NAstruct::TravelBackbone( Topology const& top, int atom, std::vector<int>& Visited ) {
  Visited[atom] = 1;
  for (Atom::bond_iterator bndat = top[atom].bondbegin();
                           bndat != top[atom].bondend(); ++bndat)
  {
    // If not all residues are selected it is possible bndat is outside
    // the scope of the Visited array.
    if ( *bndat < (int)Visited.size() && Visited[*bndat] == 0) {
      if ( top[*bndat].Element() == Atom::CARBON || 
           top[*bndat].Element() == Atom::HYDROGEN )
        Visited[*bndat] = 1;
      else if ( top[*bndat].ResNum() != top[atom].ResNum() )
        return *bndat;
      else {
        int tatom = TravelBackbone( top, *bndat, Visited );
        if (tatom != -1) return tatom;
      }
    }
  }
  return -1;
}

/** \return index in bases of 3' termini */
int Action_NAstruct::follow_base_to_3prime(Barray& Bases, unsigned int bidx, std::vector<bool>& base_is_used, int strandNum)
{
  if (base_is_used[bidx]) {
    mprinterr("Internal Error: Base %u is used already.\n", bidx);
    return -1;
  }
  base_is_used[bidx] = true;
  Bases[bidx].SetStrandNum( strandNum );
  // If no more bases in the 3' direction, this is the 3' terminal
  if (Bases[bidx].C3resIdx() == -1)
    return bidx;
  return follow_base_to_3prime(Bases, Bases[bidx].C3resIdx(), base_is_used, strandNum);
}

// Action_NAstruct::Setup()
/** Determine which bases to use in calculation and set up their reference frames.
  */
Action::RetType Action_NAstruct::Setup(ActionSetup& setup) {
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  Range actualRange;
  if (resRange_.Empty()) 
    actualRange = setup.Top().SoluteResidues();
  else 
    actualRange = resRange_;
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprintf("Warning: No residues specified for %s\n",setup.Top().c_str());
    return Action::SKIP;
  }
  if (dataname_.empty())
    dataname_ = masterDSL_->GenerateDefaultName("NA");
  // Try to fill in any custom residue maps.
  if (!nameToRef_.empty()) {
    std::vector<RefMapType::iterator> nToErase;
    for (RefMapType::iterator it = nameToRef_.begin(); it != nameToRef_.end(); ++it)
    {
      mprintf("\tAttempting to create custom map for residue %s\n", it->first.c_str());
      // Search for name in actual range.
      NameType tgtResName( it->first );
      for (Range::const_iterator resnum = actualRange.begin();
                                 resnum != actualRange.end(); ++resnum)
      {
        if (setup.Top().Res(*resnum).Name() == tgtResName) {
          mprintf("\t  Using residue %s\n", setup.Top().TruncResNameNum(*resnum).c_str());
          if (refBases_.AddCustomBase(tgtResName, setup.Top(), *resnum, it->second))
            return Action::ERR;
          else {
            // Base successfully set up. Exit this loop.
            // This mapping can be removed since it has been satisfied.
            nToErase.push_back( it );
            break;
          }
        }
      } // END loop over residues in topology
    } // END loop over residues in custom map
    // Remove names from the custom map if they were set up
    if (!nToErase.empty()) {
      for (std::vector<RefMapType::iterator>::const_iterator it = nToErase.begin();
                                                             it != nToErase.end(); ++it)
        nameToRef_.erase( *it );
    }
  }
  // DEBUG - print all residues
  //if (debug>0)
  //  actualRange.PrintRange("    NAstruct: NA res:",1);
  bool firstTimeSetup = Bases_.empty();
  unsigned int idx = 0;
  // Set up NA_base for each selected NA residue 
  for (Range::const_iterator resnum = actualRange.begin();
                             resnum != actualRange.end(); ++resnum)
  {
#   ifdef NASTRUCTDEBUG
    mprintf(" ----- Setting up %i:%s -----\n", *resnum+1, setup.Top().Res(*resnum).c_str());
#   endif
    // Set up ref coords for this base
    NA_Base currentBase;
    NA_Reference::RetType err = refBases_.SetupBaseRef( currentBase, setup.Top(), *resnum,
                                                        *masterDSL_, dataname_ );
    if (err == NA_Reference::NOT_FOUND) {
      // Residue not recognized. Print a warning.
      mprintf("Warning: Residue %i:%s not recognized as a nucleic acid residue.\n",
              *resnum+1, setup.Top().Res(*resnum).c_str());
      mprintf("Warning: For non-standard names, use the 'resmap' keyword.\n");
      continue;
    } else if (err == NA_Reference::BASE_ERROR) {
      mprinterr("Error: Could not set up residue %s for NA structure analysis.\n",
                setup.Top().TruncResNameNum(*resnum).c_str());
      return Action::ERR;
    }
    if (firstTimeSetup) {
      Bases_.push_back( currentBase );
      // Determine the largest residue for setting up frames for RMS fit later.
      maxResSize_ = std::max( maxResSize_, currentBase.InputFitMask().Nselected() );
      if (debug_>1) {
        mprintf("\tNAstruct: Res %i:%s ", *resnum+1, setup.Top().Res(*resnum).c_str());
        Bases_.back().PrintAtomNames();
        Bases_.back().InputFitMask().PrintMaskAtoms("InpMask");
        Bases_.back().RefFitMask().PrintMaskAtoms("RefMask");
      }
    } else {
      // Ensure base type has not changed. //TODO: Re-set up reference? Check # atoms etc?
      if (currentBase.Type() != Bases_[idx].Type()) {
        mprinterr("Error: Residue %s base has changed from %s to %s\n",
                  setup.Top().TruncResNameNum(*resnum).c_str(), Bases_[idx].BaseName().c_str(),
                  currentBase.BaseName().c_str());
        return Action::ERR;
      }
    }
    idx++;
  } // End Loop over NA residues
  mprintf("\tSet up %zu bases.\n", Bases_.size());
  if (Bases_.empty()) {
    mprintf("Warning: No NA bases. Skipping.\n");
    return Action::SKIP;
  }
  // Determine base connectivity (5' and 3' neighbors).
  Strands_.clear();
  std::vector<int> Visited( setup.Top().Res(Bases_.back().ResNum()).LastAtom(), 0 );
  for (Barray::iterator base = Bases_.begin(); base != Bases_.end(); ++base) {
    Residue const& res = setup.Top().Res( base->ResNum() );
    int c5neighbor = -1;
    int c3neighbor = -1;
    for (int ratom = res.FirstAtom(); ratom != res.LastAtom(); ++ratom) {
      if ( setup.Top()[ratom].Name() == "C5' " ||
           setup.Top()[ratom].Name() == "C5* " )
        c5neighbor = TravelBackbone( setup.Top(), ratom, Visited );
      else if ( setup.Top()[ratom].Name() == "C3' " ||
                setup.Top()[ratom].Name() == "C3* " )
        c3neighbor = TravelBackbone( setup.Top(), ratom, Visited ); 
    }
    std::fill( Visited.begin()+res.FirstAtom(), Visited.begin()+res.LastAtom(), 0 );
    // Convert from atom #s to residue #s
    if (c5neighbor != -1) c5neighbor = setup.Top()[c5neighbor].ResNum();
    if (c3neighbor != -1) c3neighbor = setup.Top()[c3neighbor].ResNum();
    // Find residue #s in Bases_ and set indices.
    for (idx = 0; idx != Bases_.size(); idx++) {
      if (c5neighbor == Bases_[idx].ResNum()) base->SetC5Idx( idx );
      if (c3neighbor == Bases_[idx].ResNum()) base->SetC3Idx( idx );
    }
    if (debug_ > 0)
      mprintf("DEBUG: Base res %i C3res=%i C5res=%i\n", base->ResNum()+1, base->C3resIdx(), base->C5resIdx());
  } // END loop over bases setting 5' and 3' neighbors

  // For each base, find 5' terminii, trace them to the 3' terminii, set up strands
  int strandNum = 0;
  std::vector<bool> base_is_used( Bases_.size(), false );
  for (unsigned int bidx = 0; bidx < Bases_.size(); bidx++) {
    if (!base_is_used[bidx]) {
      if (Bases_[bidx].C5resIdx() < 0) {
        if (debug_ > 0)
          mprintf("DEBUG: 5' start at residue %i\n", Bases_[bidx].ResNum()+1);
        int strandBeg = bidx;
        int strandEnd = follow_base_to_3prime(Bases_, bidx, base_is_used, strandNum);
        if (strandEnd < 0) {
          mprinterr("Internal Error: Could not follow base %i to 3' terminal.\n", Bases_[bidx].ResNum());
          return Action::ERR;
        }
        if (debug_ > 0) {
          mprintf("DEBUG: 3' end at residue %i\n", Bases_[strandEnd].ResNum()+1);
          mprintf("DEBUG: Adding strand %i from base index %i to %i\n", strandNum, strandBeg, strandEnd);
        }
        strandNum++;
        Strands_.push_back( Rpair(strandBeg, strandEnd) );
      }
    }
  } // END loop over bases adding strands

  // DEBUG - Print base connectivity.
  if (debug_ > 0) {
    mprintf("\tBase Connectivity:\n");
    for (Barray::const_iterator base = Bases_.begin(); base != Bases_.end(); ++base) {
      mprintf("\t  Residue %s", setup.Top().TruncResNameNum( base->ResNum() ).c_str());
      if (base->C5resIdx() != -1)
        mprintf(" (5'= %s)", setup.Top().TruncResNameNum( Bases_[base->C5resIdx()].ResNum() ).c_str());
      if (base->C3resIdx() != -1)
        mprintf(" (3'= %s)", setup.Top().TruncResNameNum( Bases_[base->C3resIdx()].ResNum() ).c_str());
      mprintf(" %i\n", base->StrandNum());
    }
    mprintf("\tNucleic acid strands (%zu):\n", Strands_.size());
    for (StrandArray::const_iterator st = Strands_.begin(); st != Strands_.end(); ++st)
      mprintf("\t  %li: Base index %i to %i\n", st-Strands_.begin(), st->first, st->second);
  }
  
  // Determine strand data
  if (sscalc_) {
    for (StrandArray::const_iterator strand = Strands_.begin();
                                     strand != Strands_.end(); ++strand)
    {
      if (debug_ > 0)
        mprintf("DEBUG: Strand %li, %i-%i\n", strand-Strands_.begin(), strand->first, strand->second);
      for (int b1idx = strand->first; b1idx < strand->second; b1idx++)
      {
        int b2idx = b1idx + 1;
        if (debug_ > 0)
          mprintf("DEBUG:\tStrand pair: %i to %i\n", b1idx, b2idx);

        // Get/add data set for strandpair
        Rpair strandpair(b1idx, b2idx);
        Smap::iterator entry = StrandPairs_.lower_bound( strandpair );
        if (entry == StrandPairs_.end() || entry->first != strandpair) {
          // New strand pair 
          Stype SP;
          //  MetaData md = NewStepType(BS, BP1.base1idx_, BP1.base2idx_,
          //                                BP2.base1idx_, BP2.base2idx_, Steps_.size()+1);
          MetaData md(dataname_, StrandPairs_.size()+1);
          md.SetLegend( Bases_[b1idx].BaseName() + "-" + Bases_[b2idx].BaseName() );
          md.SetAspect("dx");
          SP.dx_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
          md.SetAspect("dy");
          SP.dy_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
          md.SetAspect("dz");
          SP.dz_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
          md.SetAspect("rx");
          SP.rx_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
          md.SetAspect("ry");
          SP.ry_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
          md.SetAspect("rz");
          SP.rz_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);

          SP.strandidx_ = (unsigned int)(strand - Strands_.begin());
          //SP.base1idx_ = b1idx_;
          //SP.base2idx_ = b2idx_;

          entry = StrandPairs_.insert( entry, std::pair<Rpair,Stype>(strandpair, SP) );
        }
      } // END loop over bases in strand
    } // END loop over strands
  } // END if sscalc

  // Set up base axes pseudo-topology
  setupNframes_ = setup.Nframes();
  setupTop_ = setup.TopAddress();
  if (axesParm_ != 0) {
    if (axesParm_->Natom() > 0) {
      mprintf("\tBase axes pseudo-topology is already set up.\n");
      // Check that number of bases has not changed
      if ((unsigned int)axesParm_->Nres() != Bases_.size()) {
        mprinterr("Error: Number of bases has changed from %i to %zu.\n"
                  "Error: Base axes pseudo-topology is already set up for %i bases.\n",
                  axesParm_->Nres(), Bases_.size(), axesParm_->Nres());
        return Action::ERR;
      }
    } else {
      // Create residue information for each base axes
      std::vector<Residue> axesResidues;
      axesResidues.reserve( Bases_.size() );
      for (Barray::iterator base = Bases_.begin(); base != Bases_.end(); ++base) {
        Residue const& res = setup.Top().Res( base->ResNum() );
        axesResidues.push_back( res );
      }
      if (setup_axes_pseudoTraj( *axesParm_, *axesOut_, axesFrame_, axesResidues ))
        return Action::ERR;
    }
  }
  return Action::OK;  
}

/** Set up topology for an axes pseudo-trajectory. */
int Action_NAstruct::setup_axes_pseudoTraj(Topology& pseudo,
                                           Trajout_Single& outtraj,
                                           Frame& frame,
                                           std::vector<Residue> const& axesResidues)
const
{
  if (pseudo.Natom() > 0) {
    mprintf("\tAxes pseudo-topology '%s' is already set up.\n", pseudo.c_str());
    return 0;
  }
  // 1 pseudo bond type, Rk = 0.0, Req = 1.0 Ang., will be bond parm index 0
  BondArray bonds;
  pseudo.AddBondParm( BondParmType(0.0, 1.0) );
  int natom = 0;
  for (std::vector<Residue>::const_iterator res = axesResidues.begin();
                                            res != axesResidues.end();
                                          ++res)
  {
    // Order is origin, x, y, z
    pseudo.AddTopAtom(Atom(axisNameO_, "C"), *res);
    pseudo.AddTopAtom(Atom(axisNameX_, "H"), *res);
    pseudo.AddTopAtom(Atom(axisNameY_, "H"), *res);
    pseudo.AddTopAtom(Atom(axisNameZ_, "H"), *res);
    // Bond x y and z to origin
    pseudo.AddBond(natom, natom+1, 0);
    pseudo.AddBond(natom, natom+2, 0);
    pseudo.AddBond(natom, natom+3, 0);
    natom += 4;
  }

  pseudo.CommonSetup();
  if (!pseudo.OriginalFilename().empty()) {
    mprintf("\tWriting axes pseudo-topology to '%s'\n", pseudo.OriginalFilename().full());
    ParmFile pfile;
    if (pfile.WriteTopology(pseudo, pseudo.OriginalFilename(), ParmFile::UNKNOWN_PARM, debug_)) {
      mprinterr("Error: Could not write axes pseudo-topology '%s'\n", pseudo.c_str());
      return 1;
    }
  }

  if (outtraj.SetupTrajWrite( &pseudo, CoordinateInfo(), setupNframes_)) {
    mprinterr("Error: Could not set up axes output trajectory.\n");
    return 1;
  }
  mprintf("      "); //TODO this is a kludge; PrintInfo should be a string.
  outtraj.PrintInfo(0);

  frame = Frame(natom);
  return 0;
} 

// Action_NAstruct::CalculateHbonds()
void Action_NAstruct::CalculateHbonds() {
  for (BPmap::iterator BP = BasePairs_.begin(); BP != BasePairs_.end(); ++BP)
    BP->second.nhb_ = CalcNumHB(Bases_[BP->second.base1idx_], Bases_[BP->second.base2idx_],
                                BP->second.n_wc_hb_);
}

// Action_NAstruct::DoAction()
Action::RetType Action_NAstruct::DoAction(int frameNum, ActionFrame& frm) {
# ifdef NASTRUCTDEBUG
  mprintf("NASTRUCTDEBUG: Frame %i\n", frameNum);
# endif
  if ( findBPmode_ == REFERENCE ) {
    // Base pairs have been determined by a reference, first frame, or
    // guess. Just set up base axes and calculate number of hydrogen bonds.
    if ( SetupBaseAxes(frm.Frm()) ) return Action::ERR;
    CalculateHbonds();
  } else if ( findBPmode_ == ALL ) {
    // Base pairs determined for each individual frame. Hydrogen bonds are 
    // calculated as part of determining base pairing.
    if ( SetupBaseAxes(frm.Frm()) ) return Action::ERR;
    if ( DetermineBasePairing() ) return Action::ERR;
  } else if ( findBPmode_ == SPECIFIED ) {
    // User-specified base pairing.
    if ( SetupBaseAxes(frm.Frm()) ) return Action::ERR;
    if ( SpecifiedBasePairing() ) return Action::ERR;
    findBPmode_ = REFERENCE;
  } else if ( findBPmode_ == FIRST) {
    // Base pairs need to be determined from first frame.
#   ifdef MPI
    // Ensure all processes set up base axes for base pair determinination
    // from the first frame on master.
    int err = 0;
    if (trajComm_.Master()) { // TODO MasterBcast?
      for (int rank = 1; rank < trajComm_.Size(); rank++)
        frm.Frm().SendFrame( rank, trajComm_);
      err = SetupBaseAxes(frm.Frm());
      if (err==0) err = DetermineBasePairing();
    } else {
      Frame refFrame = frm.Frm();
      refFrame.RecvFrame(0, trajComm_);
      err = SetupBaseAxes( refFrame );
      if (err != 0)
        rprinterr("Error: Could not sync nastruct reference first frame.\n");
      else {
        err = DetermineBasePairing();
        // Now re-set up base axes from current frame on non-master
        if (err == 0) err = SetupBaseAxes( frm.Frm() );
        // # hbonds currently based on reference; re-calc for current frame.
        if (err == 0) CalculateHbonds();
      }
    }
    if (trajComm_.CheckError( err )) return Action::ERR;
#   else
    if ( SetupBaseAxes(frm.Frm()) ) return Action::ERR;
    if ( DetermineBasePairing() ) return Action::ERR;
#   endif
    findBPmode_ = REFERENCE;
  }

  // Set up base pair axes pseudo-topology. This is done inside DoAction
  // because we may not know about base pairing until the first frame.
  if (bpAxesParm_ != 0) {
    if (bpAxesParm_->Natom() > 0) {
      //mprintf("\tBase pair axes pseudo-topology is already set up.\n");
      // Check that number of base pairs has not changed
      if ((unsigned int)bpAxesParm_->Nres() != BasePairs_.size()) {
        mprinterr("Error: Number of base pairs has changed from %i to %zu.\n"
                  "Error: Base pair axes pseudo-topology is already set up for %i bases.\n",
                  bpAxesParm_->Nres(), BasePairs_.size(), bpAxesParm_->Nres());
        return Action::ERR;
      }
    } else {
      // Create residue information for each base pair axes
      std::vector<Residue> bpAxesResidues;
      bpAxesResidues.reserve( BasePairs_.size() );
      for (BPmap::const_iterator it = BasePairs_.begin(); it != BasePairs_.end(); ++it)
      {
        //if (it->second.nhb_ < 1 && skipIfNoHB_) continue;
        BPtype const& BP = it->second;
        int b1 = BP.base1idx_;
        //int b2 = BP.base2idx_;
        NA_Base const& base1 = Bases_[b1];
        //NA_Base const& base2 = Bases_[b2];
        Residue const& res = setupTop_->Res( base1.ResNum() );
        bpAxesResidues.push_back( res );
      }
      if (setup_axes_pseudoTraj( *bpAxesParm_, *bpAxesOut_, bpAxesFrame_, bpAxesResidues ))
        return Action::ERR;
    }
  }

  // Determine strand parameters if desired
  if (sscalc_)
    DetermineStrandParameters(frameNum);

  // Output base axes if needed. Do it here because DeterminePairParameters()
  // can flip base axes.
  if (axesOut_ != 0) {
    axesFrame_.ClearAtoms();
    for (std::vector<NA_Base>::const_iterator base = Bases_.begin(); 
                                              base != Bases_.end(); ++base)
      axesToFrame( axesFrame_, base->Axis() );
    axesOut_->WriteSingle(frm.TrajoutNum(), axesFrame_);
  }

  // Determine base parameters
  DeterminePairParameters(frameNum);

  // Determine base pair step parameters
  DetermineStepParameters(frameNum);

  // Set up base pair step axes pseudo-topology. This is done inside DoAction
  // because we may not know about base pairing until the first frame.
  if (stepAxesParm_ != 0) {
    if (stepAxesParm_->Natom() > 0) {
      //mprintf("\tBase pair step axes pseudo-topology is already set up.\n");
      // Check that number of base pair steps has not changed
      if ((unsigned int)stepAxesParm_->Nres() != Steps_.size()) {
        mprinterr("Error: Number of base pair steps has changed from %i to %zu.\n"
                  "Error: Base pair step axes pseudo-topology is already set up for %i bases.\n",
                  stepAxesParm_->Nres(), Steps_.size(), stepAxesParm_->Nres());
        return Action::ERR;
      }
    } else {
      // Create residue information for each base pair step axes
      std::vector<Residue> stepAxesResidues;
      stepAxesResidues.reserve( Steps_.size() );
      for (StepMap::const_iterator it = Steps_.begin(); it != Steps_.end(); ++it)
      {
        StepType const& BS = it->second;
        NA_Base const& base1 = Bases_[BS.b1idx_];
        Residue const& res = setupTop_->Res( base1.ResNum() );
        stepAxesResidues.push_back( res );
      }
      if (setup_axes_pseudoTraj( *stepAxesParm_, *stepAxesOut_, stepAxesFrame_, stepAxesResidues ))
        return Action::ERR;
    }
  }

  // Output base pair axes if needed
  if (bpAxesOut_ != 0) {
    bpAxesFrame_.ClearAtoms();
    for (BPmap::const_iterator it = BasePairs_.begin(); it != BasePairs_.end(); ++it) {
      //if (it->second.nhb_ < 1 && skipIfNoHB_) continue;
      axesToFrame( bpAxesFrame_, it->second.bpaxis_ );
    }
    bpAxesOut_->WriteSingle(frm.TrajoutNum(), bpAxesFrame_);
  }
  // Output base pair step axes if needed
  if (stepAxesOut_ != 0) {
    stepAxesFrame_.ClearAtoms();
    for (StepMap::const_iterator it = Steps_.begin(); it != Steps_.end(); ++it)
      axesToFrame( stepAxesFrame_, it->second.stepaxis_ );
    stepAxesOut_->WriteSingle(frm.TrajoutNum(), stepAxesFrame_);
  }

  nframes_++;
  return Action::OK;
} 

// UpdateTimeSeries()
static inline void UpdateTimeSeries(unsigned int nframes_, DataSet* ds) {
  if (ds != 0) {
    if (ds->Type() == DataSet::FLOAT) {
      const float ZERO = 0.0;
      if (ds->Size() < nframes_) ds->Add(nframes_ - 1, &ZERO);
    } else if (ds->Type() == DataSet::INTEGER) {
      const int ZERO = 0;
      if (ds->Size() < nframes_) ds->Add(nframes_ - 1, &ZERO);
    } else if (ds->Type() == DataSet::VECTOR) {
      double ZERO[6];
      for (int i = 0; i < 6; i++)
        ZERO[i] = 0.0;
      if (ds->Size() < nframes_) ds->Add(nframes_ - 1, ZERO);
    }
  }
}

// Action_NAstruct::NewStepType()
MetaData Action_NAstruct::NewStepType( StepType& BS, int BP1_1, int BP1_2, int BP2_1, int BP2_2,
                                       int idx ) const
{
  // New base pair step
  MetaData md(dataname_, idx);
  md.SetLegend( Bases_[BP1_1].BaseName() + Bases_[BP1_2].BaseName() + "-" +
                Bases_[BP2_1].BaseName() + Bases_[BP2_2].BaseName() );
  md.SetAspect("shift");
  BS.shift_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("slide");
  BS.slide_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("rise");
  BS.rise_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("tilt");
  BS.tilt_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("roll");
  BS.roll_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("twist");
  BS.twist_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("xdisp");
  BS.xdisp_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("ydisp");
  BS.ydisp_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("hrise");
  BS.hrise_  = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("incl");
  BS.incl_   = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("tip");
  BS.tip_    = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("htwist");
  BS.htwist_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  md.SetAspect("zp");
  BS.Zp_     = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
  BS.b1idx_ = BP1_1;
  BS.b2idx_ = BP1_2;
  BS.b3idx_ = BP2_1;
  BS.b4idx_ = BP2_2;
  BS.majGroove_ = 0;
  BS.minGroove_ = 0;
  return md;
}

#ifdef MPI
void Action_NAstruct::NA_Sync( DataSet_1D* dsIn, std::vector<int> const& rank_offsets,
                               std::vector<int> const& rank_frames, int rank, int in ) const
{
  if (dsIn == 0) return;
  //rprintf("Calling sync for '%s'\n", dsIn->legend());
  DataSet_float& DS = static_cast<DataSet_float&>( *dsIn );
  if (trajComm_.Master()) { //TODO put this inside DataSet_float?
    DS.Resize( nframes_ );
    float* d_beg = DS.Ptr() + rank_offsets[ rank ];
    //mprintf("\tResizing nastruct series data to %i, starting frame %i, # frames %i\n",
    //        nframes_, rank_offsets[rank], rank_frames[rank]);
    trajComm_.Recv( d_beg,    rank_frames[ rank ], MPI_FLOAT, rank, 1501 + in );
  } else {
    trajComm_.Send( DS.Ptr(), DS.Size(),           MPI_FLOAT, 0,    1501 + in );
  }
  dsIn->SetNeedsSync( false );
}

// Action_NAstruct::SyncAction()
int Action_NAstruct::SyncAction() {
  // Make sure all time series are updated at this point.
  UpdateSeries();
  // TODO consolidate # frames / offset calc code with Action_Hbond
  // Get total number of frames
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &nframes_, 1, MPI_INT, &rank_frames[0] );
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      nframes_ += rank_frames[rank];
  }
  // Convert rank frames to offsets.
  std::vector<int> rank_offsets( trajComm_.Size(), 0 );
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      rank_offsets[rank] = rank_offsets[rank-1] + rank_frames[rank-1];
  }
  //rprinterr("DEBUG: Number base pairs: %zu\n", BasePairs_.size());
  //rprinterr("DEBUG: Number base pair steps: %zu\n", Steps_.size());
  // Since base pair steps are generated after base pair parameters are
  // calculated, the number of steps and/or the set ordering in the master
  // DataSetList may be different on different ranks. Assume that Bases_
  // and BasePairs_ are the same.
  // Need to know how many steps on each process 
  int num_steps = (int)Steps_.size();
  std::vector<int> nsteps_on_rank;
  if (trajComm_.Master())
    nsteps_on_rank.resize( trajComm_.Size() );
  trajComm_.GatherMaster( &num_steps, 1, MPI_INT, &nsteps_on_rank[0] );
  // bpidx1, bpidx2, b1idx, b2idx, b3idx, b4idx, groove(0=none, 1=maj, 2=min, 3=both)
  const int iSize = 7;
  std::vector<int> iArray;
  if (trajComm_.Master()) {
    for (int rank = 1; rank < trajComm_.Size(); rank++) {
      if (nsteps_on_rank[rank] > 0) {
        //mprintf("DEBUG:\tReceiving %i steps from rank %i.\n", nsteps_on_rank[rank], rank);
        iArray.resize( iSize * nsteps_on_rank[rank] );
        trajComm_.Recv( &iArray[0], iArray.size(), MPI_INT, rank, 1500 );
        int ii = 0;
        for (int in = 0; in != nsteps_on_rank[rank]; in++, ii += iSize) {
          Rpair steppair( iArray[ii], iArray[ii+1] );
          StepMap::iterator entry = Steps_.lower_bound( steppair );
          if (entry == Steps_.end() || entry->first != steppair) {
            // New base pair step
            //mprintf("NEW BASE PAIR STEP: %i %i\n", iArray[0], iArray[1]);
            StepType BS;
            MetaData md = NewStepType( BS, iArray[2], iArray[3],
                                           iArray[4], iArray[5], Steps_.size()+1 );
            // Groove width: iArray[6]=(0=none, 1=maj, 2=min, 3=both)
            if (grooveCalcType_ == HASSAN_CALLADINE && iArray[6] > 0) {
              MetaData md = BS.shift_->Meta();
              if (iArray[6] == 1 || iArray[6] == 3) {
                md.SetAspect("major");
                BS.majGroove_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
              }
              if (iArray[6] == 2 || iArray[6] == 3) {
                md.SetAspect("minor");
                BS.minGroove_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::FLOAT, md);
              }
            }
            entry = Steps_.insert( entry, std::pair<Rpair, StepType>(steppair, BS) );
          }
          //else mprintf("EXISTING BASE PAIR STEP: %i %i\n", entry->first.first, entry->first.second);
          // Synchronize all step data sets from rank.
          NA_Sync( entry->second.shift_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.slide_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.rise_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.tilt_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.roll_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.twist_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.xdisp_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.ydisp_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.hrise_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.incl_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.tip_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.htwist_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.Zp_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.majGroove_, rank_offsets, rank_frames, rank, in );
          NA_Sync( entry->second.minGroove_, rank_offsets, rank_frames, rank, in );
        } // END master loop over steps from rank
      }
    } // END master loop over ranks
    // At this point we have all step sets from all ranks. Mark all step sets
    // smaller than nframes_ as synced and ensure the time series has been
    // updated to reflect overall # frames.
    for (StepMap::iterator step = Steps_.begin(); step != Steps_.end(); ++step) {
      if ((int)step->second.shift_->Size() < nframes_) {
        step->second.shift_->SetNeedsSync(false);
        step->second.slide_->SetNeedsSync(false);
        step->second.rise_->SetNeedsSync(false);
        step->second.tilt_->SetNeedsSync(false);
        step->second.roll_->SetNeedsSync(false);
        step->second.twist_->SetNeedsSync(false);
        step->second.xdisp_->SetNeedsSync(false);
        step->second.ydisp_->SetNeedsSync(false);
        step->second.hrise_->SetNeedsSync(false);
        step->second.incl_->SetNeedsSync(false);
        step->second.tip_->SetNeedsSync(false);
        step->second.htwist_->SetNeedsSync(false);
        step->second.Zp_->SetNeedsSync(false);
        UpdateTimeSeries(nframes_, step->second.shift_);
        UpdateTimeSeries(nframes_, step->second.slide_);
        UpdateTimeSeries(nframes_, step->second.rise_);
        UpdateTimeSeries(nframes_, step->second.tilt_);
        UpdateTimeSeries(nframes_, step->second.roll_);
        UpdateTimeSeries(nframes_, step->second.twist_);
        UpdateTimeSeries(nframes_, step->second.xdisp_);
        UpdateTimeSeries(nframes_, step->second.ydisp_);
        UpdateTimeSeries(nframes_, step->second.hrise_);
        UpdateTimeSeries(nframes_, step->second.incl_);
        UpdateTimeSeries(nframes_, step->second.tip_);
        UpdateTimeSeries(nframes_, step->second.htwist_);
        UpdateTimeSeries(nframes_, step->second.Zp_);
        if (step->second.majGroove_!=0) {
          step->second.majGroove_->SetNeedsSync(false);
          UpdateTimeSeries(nframes_, step->second.majGroove_);
        }
        if (step->second.minGroove_!=0) {
          step->second.minGroove_->SetNeedsSync(false);
          UpdateTimeSeries(nframes_, step->second.minGroove_);
        }
      }
    }
  } else {
    if (Steps_.size() > 0) {
      iArray.reserve( iSize * Steps_.size() );
      for (StepMap::const_iterator step = Steps_.begin(); step != Steps_.end(); ++step) {
        iArray.push_back( step->first.first );
        iArray.push_back( step->first.second );
        iArray.push_back( step->second.b1idx_ );
        iArray.push_back( step->second.b2idx_ );
        iArray.push_back( step->second.b3idx_ );
        iArray.push_back( step->second.b4idx_ );
        // groove(0=none, 1=maj, 2=min, 3=both)
        if ( step->second.majGroove_ != 0 && step->second.minGroove_ != 0 )
          iArray.push_back( 3 );
        else if ( step->second.minGroove_ != 0)
          iArray.push_back( 2 );
        else if ( step->second.majGroove_ != 0)
          iArray.push_back( 1 );
        else
          iArray.push_back( 0 );
      }
      trajComm_.Send( &iArray[0], iArray.size(), MPI_INT, 0, 1500 );
      // Send step data to master.
      int in = 0;
      for (StepMap::const_iterator step = Steps_.begin(); step != Steps_.end(); ++step, ++in) {
        NA_Sync( step->second.shift_, rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.slide_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.rise_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.tilt_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.roll_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.twist_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.xdisp_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.ydisp_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.hrise_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.incl_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.tip_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.htwist_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.Zp_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.majGroove_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
        NA_Sync( step->second.minGroove_ , rank_offsets, rank_frames, trajComm_.Rank(), in );
      }
    }
  }
  return 0;
}
#endif

// Action_NAstruct::UpdateSeries()
void Action_NAstruct::UpdateSeries() {
  if (seriesUpdated_) return;
  if (nframes_ > 0) {
    // Base pair data
    for (BPmap::const_iterator it = BasePairs_.begin(); it != BasePairs_.end(); ++it) {
      BPtype const& BP = it->second;
      UpdateTimeSeries( nframes_, BP.shear_ );
      UpdateTimeSeries( nframes_, BP.stretch_ );
      UpdateTimeSeries( nframes_, BP.stagger_ );
      UpdateTimeSeries( nframes_, BP.buckle_ );
      UpdateTimeSeries( nframes_, BP.prop_ );
      UpdateTimeSeries( nframes_, BP.opening_ );
      UpdateTimeSeries( nframes_, BP.hbonds_ );
      UpdateTimeSeries( nframes_, BP.isBP_ );
      UpdateTimeSeries( nframes_, BP.major_ );
      UpdateTimeSeries( nframes_, BP.minor_ );
#     ifdef NASTRUCTDEBUG
      UpdateTimeSeries( nframes_, BP.axes_nxyz_ );
#     endif
      UpdateTimeSeries( nframes_, Bases_[BP.base1idx_].Pucker() );
      UpdateTimeSeries( nframes_, Bases_[BP.base2idx_].Pucker() );
    }
    // Step and helix data
    for (StepMap::const_iterator it = Steps_.begin(); it != Steps_.end(); ++it) {
      StepType const& BS = it->second;
      UpdateTimeSeries( nframes_, BS.shift_ );
      UpdateTimeSeries( nframes_, BS.slide_ );
      UpdateTimeSeries( nframes_, BS.rise_ );
      UpdateTimeSeries( nframes_, BS.tilt_ );
      UpdateTimeSeries( nframes_, BS.roll_ );
      UpdateTimeSeries( nframes_, BS.twist_ );
      UpdateTimeSeries( nframes_, BS.xdisp_ );
      UpdateTimeSeries( nframes_, BS.ydisp_ );
      UpdateTimeSeries( nframes_, BS.hrise_ );
      UpdateTimeSeries( nframes_, BS.incl_ );
      UpdateTimeSeries( nframes_, BS.tip_ );
      UpdateTimeSeries( nframes_, BS.htwist_ );
      UpdateTimeSeries( nframes_, BS.Zp_ );
      UpdateTimeSeries( nframes_, BS.majGroove_ );
      UpdateTimeSeries( nframes_, BS.minGroove_ );
    }
    // Strand data
    for (Smap::const_iterator it = StrandPairs_.begin(); it != StrandPairs_.end(); ++it) {
      Stype const& SP = it->second;
      UpdateTimeSeries( nframes_, SP.dx_ );
      UpdateTimeSeries( nframes_, SP.dy_ );
      UpdateTimeSeries( nframes_, SP.dz_ );
      UpdateTimeSeries( nframes_, SP.rx_ );
      UpdateTimeSeries( nframes_, SP.ry_ );
      UpdateTimeSeries( nframes_, SP.rz_ );
    }
  }
  // Should only be called once.
  seriesUpdated_ = true;
}

// Output Format Strings
static const char* BP_OUTPUT_FMT = "%8i %8i %8i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %2.0f %2.0f";
static const char* SS_OUTPUT_FMT = "%8i %8i %8i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f";
static const char* GROOVE_FMT = " %10.4f %10.4f";
static const char* HELIX_OUTPUT_FMT = "%8i %4i-%-4i %4i-%-4i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f";
static const char* STEP_OUTPUT_FMT = "%8i %4i-%-4i %4i-%-4i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f";

// Action_NAstruct::Print()
void Action_NAstruct::Print() {
  // Ensure all series have been updated for all frames.
  UpdateSeries();

  // ---------- Base pair parameters -------------
  if (bpout_ != 0) {
    // Check that there is actually data
    if ( BasePairs_.empty() || nframes_ < 1)
      mprintf("Warning: Could not write BP file %s: No BP data.\n", bpout_->Filename().full());
    else {
      mprintf("\tBase pair output file %s; %i frames, %zu base pairs.\n", 
              bpout_->Filename().full(), nframes_, BasePairs_.size());
      //  File header
      if (printheader_) {
        bpout_->Printf("%-8s %8s %8s %10s %10s %10s %10s %10s %10s %2s %2s",
                       "#Frame","Base1","Base2", "Shear","Stretch","Stagger",
                       "Buckle","Propeller","Opening", "BP", "HB");
        if (grooveCalcType_ == PP_OO)
          bpout_->Printf(" %10s %10s", "Major", "Minor");
        bpout_->Printf("\n");
      }
      // Loop over all frames
      for (int frame = 0; frame < nframes_; ++frame) {
        for (BPmap::const_iterator it = BasePairs_.begin();
                                   it != BasePairs_.end(); ++it)
        {
          BPtype const& BP = it->second;
          bpout_->Printf(BP_OUTPUT_FMT, frame+1, 
                         Bases_[BP.base1idx_].ResNum()+1, Bases_[BP.base2idx_].ResNum()+1,
                         BP.shear_->Dval(frame),   BP.stretch_->Dval(frame),
                         BP.stagger_->Dval(frame), BP.buckle_->Dval(frame),
                         BP.prop_->Dval(frame),    BP.opening_->Dval(frame),
                         BP.isBP_->Dval(frame),    BP.hbonds_->Dval(frame));
          if (grooveCalcType_ == PP_OO) 
            bpout_->Printf(GROOVE_FMT, BP.major_->Dval(frame), BP.minor_->Dval(frame));
          bpout_->Printf("\n");
        }
        if (spaceBetweenFrames_) bpout_->Printf("\n");
      }
    }
  }

  // ---------- Single strand parameters ---------
  if (sscalc_ && ssout_ != 0) {
    // Check that there is actually data
    if ( StrandPairs_.empty() || nframes_ < 1)
      mprintf("Warning: Could not write single strand file %s: No SS data.\n",
              ssout_->Filename().full());
    else {
      mprintf("\tSingle strand output file %s; %i frames, %zu base pairs.\n", 
              ssout_->Filename().full(), nframes_, StrandPairs_.size());
      //  File header
      if (printheader_) {
        ssout_->Printf("%-8s %8s %8s %10s %10s %10s %10s %10s %10s\n",
                       "#Frame","Base1","Base2", "DX","DY","DZ", "RX","RY","RZ");
      }
      // Loop over all frames
      for (int frame = 0; frame < nframes_; ++frame) {
        for (Smap::const_iterator it = StrandPairs_.begin();
                                  it != StrandPairs_.end(); ++it)
        {
          Stype const& SP = it->second;
          ssout_->Printf(SS_OUTPUT_FMT, frame+1, 
                         Bases_[it->first.first].ResNum()+1, Bases_[it->first.second].ResNum()+1,
                         SP.dx_->Dval(frame), SP.dy_->Dval(frame), SP.dz_->Dval(frame),
                         SP.rx_->Dval(frame), SP.ry_->Dval(frame), SP.rz_->Dval(frame));
          ssout_->Printf("\n");
        }
        if (spaceBetweenFrames_) ssout_->Printf("\n");
      }
    }
  }

  // ---------- Base pair step parameters --------
  // Check that there is actually data
  if ( Steps_.empty() || nframes_ < 1 ) {
    if (stepout_ != 0 || helixout_ != 0)
      mprintf("Warning: Could not write BPstep / helix files: No data.\n");
  } else {
    if (stepout_ != 0) {
      mprintf("\tBase pair step output file %s\n\tHelix output file %s:\n"
              "\t  %i frames, %zu base pair steps.\n", 
              stepout_->Filename().full(), helixout_->Filename().full(),
              nframes_, Steps_.size());
      // -- Base pair step frames ------
      if (printheader_) {
        stepout_->Printf("%-8s %-9s %-9s %10s %10s %10s %10s %10s %10s %10s","#Frame","BP1","BP2",
                       "Shift","Slide","Rise","Tilt","Roll","Twist","Zp");
        if (grooveCalcType_ == HASSAN_CALLADINE)
          stepout_->Printf(" %10s %10s\n", "Major", "Minor");
        stepout_->Printf("\n");
      }
      for (int frame = 0; frame < nframes_; ++frame) {
        for (StepMap::const_iterator it = Steps_.begin(); it != Steps_.end(); ++it)
        {
          StepType const& BS = it->second;
          // BPstep write
          stepout_->Printf(STEP_OUTPUT_FMT, frame+1, 
                         Bases_[BS.b1idx_].ResNum()+1, Bases_[BS.b2idx_].ResNum()+1,
                         Bases_[BS.b3idx_].ResNum()+1, Bases_[BS.b4idx_].ResNum()+1,
                         BS.shift_->Dval(frame), BS.slide_->Dval(frame),
                         BS.rise_->Dval(frame),  BS.tilt_->Dval(frame),
                         BS.roll_->Dval(frame),  BS.twist_->Dval(frame),
                         BS.Zp_->Dval(frame));
          if (grooveCalcType_ == HASSAN_CALLADINE) {
            if (BS.majGroove_ == 0)
              stepout_->Printf(" %10s", "----");
            else
              stepout_->Printf(" %10.4f", BS.majGroove_->Dval(frame));
            if (BS.minGroove_ == 0)
              stepout_->Printf(" %10s", "----");
            else
              stepout_->Printf(" %10.4f", BS.minGroove_->Dval(frame));
          }
          stepout_->Printf("\n");
        }
        if (spaceBetweenFrames_) stepout_->Printf("\n");
      }
    } // END stepout_ file defined
    // -- Helix frames -----------------
    if (helixout_ != 0) {
      if (printheader_)
        helixout_->Printf("%-8s %-9s %-9s %10s %10s %10s %10s %10s %10s\n","#Frame","BP1","BP2",
                        "X-disp","Y-disp","Rise","Incl.","Tip","Twist");
      for (int frame = 0; frame < nframes_; ++frame) {
        for (StepMap::const_iterator it = Steps_.begin(); it != Steps_.end(); ++it)
        {
          StepType const& BS = it->second;
          // Helix write
          helixout_->Printf(HELIX_OUTPUT_FMT, frame+1,
                            Bases_[BS.b1idx_].ResNum()+1, Bases_[BS.b2idx_].ResNum()+1,
                            Bases_[BS.b3idx_].ResNum()+1, Bases_[BS.b4idx_].ResNum()+1,
                            BS.xdisp_->Dval(frame), BS.ydisp_->Dval(frame),
                            BS.hrise_->Dval(frame), BS.incl_->Dval(frame),
                            BS.tip_->Dval(frame),   BS.htwist_->Dval(frame));
          helixout_->Printf("\n");
        }
        if (spaceBetweenFrames_) helixout_->Printf("\n");
      } // END loop over frames for helix data
    } // END helixout_ file defined
  } // END Steps_ contains data
}
