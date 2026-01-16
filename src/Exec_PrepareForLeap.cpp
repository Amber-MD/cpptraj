#include "Exec_PrepareForLeap.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "LeapInterface.h"
#include "ParmFile.h"
#include "Remote.h"
#include "StringRoutines.h"
#include "Structure/Disulfide.h"
#include "Structure/HisProt.h"
#include "Structure/MetalCenterFinder.h"
#include "Structure/PdbCleaner.h"
#include "Structure/ResStatArray.h"
#include "Structure/SugarBuilder.h"
#include "Structure/Sugar.h"
#include "Trajout_Single.h"
#include <stack> // FindTerByBonds
#include <cctype> // tolower
#include <algorithm> // unique

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Exec_PrepareForLeap::Exec_PrepareForLeap() : Exec(COORDS),
  errorsAreFatal_(true),
  downloadParams_(true),
  bondUnknownResidues_(false),
  debug_(0)
{
  SetHidden(false);
}

/** If file not present, use a default set of residue names. */
void Exec_PrepareForLeap::SetPdbResNames() {
  //Protein
  pdb_res_names_.insert("ACE");
  pdb_res_names_.insert("ALA");
  pdb_res_names_.insert("ARG");
  pdb_res_names_.insert("ASH");
  pdb_res_names_.insert("ASN");
  pdb_res_names_.insert("ASP");
  pdb_res_names_.insert("CYM");
  pdb_res_names_.insert("CYS");
  pdb_res_names_.insert("CYX");
  pdb_res_names_.insert("GLH");
  pdb_res_names_.insert("GLN");
  pdb_res_names_.insert("GLU");
  pdb_res_names_.insert("GLY");
  pdb_res_names_.insert("HIE");
  pdb_res_names_.insert("HIP");
  pdb_res_names_.insert("HIS");
  pdb_res_names_.insert("HYP"); // Recognized by Glycam
  pdb_res_names_.insert("ILE");
  pdb_res_names_.insert("LEU");
  pdb_res_names_.insert("LYN");
  pdb_res_names_.insert("LYS");
  pdb_res_names_.insert("MET");
  pdb_res_names_.insert("NME");
  pdb_res_names_.insert("PHE");
  pdb_res_names_.insert("PRO");
  pdb_res_names_.insert("SER");
  pdb_res_names_.insert("THR");
  pdb_res_names_.insert("TRP");
  pdb_res_names_.insert("TYR");
  pdb_res_names_.insert("VAL");
  // DNA
  pdb_res_names_.insert("DA");
  pdb_res_names_.insert("DC");
  pdb_res_names_.insert("DG");
  pdb_res_names_.insert("DT");
  // RNA
  pdb_res_names_.insert("A");
  pdb_res_names_.insert("C");
  pdb_res_names_.insert("G");
  pdb_res_names_.insert("U");
}

/** Load PDB residue names recognized by Amber FFs from file. */
int Exec_PrepareForLeap::LoadPdbResNames(std::string const& fnameIn)
{
  std::string fname = fnameIn;
  if (fnameIn.empty()) {
    // Check CPPTRAJHOME
    const char* env = getenv("CPPTRAJHOME");
    if (env != 0) {
      fname.assign(env);
      fname += "/dat/PDB_ResidueNames.txt";
      mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", fname.c_str());
    } else {
      // Check AMBERHOME
      env = getenv("AMBERHOME");
      if (env != 0) {
        fname.assign(env);
        fname += "/AmberTools/src/cpptraj/dat/PDB_ResidueNames.txt";
        mprintf("Info: Parameter file path from AMBERHOME variable: '%s'\n", fname.c_str());
      }
    }
  }
  if (fname.empty()) {
    mprintf("Warning: No PDB residue name file specified and/or CPPTRAJHOME not set.\n"
            "Warning: Using standard set of PDB residue names.\n");
    SetPdbResNames();
    return 0;
  }
  mprintf("\tReading PDB residue names from '%s'\n", fname.c_str());

  CpptrajFile infile;
  if (infile.OpenRead(fname)) {
    mprinterr("Error: Could not open PDB residue name file.\n");
    return 1;
  }
  const char* ptr = 0;
  while ( (ptr = infile.NextLine()) != 0 ) {
    ArgList argline( ptr, " " );
    if (argline.Nargs() > 0) {
      if (argline[0][0] != '#') {
        pdb_res_names_.insert( argline[0] );
      }
    }
  }
  infile.CloseFile();

  return 0;
}

/** \return True if residue name is in pdb_to_glycam_ or pdb_res_names_,
  *               or is solvent.
  */
bool Exec_PrepareForLeap::IsRecognizedPdbRes(NameType const& rname,
                                             SugarBuilder const& sugarBuilder)
const
{
  if (sugarBuilder.IsRecognizedPdbSugar(rname))
    return true;
  SetType::const_iterator amberIt = pdb_res_names_.find( rname );
  if (amberIt != pdb_res_names_.end())
    return true;
  if (rname == solventResName_)
    return true;
  return false;
}

/** \return Array of residue numbers with unrecognized PDB res names. */
Exec_PrepareForLeap::Iarray
  Exec_PrepareForLeap::GetUnrecognizedPdbResidues(Topology const& topIn,
                                                  SugarBuilder const& sugarBuilder)
const
{
  Iarray rnums;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    if (!IsRecognizedPdbRes( topIn.Res(ires).Name(), sugarBuilder ))
    {
      mprintf("\t%s is unrecognized.\n", topIn.TruncResNameOnumId(ires).c_str());
      rnums.push_back( ires );
    }
  }
  return rnums;
}

/** Given an array of residue numbers with unrecognized PDB res names,
  * generate an array with true for unrecognized residues that are
  * either isolated or only bound to other unrecognized residues.
  */
Exec_PrepareForLeap::Iarray
  Exec_PrepareForLeap::GetIsolatedUnrecognizedResidues(Topology const& topIn,
                                                       Iarray const& rnums)
const
{
  typedef std::vector<bool> Barray;
  Barray isRecognized(topIn.Nres(), true);
  for (Iarray::const_iterator it = rnums.begin(); it != rnums.end(); ++it)
    isRecognized[ *it ] = false;

  Iarray isolated;
  for (Iarray::const_iterator it = rnums.begin(); it != rnums.end(); ++it)
  {
    bool isIsolated = true;
    Residue const& res = topIn.Res( *it );
    for (int at = res.FirstAtom(); at != res.LastAtom(); ++at)
    {
      for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat)
      {
        if (topIn[*bat].ResNum() != *it)
        {
          // This bonded atom is in another residue. Is that residue recognized?
          if ( isRecognized[ topIn[*bat].ResNum() ] ) {
            // Residue *it is bonded to a recognized residue. Not isolated.
            isIsolated = false;
            break;
          }
        }
      } // END loop over residue atoms bonded atoms
      if (!isIsolated) break;
    } // END loop over residue atoms
    if (isIsolated) {
      mprintf("\t%s is isolated and unrecognized.\n", topIn.TruncResNameOnumId(*it).c_str());
      isolated.push_back( *it );
    }
  } // END loop over unrecognized residues

  return isolated;
}

// -----------------------------------------------------------------------------
/** Determine where molecules end based on connectivity. */
int Exec_PrepareForLeap::FindTerByBonds(Topology& topIn, CharMask const& maskIn)
const
{
  // NOTE: this code is the same algorithm from Topology::NonrecursiveMolSearch
  // TODO use a common molecule search backend
  std::stack<unsigned int> nextAtomToSearch;
  bool unassignedAtomsRemain = true;
  unsigned int currentAtom = 0;
  unsigned int currentMol = 0;
  unsigned int lowestUnassignedAtom = 0;
  Iarray atomMolNum( topIn.Natom(), -1 );
  while (unassignedAtomsRemain) {
    // This atom is in molecule.
    atomMolNum[currentAtom] = currentMol;
    //mprintf("DEBUG:\tAssigned atom %u to mol %u\n", currentAtom, currentMol);
    // All atoms bonded to this one are in molecule.
    for (Atom::bond_iterator batom = topIn[currentAtom].bondbegin();
                             batom != topIn[currentAtom].bondend(); ++batom)
    {
      if (atomMolNum[*batom] == -1) { // -1 is no molecule
        if (topIn[*batom].Nbonds() > 1)
          // Bonded atom has more than 1 bond; needs to be searched.
          nextAtomToSearch.push( *batom );
        else {
          // Bonded atom only bonded to current atom. No more search needed.
          atomMolNum[*batom] = currentMol;
          //mprintf("DEBUG:\t\tAssigned bonded atom %i to mol %u\n", *batom, currentMol);
        }
      }
    }
    if (nextAtomToSearch.empty()) {
      //mprintf("DEBUG:\tNo atoms left in stack. Searching for next unmarked atom.\n");
      // No more atoms to search. Find next unmarked atom.
      currentMol++;
      unsigned int idx = lowestUnassignedAtom;
      for (; idx != atomMolNum.size(); idx++)
        if (atomMolNum[idx] == -1) break;
      if (idx == atomMolNum.size())
        unassignedAtomsRemain = false;
      else {
        currentAtom = idx;
        lowestUnassignedAtom = idx + 1;
      }
    } else {
      currentAtom = nextAtomToSearch.top();
      nextAtomToSearch.pop();
      //mprintf("DEBUG:\tNext atom from stack: %u\n", currentAtom);
    }
  }
  //t_nostack.Stop();
  //t_nostack.WriteTiming(1, "Non-recursive mol search:");
  //return (int)currentMol;
  // For each selected atom, find last atom in corresponding molecule,
  // set corresponding residue as TER.
  int at = 0;
  while (at < topIn.Natom()) {
    // Find the next selected atom
    while (at < topIn.Natom() && !maskIn.AtomInCharMask(at)) at++;
    if (at < topIn.Natom()) {
      int currentMol = atomMolNum[at];
      // Seek to end of molecule
      while (at < topIn.Natom() && currentMol == atomMolNum[at]) at++;
      // The previous atom is the end
      int lastRes = topIn[at-1].ResNum();
      mprintf("\tSetting residue %s as terminal.\n",
              topIn.TruncResNameOnumId(lastRes).c_str());
      topIn.SetRes(lastRes).SetTerminal( true );
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------

/// \return index of oxygen atom bonded to this atom but not in same residue
static inline int getLinkOxygenIdx(Topology const& leaptop, int at, int rnum) {
  int o_idx = -1;
  for (Atom::bond_iterator bat = leaptop[at].bondbegin();
                           bat != leaptop[at].bondend(); ++bat)
  {
    if (leaptop[*bat].Element() == Atom::OXYGEN && leaptop[*bat].ResNum() != rnum) {
      o_idx = *bat;
      break;
    }
  }
  return o_idx;
}

/// \return index of carbon bonded to link oxygen not in same residue
static inline int getLinkCarbonIdx(Topology const& leaptop, int at, int rnum)
{
  int o_idx = getLinkOxygenIdx(leaptop, at, rnum);
  if (o_idx == -1) return o_idx;
  int c_idx = -1;
  for (Atom::bond_iterator bat = leaptop[o_idx].bondbegin();
                           bat != leaptop[o_idx].bondend(); ++bat)
  {
    if (leaptop[*bat].Element() == Atom::CARBON && leaptop[*bat].ResNum() != rnum) {
      c_idx = *bat;
      break;
    }
  }
  return c_idx;
}

/** Run leap to generate topology. Modify the topology if needed. */
int Exec_PrepareForLeap::RunLeap(std::string const& ff_file,
                                 std::string const& leapfilename) const
{
  if (leapfilename.empty()) {
    mprintf("Warning: No leap input file name was specified, not running leap.\n");
    return 0;
  }
  if (ff_file.empty()) {
    mprintf("Warning: No leap input file with force fields was specified, not running leap.\n");
    return 0;
  }
  mprintf("\tExecuting leap.\n");

  std::string topname = leapunitname_ + ".parm7";
  std::string rstname = leapunitname_ + ".rst7";

  Cpptraj::LeapInterface LEAP(debug_);
  LEAP.AddInputFile( ff_file );
  LEAP.AddInputFile( leapfilename );
  LEAP.AddCommand("saveamberparm " + leapunitname_ + " " +
                  topname + " " + rstname);

  if (LEAP.RunLeap()) {
    mprinterr("Error: Leap failed.\n");
    return 1;
  }

  // Load the leap topology;
  Topology leaptop;
  ParmFile parm;
  if (parm.ReadTopology(leaptop, topname, debug_)) return 1;

  bool top_is_modified = false;
  // Go through each residue. Find ones that need to be adjusted.
  // NOTE: If deoxy carbons are ever handled, need to add H1 hydrogen and
  //       add the former -OH charge to the carbon.
  for (int rnum = 0; rnum != leaptop.Nres(); rnum++)
  {
    Residue const& res = leaptop.Res(rnum);
    if (res.Name() == "SO3") {
      int o_idx = -1;
      // Need to adjust the charge on the bonded oxygen by +0.031
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        if (leaptop[at].Element() == Atom::SULFUR) {
          o_idx = getLinkOxygenIdx( leaptop, at, rnum );
          if (o_idx != -1) break;
        }
      }
      if (o_idx == -1) {
        mprinterr("Error: Could not find oxygen link atom for '%s'\n",
                  leaptop.TruncResNameOnumId(rnum).c_str());
        return 1;
      }
      double newcharge = leaptop[o_idx].Charge() + 0.031;
      mprintf("\tFxn group '%s'; changing charge on %s from %f to %f\n", *(res.Name()),
              leaptop.AtomMaskName(o_idx).c_str(), leaptop[o_idx].Charge(), newcharge);
      leaptop.SetAtom(o_idx).SetCharge( newcharge );
      top_is_modified = true;
    } else if (res.Name() == "MEX") {
      int c_idx = -1;
      // Need to adjust the charge on the carbon bonded to link oxygen by -0.039
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        if (leaptop[at].Element() == Atom::CARBON) {
          c_idx = getLinkCarbonIdx( leaptop, at, rnum );
          if (c_idx != -1) break;
        }
      }
      if (c_idx == -1) {
        mprinterr("Error: Could not find carbon bonded to oxygen link atom for '%s'\n",
                  leaptop.TruncResNameOnumId(rnum).c_str());
        return 1;
      }
      double newcharge = leaptop[c_idx].Charge() - 0.039;
      mprintf("\tFxn group '%s'; changing charge on %s from %f to %f\n", *(res.Name()),
              leaptop.AtomMaskName(c_idx).c_str(), leaptop[c_idx].Charge(), newcharge);
      leaptop.SetAtom(c_idx).SetCharge( newcharge );
      top_is_modified = true;
    } else if (res.Name() == "ACX") {
      int c_idx = -1;
      // Need to adjust the charge on the carbon bonded to link oxygen by +0.008
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        if (leaptop[at].Element() == Atom::CARBON) {
          // This needs to be the acetyl carbon, ensure it is bonded to an oxygen
          for (Atom::bond_iterator bat = leaptop[at].bondbegin();
                                   bat != leaptop[at].bondend(); ++bat)
          {
            if (leaptop[*bat].Element() == Atom::OXYGEN) {
              c_idx = getLinkCarbonIdx( leaptop, at, rnum );
              if (c_idx != -1) break;
            }
          }
          if (c_idx != -1) break;
        }
      }
      if (c_idx == -1) {
        mprinterr("Error: Could not find carbon bonded to oxygen link atom for '%s'\n",
                  leaptop.TruncResNameOnumId(rnum).c_str());
        return 1;
      }
      double newcharge = leaptop[c_idx].Charge() + 0.008;
      mprintf("\tFxn group '%s'; changing charge on %s from %f to %f\n", *(res.Name()),
              leaptop.AtomMaskName(c_idx).c_str(), leaptop[c_idx].Charge(), newcharge);
      leaptop.SetAtom(c_idx).SetCharge( newcharge );
      top_is_modified = true;
    }
  }

  // DEBUG: Print out total charge on each residue
  double total_q = 0;
  for (Topology::res_iterator res = leaptop.ResStart(); res != leaptop.ResEnd(); ++res)
  {
    double tcharge = 0;
    for (int at = res->FirstAtom(); at != res->LastAtom(); ++at) {
      total_q += leaptop[at].Charge();
      tcharge += leaptop[at].Charge();
    }
    if (debug_ > 0) {
      mprintf("DEBUG:\tResidue %10s charge= %12.5f\n",
              leaptop.TruncResNameOnumId(res-leaptop.ResStart()).c_str(), tcharge);
    }
  }
  mprintf("\tTotal charge: %16.8f\n", total_q);

  // If topology was modified, write it back out
  if (top_is_modified) {
    mprintf("\tWriting modified topology back to '%s'\n", topname.c_str());
    parm.WriteTopology(leaptop, topname, parm.CurrentFormat(), debug_);
  }

  return 0;
}

/** Try to download missing parameters. */
int Exec_PrepareForLeap::DownloadParameters(ResStatArray& resStat, RmapType const& resNames,
                                            Topology const& topIn, CpptrajFile* leapInput,
                                            std::vector<BondType>& LeapBonds)
const
{
  Cpptraj::Remote remote( parameterURL_ );
  remote.SetDebug(debug_);
  for (RmapType::const_iterator it = resNames.begin(); it != resNames.end(); ++it)
  {
    std::string rname = it->first.Truncated();
    mprintf("\t\tSearching for parameters for residue '%s'", rname.c_str());
    for (Iarray::const_iterator rn = it->second.begin(); rn != it->second.end(); ++rn)
      mprintf(" %i", *rn + 1);
    mprintf("\n");
    // Assume parameters are in a subdirectory starting with lowercase version
    // of the first letter of the residue.
    char lcase = tolower( rname[0] );
    std::string rfbase = std::string(1, lcase) + "/" + rname;
    if (debug_ > 0)
      mprintf("DEBUG: Base name: %s\n", rfbase.c_str());
    int err = 0;
    err += remote.DownloadFile( rfbase + ".mol2" );
    if (err == 0) {
      err += remote.DownloadFile( rfbase + ".frcmod" );
      if (err != 0)
        mprintf("Warning: Could not download %s.frcmod\n", rname.c_str());
    } else {
      mprintf("Warning: Could not download %s.mol2\n", rname.c_str());
    }
    if (err != 0)
      mprintf("Warning: Could not download parameter files for '%s'\n", rname.c_str());
    else {
      // Sanity check - make sure the files are there.
      if (!File::Exists(rname + ".mol2") || !File::Exists(rname + ".frcmod")) {
        mprinterr("Error: Problem downloading parameter files for '%s'\n", rname.c_str());
        return 1;
      }
      // Add leap input
      if (leapInput != 0) {
        leapInput->Printf("%s = loadmol2 %s.mol2\n", rname.c_str(), rname.c_str());
        leapInput->Printf("parm%s = loadamberparams %s.frcmod\n", rname.c_str(), rname.c_str());
      }
      // Downloaded mol2 files do not have connectivity, so we need to add bonds.
      for (Iarray::const_iterator rn = it->second.begin(); rn != it->second.end(); ++rn)
      {
        // Assume we are all good with parameters for these residues now.
        resStat[*rn] = ResStatArray::VALIDATED;
        Residue const& currentRes = topIn.Res( *rn );
        // Check each atom for bonds to another residue
        for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); ++at) {
          for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat) {
            if ( topIn[*bat].ResNum() != *rn ) {
              if (bondUnknownResidues_) {
                mprintf("\t\t\tRes %s is bonded to res %s\n",
                        topIn.TruncResAtomName(at).c_str(),
                        topIn.TruncResAtomName(*bat).c_str());
                LeapBonds.push_back( BondType(at, *bat, -1) );
              } else {
                mprintf("Info:\t\t\tRes %s might be bonded to res %s\n",
                        topIn.TruncResAtomName(at).c_str(),
                        topIn.TruncResAtomName(*bat).c_str());

              }
            }
          } // END loop over bonded atoms
        } // END loop over atoms in residue
      } // END loop over residue numbers
    }
  } // END loop over residue names to get parameters for
  return 0;
}

// Exec_PrepareForLeap::Help()
void Exec_PrepareForLeap::Help() const
{
  mprintf("\tcrdset <coords set> [frame <#>] name <out coords set>\n"
          "\t[pdbout <pdbfile> [terbymol]] [problemout <file>]\n"
          "\t[leapunitname <unit>] [out <leap input file> [runleap <ff file>]]\n"
          "\t[skiperrors] [{dlparams|nodlparams}] [{bondunknown|nobondunknown}]\n"
          "\t[nowat [watermask <watermask>] [noh]\n"
          "\t[keepaltloc {<alt loc ID>|highestocc}]\n"
          "\t[stripmask <stripmask>] [solventresname <solventresname>]\n"
          "\t[molmask <molmask> ...] [determinemolmask <mask>]\n"
          "%s"
          "%s"
          "\t[{nosugars |\n"
          "\t  sugarmask <sugarmask> [noc1search] [nosplitres]\n"
          "\t  [rescut <residue cutoff>] [bondoffset <offset>]\n"
          "\t  [resmapfile <file>]\n"
          "\t  [hasglycam] [determinesugarsby {geom|name}]\n"
          "\t }]\n"
          "  Prepare the structure in the given coords set for easier processing\n"
          "  with the LEaP program from AmberTools. Any existing/potential\n"
          "  disulfide bonds will be identified and the residue names changed\n"
          "  to <name> (CYX by default), and if specified any sugars\n"
          "  recognized in the <sugarmask> region will be identified and have\n"
          "  their names changed to Glycam names. Disulfides and sugars will\n"
          "  have any inter-residue bonds removed, and the appropriate LEaP\n"
          "  input to add the bonds back once the structure has been loaded\n"
          "  into LEaP will be written to <leap input file>.\n"
          "  The command will attempt to download parameters for unknown\n"
          "  residues unless 'nodlparams' is specified.\n",
          HisProt::keywords_,
          Disulfide::keywords_
         );
}

// Exec_PrepareForLeap::Execute()
Exec::RetType Exec_PrepareForLeap::Execute(CpptrajState& State, ArgList& argIn)
{
  mprintf("\tPREPAREFORLEAP:\n");
  mprintf("# Citation: Roe, D.R.; Bergonzo, C.; \"PrepareForLeap: An Automated Tool for\n"
          "#           Fast PDB-to-Parameter Generation.\"\n"
          "#           J. Comp. Chem. (2022), V. 43, I. 13, pp 930-935.\n" );
  debug_ = State.Debug();
  errorsAreFatal_ = !argIn.hasKey("skiperrors");
  if (argIn.hasKey("dlparams"))
    downloadParams_ = true;
  else if (argIn.hasKey("nodlparams"))
    downloadParams_ = false;
  else
    downloadParams_ = true;
  if (argIn.hasKey("bondunknown"))
    bondUnknownResidues_ = true;
  else if (argIn.hasKey("nobondunknown"))
    bondUnknownResidues_ = false;
  else
    bondUnknownResidues_ = false;
  // Get input coords
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify input COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (ds == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)ds) );
  // Get frame from input coords
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);

  // Copy input topology, may be modified.
  Topology topIn = coords.Top();

  // Allocate output COORDS data set
  std::string outname = argIn.GetStringKey("name");
  if (outname.empty()) {
    mprinterr("Error: Must specify output COORDS set with 'name'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords_CRD* outCoords = (DataSet_Coords_CRD*)State.DSL().AddSet( DataSet::COORDS, outname );
  if (outCoords == 0) {
    mprinterr("Error: Could not allocate output COORDS set.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tPrepared system will be saved to COORDS set '%s'\n", outCoords->legend());

  std::string leapffname = argIn.GetStringKey("runleap");
  if (!leapffname.empty()) {
#   ifdef _MSC_VER
    mprinterr("Error: Cannot use LEaP interface on windows.\n");
    return CpptrajState::ERR;
#   else
    mprintf("\tWill attempt to run leap with force fields specified in file '%s'\n",
            leapffname.c_str());
#   endif
  }

  std::string pdbout = argIn.GetStringKey("pdbout");
  if (!pdbout.empty())
    mprintf("\tPDB will be written to %s\n", pdbout.c_str());
  else {
    if (!leapffname.empty()) {
      mprinterr("Error: Must specify PDB file name with 'pdbout' if 'runleap' specified.\n");
      return CpptrajState::ERR;
    }
  }
  std::string pdb_ter_arg;
  if (!argIn.hasKey("terbymol")) {
    mprintf("\tUsing original TER cards where possible.\n");
    pdb_ter_arg.assign("pdbter");
  } else
    mprintf("\tGenerating TER cards based on molecular connectivity.\n");

  std::string leapfilename = argIn.GetStringKey("out");
  if (!leapfilename.empty())
    mprintf("\tWriting leap input to '%s'\n", leapfilename.c_str());
  else {
    if (!leapffname.empty()) {
      mprinterr("Error: Must specify leap input file name with 'out' if 'runleap' specified.\n");
      return CpptrajState::ERR;
    }
  }
  leapunitname_ = argIn.GetStringKey("leapunitname", "m");
  mprintf("\tUsing leap unit name: %s\n", leapunitname_.c_str());
  if (validDouble(leapunitname_))
    mprintf("Warning: LEaP unit name '%s' is a valid number; this may confuse some LEaP commands.\n",
             leapunitname_.c_str());
  solventResName_ = argIn.GetStringKey("solventresname", "HOH");
  mprintf("\tSolvent residue name: %s\n", solventResName_.c_str());
  // TODO functional group stuff should be in a file

  bool prepare_sugars = !argIn.hasKey("nosugars");
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");

  // Load PDB residue names recognized by amber
  if (LoadPdbResNames( argIn.GetStringKey("resnamefile" ) )) {
    mprinterr("Error: PDB residue name file load failed.\n");
    return CpptrajState::ERR;
  }
  mprintf("\t%zu PDB residue names recognized by Amber FFs.\n", pdb_res_names_.size());
  // DEBUG
  if (debug_ > 0) {
    mprintf("\tPDB residue names recognized by Amber FFs:\n");
    for (SetType::const_iterator it = pdb_res_names_.begin(); it != pdb_res_names_.end(); ++it)
      mprintf("\t  %s\n", *(*it));
  }
  if (downloadParams_)
    mprintf("\tWill attempt to download parameters for unknown residues.\n");
  else
    mprintf("\tWill not attempt to download paramters for unknown residues.\n");
  if (bondUnknownResidues_)
    mprintf("\tWill attempt to bond unknown residues.\n");
  else
    mprintf("\tWill not attempt to bond unknown residues.\n");
  std::string problemoutname = argIn.GetStringKey("problemout");
  mprintf("\tWriting detected problems to %s\n", problemoutname.c_str());

  // Load PDB to glycam residue name map
  SugarBuilder sugarBuilder(debug_);
  if (prepare_sugars) {
    // Init options
    if (sugarBuilder.InitOptions( argIn.hasKey("hasglycam"),
                                  argIn.getKeyDouble("rescut", 8.0),
                                  argIn.getKeyDouble("bondoffset", 0.2),
                                  argIn.GetStringKey("sugarmask"),
                                  argIn.GetStringKey("determinesugarsby", "geometry"),
                                  argIn.GetStringKey("resmapfile") ))
    {
      mprinterr("Error: Sugar options init failed.\n");
      return CpptrajState::ERR;
    }
  }

  // Do histidine detection before H atoms are removed
  if (!argIn.hasKey("nohisdetect")) {
    Cpptraj::Structure::HisProt hisProt;
    if (hisProt.InitHisProt( argIn, debug_ )) {
      mprinterr("Error: Could not initialize histidine detection.\n");
      return CpptrajState::ERR;
    }
    hisProt.HisProtInfo();
    // Add epsilon, delta, and double-protonated names as recognized.
    pdb_res_names_.insert( hisProt.EpsilonProtHisName() );
    pdb_res_names_.insert( hisProt.DeltaProtHisName() );
    pdb_res_names_.insert( hisProt.DoubleProtHisName() );
    if (hisProt.DetermineHisProt( topIn )) {
      mprinterr("Error: HIS protonation detection failed.\n");
      return CpptrajState::ERR;
    }
  }

  Iarray pdbResToRemove;
  std::string removeArg = argIn.GetStringKey("remove");
  if (!removeArg.empty()) {
    if (removeArg == "unrecognized") {
      mprintf("\tRemoving unrecognized PDB residues.\n");
      pdbResToRemove = GetUnrecognizedPdbResidues( topIn, sugarBuilder );
    } else if (removeArg == "isolated") {
      mprintf("\tRemoving unrecognized and isolated PDB residues.\n");
      Iarray unrecognizedPdbRes = GetUnrecognizedPdbResidues( topIn, sugarBuilder );
      pdbResToRemove = GetIsolatedUnrecognizedResidues( topIn, unrecognizedPdbRes );
    } else {
      mprinterr("Error: Unrecognized keyword for 'remove': %s\n", removeArg.c_str());
      return CpptrajState::ERR;
    }
  }

  // Deal with any coordinate modifications
  Cpptraj::Structure::PdbCleaner pdbCleaner;
  pdbCleaner.SetDebug( debug_ );
  if (pdbCleaner.InitPdbCleaner( argIn, solventResName_, pdbResToRemove )) {
    mprinterr("Error: Could not init PDB cleaner.\n");
    return CpptrajState::ERR;
  }
  if (pdbCleaner.SetupPdbCleaner( topIn )) {
    mprinterr("Error: Could not set up PDB cleaner.\n");
    return CpptrajState::ERR;
  }
  pdbCleaner.PdbCleanerInfo();
  if (pdbCleaner.ModifyCoords(topIn, frameIn)) {
    mprinterr("Error: Could not clean PDB.\n");
    return CpptrajState::ERR;
  }

  // If preparing sugars, need to set up an atom map and potentially
  // search for terminal sugars/missing bonds. Do this here after all atom
  // modifications have been done.
  if (prepare_sugars) {
    bool splitres = !argIn.hasKey("nosplitres");
    if (splitres)
      mprintf("\tWill split off recognized sugar functional groups into separate residues.\n");
    else
      mprintf("\tNot splitting recognized sugar functional groups into separate residues.\n");
    bool c1bondsearch = !argIn.hasKey("noc1search");
    if (c1bondsearch)
      mprintf("\tWill search for missing bonds to sugar anomeric atoms.\n");
    else
      mprintf("\tNot searching for missing bonds to sugar anomeric atoms.\n");
    // May need to modify sugar structure/topology, either by splitting
    // C1 hydroxyls of terminal sugars into ROH residues, and/or by
    // adding missing bonds to C1 atoms.
    // This is done before any identification takes place since we want
    // to identify based on the most up-to-date topology.
    // NOTE: The c1 bonds will remain, so put them in a separate array
    std::vector<BondType> C1Bonds;
    if (sugarBuilder.FixSugarsStructure(topIn, frameIn,
                                        c1bondsearch, splitres, solventResName_,
                                        C1Bonds))
    {
      mprinterr("Error: Sugar structure modification failed.\n");
      return CpptrajState::ERR;
    }
    // NOTE: If IdSugarRing() is to be used after this point, the map
    //       will need to be recreated.
    // Since FixSugarsStructure() can re-order atoms, need
    // to recreate the map.
    //myMap_.ClearMap();
    //if (myMap_.Setup(topIn, frameIn)) {
    //  mprinterr("Error: Atom map second setup failed\n");
    //  return CpptrajState::ERR;
    //}
    //myMap_.DetermineAtomIDs();
  }

  // ----- Below here, no more removing/reordering atoms. ------------ 

  // Each residue starts out unknown.
  ResStatArray resStat( topIn.Nres() );

  // Get masks for molecules now since bond info in topology may be modified later.
  std::vector<AtomMask> molMasks;
  std::string mstr = argIn.GetStringKey("molmask");
  while (!mstr.empty()) {
    mprintf("\tAll atoms selected by '%s' will be in same molecule.\n", mstr.c_str());
    molMasks.push_back( AtomMask() );
    if (molMasks.back().SetMaskString( mstr )) {
      mprinterr("Error: Invalid mask.\n");
      return CpptrajState::ERR;
    }
    if (topIn.SetupIntegerMask( molMasks.back() )) return CpptrajState::ERR;
    molMasks.back().MaskInfo();
    if (molMasks.back().None()) {
      mprinterr("Error: Nothing selected by mask.\n");
      return CpptrajState::ERR;
    }
    mstr = argIn.GetStringKey("molmask");
  }
  CharMask determineMolMask;
  mstr = argIn.GetStringKey("determinemolmask");
  if (!mstr.empty()) {
    mprintf("\tAtoms in mask '%s' will determine molecules by bonds.\n", mstr.c_str());
    if (determineMolMask.SetMaskString(mstr)) {
      mprinterr("Error: Invalid mask.\n");
      return CpptrajState::ERR;
    }
    if (topIn.SetupCharMask( determineMolMask )) return CpptrajState::ERR;
    determineMolMask.MaskInfo();
    if (determineMolMask.None()) {
      mprinterr("Error: Nothing selected by mask.\n");
      return CpptrajState::ERR;
    }
  }

  //CpptrajFile* outfile = State.DFL().AddCpptrajFile(leapfilename,
  //                                                  "LEaP Input", DataFileList::TEXT, true);
  // NOTE: This needs to contain ONLY leap input, so dont put it on the master file list
  CpptrajFile LEAPOUT;
  if (LEAPOUT.OpenWrite(leapfilename)) return CpptrajState::ERR;
  CpptrajFile* outfile = &LEAPOUT;
  if (outfile == 0) return CpptrajState::ERR;
  mprintf("\tLEaP input containing 'loadpdb' and bond commands for disulfides,\n"
          "\t  sugars, etc will be written to '%s'\n", outfile->Filename().full());

  // Array that will hold bonds that need to be made in LEaP
  std::vector<BondType> LeapBonds;

  // Disulfide search
  if (!argIn.hasKey("nodisulfides")) {
    Cpptraj::Structure::Disulfide disulfide;
    if (disulfide.InitDisulfide( argIn, Disulfide::NO_ADD_BONDS, debug_ )) {
      mprinterr("Error: Could not init disulfide search.\n");
      return CpptrajState::ERR;
    }
    if (disulfide.SearchForDisulfides( resStat, topIn, frameIn, LeapBonds ))
    {
      mprinterr("Error: Disulfide search failed.\n");
      return CpptrajState::ERR;
    }
  } else {
    mprintf("\tNot searching for disulfides.\n");
  }

  // Metal center search
  if (!argIn.hasKey("nometals")) {
    Cpptraj::Structure::MetalCenterFinder MC;
    if (MC.InitMetalCenters( argIn, debug_ )) {
      mprinterr("Error: Could not init metal center search.\n");
      return CpptrajState::ERR;
    }
    if (MC.FindMetalCenters( topIn, frameIn )) {
      mprinterr("Error: Metal center search failed.\n");
      return CpptrajState::ERR;
    }
    MC.PrintMetalCenters(topIn);
  } else {
    mprintf("\tNot searching for metal centers.\n");
  }

  // Prepare sugars
  if (prepare_sugars) {
    std::vector<BondType> SugarBonds; 
    if (sugarBuilder.PrepareSugars(errorsAreFatal_, resStat, topIn, frameIn, SugarBonds))
    {
      mprinterr("Error: Sugar preparation failed.\n");
      return CpptrajState::ERR;
    }
    // Remove bonds to sugar
    if (!SugarBonds.empty()) {
      for (std::vector<BondType>::const_iterator it = SugarBonds.begin();
                                                 it != SugarBonds.end(); ++it)
      {
        LeapBonds.push_back( *it );
        topIn.RemoveBond(it->A1(), it->A2());
      }
      // Bonds to sugars have been removed, so regenerate molecule info
      topIn.DetermineMolecules();
    }
    // Set each sugar residue as a terminal residue
    sugarBuilder.SetEachSugarAsTerminal(topIn);
  } else {
    mprintf("\tNot preparing sugars.\n");
  }

  // Determine unknown residues we may want to find parameters for
  NameType solvName(solventResName_);
  RmapType residuesToFindParamsFor;
  for (ResStatArray::iterator it = resStat.begin(); it != resStat.end(); ++it)
  {
    NameType const& residueName = topIn.Res(it-resStat.begin()).Name();
    // Skip water if not removing water
    if (!pdbCleaner.RemoveWater() && residueName == solvName) continue;
    // If status is unknown, see if this is a common residue name
    if ( *it == ResStatArray::UNKNOWN ) {
      SetType::const_iterator pname = pdb_res_names_.find( residueName );
      if (pname == pdb_res_names_.end()) {
        mprintf("\t%s is an unrecognized name and may not have parameters.\n",
                topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        RmapType::iterator ret = residuesToFindParamsFor.lower_bound( residueName );
        if (ret == residuesToFindParamsFor.end() || ret->first != residueName)
        {
          // New residue to get params for
          ret = residuesToFindParamsFor.insert(ret, RpairType(residueName, Iarray()));
        }
        ret->second.push_back( it - resStat.begin() );
      } else
        *it = ResStatArray::VALIDATED;
    }
  }
  if (downloadParams_ && !residuesToFindParamsFor.empty()) {
    mprintf("\tResidues to find parameters for:");
    for (RmapType::const_iterator it = residuesToFindParamsFor.begin();
                                  it != residuesToFindParamsFor.end(); ++it)
    {
      mprintf(" %s", it->first.Truncated().c_str());
      for (Iarray::const_iterator rn = it->second.begin(); rn != it->second.end(); ++rn)
        mprintf(" %i", *rn + 1);
    }
    mprintf("\n");
    // Set default parameter URL if not yet set.
    if (parameterURL_.empty())
      parameterURL_.assign("https://github.com/phenix-project/geostd/raw/refs/heads/master");
    if (DownloadParameters(resStat, residuesToFindParamsFor, topIn, outfile, LeapBonds)) {
      mprinterr("Error: Download parameters failed.\n");
      return CpptrajState::ERR;
    }
  }

  // Add the loadpdb command if we are writing a PDB file.
  // TODO add 'addPdbResMap { { 1 "NH2" "NHE" } }' to recognize NHE?
  if (!pdbout.empty())
    outfile->Printf("%s = loadpdb %s\n", leapunitname_.c_str(), pdbout.c_str());

  // Remove any duplicate bonds
  std::sort( LeapBonds.begin(), LeapBonds.end() );
  std::vector<BondType>::const_iterator it = std::unique( LeapBonds.begin(), LeapBonds.end() );
  LeapBonds.resize( it - LeapBonds.begin() );

  // Create LEaP input for bonds that need to be made in LEaP
  for (std::vector<BondType>::const_iterator bnd = LeapBonds.begin();
                                             bnd != LeapBonds.end(); ++bnd)
    outfile->Printf("%s\n", Cpptraj::LeapInterface::LeapBond(bnd->A1(), bnd->A2(), leapunitname_, topIn).c_str());

  // Count any solvent molecules
  if (!pdbCleaner.RemoveWater()) {
    unsigned int nsolvent = 0;
    for (Topology::res_iterator res = topIn.ResStart(); res != topIn.ResEnd(); ++res) {
      if ( res->Name() == solvName) {
        nsolvent++;
        resStat[res-topIn.ResStart()] = ResStatArray::VALIDATED;
        // Set as terminal; TODO is this needed? Leap seems ok with not having TER for HOH
        topIn.SetRes(res-topIn.ResStart()).SetTerminal(true);
      }
    }
    if (nsolvent > 0) mprintf("\t%u solvent residues.\n", nsolvent);
  }

  // Residue validation.
  bool structure_has_problems = false;
  for (ResStatArray::const_iterator it = resStat.begin(); it != resStat.end(); ++it)
  {
    if (*it != ResStatArray::VALIDATED) {
      structure_has_problems = true;
      break;
    }
  }
  CpptrajFile* problemout = 0;
  if (problemoutname.empty() || !structure_has_problems) {
    problemout = State.DFL().AddCpptrajFile("", "Structure Problems", DataFileList::TEXT, true);
  } else if (structure_has_problems) {
    // Only open the file if there are problems
    problemout = State.DFL().AddCpptrajFile(problemoutname, "Structure Problems",
                                            DataFileList::TEXT, true);
  }
  if (problemout == 0) {
    mprinterr("Internal Error: Could not allocate problemout file.\n");
    return CpptrajState::ERR;
  }
  // Print warnings related to functional groups
  sugarBuilder.PrintFxnGroupWarnings(topIn);
  //mprintf("\tResidues with potential problems:\n");
  int fatal_errors = 0;
  int potential_errors = 0;
  static const char* msg1 = "Potential problem : ";
  static const char* msg2 = "Fatal problem     : ";
  for (ResStatArray::const_iterator it = resStat.begin(); it != resStat.end(); ++it)
  {
    //if ( *it == VALIDATED )
    //  mprintf("\t\t%s VALIDATED\n", topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    //else
    //  mprintf("\t\t%s UNKNOWN\n", topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    // ----- Warnings --------
    if ( *it == ResStatArray::UNKNOWN ) {
      //SetType::const_iterator pname = pdb_res_names_.find( topIn.Res(it-resStat.begin()).Name() );
      //if (pname == pdb_res_names_.end())
        problemout->Printf("\t%s%s is an unrecognized name and may not have parameters.\n",
                msg1, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        potential_errors++;
      //else
      //  *it = ResStatArray::VALIDATED;
    } else if ( *it == ResStatArray::SUGAR_NAME_MISMATCH ) {
        problemout->Printf("\t%s%s sugar anomer type and/or configuration is not consistent with name.\n",
                msg1, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        potential_errors++;
    // ----- Fatal Errors ----
    } else if ( *it == ResStatArray::SUGAR_UNRECOGNIZED_LINK_RES ) {
        problemout->Printf("\t%s%s is linked to a sugar but has no sugar-linkage form.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    } else if ( *it == ResStatArray::SUGAR_UNRECOGNIZED_LINKAGE ) {
        problemout->Printf("\t%s%s is a sugar with an unrecognized linkage.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    } else if ( *it == ResStatArray::SUGAR_NO_LINKAGE ) {
        problemout->Printf("\t%s%s is an incomplete sugar with no linkages.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    } else if ( *it == ResStatArray::SUGAR_NO_CHAIN_FOR_LINK ) {
        problemout->Printf("\t%s%s could not identify chain atoms for determining linkages.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
/*    } else if ( *it == SUGAR_MISSING_C1X ) { // TODO should this be a warning
        mprintf("\t%s%s Sugar is missing anomeric carbon substituent.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;*/
    } else if ( *it == ResStatArray::SUGAR_SETUP_FAILED ) {
        problemout->Printf("\t%s%s Sugar setup failed and could not be identified.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    }
  }
  if ((fatal_errors + potential_errors) > 0 && !problemout->IsStream())
    mprintf("Warning: Check '%s' for structure problems.\n", problemout->Filename().full());

  // Try to set terminal residues
  if (!molMasks.empty() || determineMolMask.MaskStringSet()) {
    // Reset terminal status
    for (int rnum = 0; rnum != topIn.Nres(); rnum++)
      topIn.SetRes(rnum).SetTerminal(false);
    // The final residue of each molMask is terminal
    for (std::vector<AtomMask>::const_iterator mask = molMasks.begin();
                                               mask != molMasks.end(); ++mask)
    {
      //std::vector<int> Rnums = coords.Top().ResnumsSelectedBy( *mask );
      int lastAtom = mask->back();
      int lastRes = topIn[lastAtom].ResNum();
      mprintf("\tSetting residue %s as terminal.\n",
        topIn.TruncResNameOnumId(lastRes).c_str());
      topIn.SetRes(lastRes).SetTerminal( true );
    }
    // Set ter based on connectivity
    if (determineMolMask.MaskStringSet()) {
      if (FindTerByBonds(topIn, determineMolMask)) {
        mprinterr("Error: Could not set TER by connectivity.\n");
        return CpptrajState::ERR;
      }
    }
  }

  // Setup output COORDS
  outCoords->CoordsSetup( topIn, coords.CoordsInfo() );
  outCoords->AddFrame( frameIn );

  if (!pdbout.empty()) {
    Trajout_Single PDB;
    PDB.SetDebug( debug_ );
    if (PDB.InitTrajWrite( pdbout, "topresnum " + pdb_ter_arg, State.DSL(), TrajectoryFile::PDBFILE)) {
      mprinterr("Error: Could not initialize output PDB\n");
      return CpptrajState::ERR;
    }
    if (PDB.SetupTrajWrite(outCoords->TopPtr(), outCoords->CoordsInfo(), 1)) {
      mprinterr("Error: Could not set up output PDB\n");
      return CpptrajState::ERR;
    }
    PDB.PrintInfo(1);
    PDB.WriteSingle(0, frameIn);
    PDB.EndTraj();
  }

  outfile->CloseFile();

  if (fatal_errors > 0) {
    if (errorsAreFatal_) {
      mprinterr("Error: %i errors were encountered that will prevent LEaP from running successfully.\n", fatal_errors);
      return CpptrajState::ERR;
    } else {
      mprintf("Warning: %i errors were encountered that will prevent LEaP from running successfully.\n", fatal_errors);
      mprintf("Warning: Continuing on anyway, but final structure **NEEDS VALIDATION**.\n");
    }
  }
  // Run leap if needed
  if (!leapffname.empty()) {
    if (RunLeap( leapffname, leapfilename )) {
      mprinterr("Error: Running leap failed.\n");
      return CpptrajState::ERR;
    }
  }

  return CpptrajState::OK;
}
