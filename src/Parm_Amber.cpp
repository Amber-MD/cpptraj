#include <cstdio> // sscanf
#include <cstring> // strlen, strncmp
#include <ctime> // time for parm write
#include <locale> // isspace
#include <cstdlib> // atoi, atof
#include "Parm_Amber.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // RemoveTrailingWhitespace, SetIntegerFormatString etc
#include "Box.h"
#include "Constants.h" // ELECTOAMBER, AMBERTOELEC

// ---------- Constants and Enumerated types -----------------------------------
/// Enumerated type for FLAG_POINTERS section
const int Parm_Amber::AMBERPOINTERS=31;
enum topValues {
//0         1       2      3       4       5       6       7      8       9
  NATOM=0,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM,
  NNB,      NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,
  IFPERT,   NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS,  IFCAP,
  NEXTRA
};
/// Number of unique amber parm FLAGs
const int Parm_Amber::NUMAMBERPARMFLAGS=44;
/// Constant strings for fortran formats corresponding to Amber parm flags
const char* Parm_Amber::AmberParmFmt[NUMAMBERPARMFLAGS] = {
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(3I8)",
"%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(20a4)",   "%FORMAT(1a80)"
};
/// Constant strings for Amber parm flags
const char* Parm_Amber::AmberParmFlag[NUMAMBERPARMFLAGS] = {
  "POINTERS",
  "ATOM_NAME",
  "CHARGE",
  "MASS",
  "RESIDUE_LABEL",
  "RESIDUE_POINTER",
  "AMBER_ATOM_TYPE",
  "BONDS_INC_HYDROGEN",
  "BONDS_WITHOUT_HYDROGEN",
  "SOLVENT_POINTERS",
  "ATOMS_PER_MOLECULE",
  "BOX_DIMENSIONS",
  "ATOM_TYPE_INDEX",
  "NUMBER_EXCLUDED_ATOMS",
  "NONBONDED_PARM_INDEX",
  "LENNARD_JONES_ACOEF",
  "LENNARD_JONES_BCOEF",
  "EXCLUDED_ATOMS_LIST",
  "RADII",
  "SCREEN",
  "BOND_FORCE_CONSTANT",
  "BOND_EQUIL_VALUE",
  "ANGLE_FORCE_CONSTANT",
  "ANGLE_EQUIL_VALUE",
  "DIHEDRAL_FORCE_CONSTANT",
  "DIHEDRAL_PERIODICITY",
  "DIHEDRAL_PHASE",
  "SCEE_SCALE_FACTOR",
  "SCNB_SCALE_FACTOR",
  "SOLTY",
  "ANGLES_INC_HYDROGEN",
  "ANGLES_WITHOUT_HYDROGEN",
  "DIHEDRALS_INC_HYDROGEN",
  "DIHEDRALS_WITHOUT_HYDROGEN",
  "HBOND_ACOEF",
  "HBOND_BCOEF",
  "HBCUT",
  "TREE_CHAIN_CLASSIFICATION",
  "JOIN_ARRAY",
  "IROTAT",
  "ATOMIC_NUMBER",
  "TITLE",
  "CTITLE",
  "RADIUS_SET"
};

// -----------------------------------------------------------------------------
// CONSTRUCTOR
Parm_Amber::Parm_Amber() :
  debug_(0), 
  newParm_(false),
  ftype_(UNKNOWN_FTYPE),
  fncols_(0),
  fprecision_(0),
  fwidth_(0),
  error_count_(0),
  buffer_(0),
  buffer_size_(0),
  buffer_max_size_(0)
{ }

// DESTRUCTOR
Parm_Amber::~Parm_Amber() {
  if (buffer_!=0) delete[] buffer_;
}

// Parm_Amber::ID_ParmFormat()
bool Parm_Amber::ID_ParmFormat(CpptrajFile& fileIn) {
  int iamber[12];
  char lineBuf[BUF_SIZE];
  // Assumes already set up for READ
  if (fileIn.OpenFile()) return false;
  fileIn.Gets(lineBuf, BUF_SIZE);
  // Check for %VERSION
  if (strncmp(lineBuf,"%VERSION",8)==0) {
    fileIn.Gets(lineBuf, BUF_SIZE);
    // Check for %FLAG
    if (strncmp(lineBuf,"%FLAG",5)==0) {
      if (debug_>0) mprintf("  AMBER TOPOLOGY file\n");
      newParm_ = true;
      fileIn.CloseFile();
      return true;
    }
  } else {
    // Since %VERSION not present, If first line is 81 bytes and the second 
    // line has 12 numbers in 12I6 format, assume old-style Amber topology
    // NOTE: Could also be less than 81? Only look for 12 numbers?
    int line1size = (int)strlen(lineBuf);
    if (line1size == (81 + fileIn.IsDos())) {
      fileIn.Gets(lineBuf, BUF_SIZE);
      if ( sscanf(lineBuf,"%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i", 
                  iamber,   iamber+1, iamber+2,  iamber+3, 
                  iamber+4, iamber+5, iamber+6,  iamber+7, 
                  iamber+8, iamber+9, iamber+10, iamber+11) == 12 )
      {
        if (debug_>0) mprintf("  AMBER TOPOLOGY, OLD FORMAT\n");
        fileIn.CloseFile();
        return true;
      }
    }
  }
  fileIn.CloseFile();
  return false;
}

// Parm_Amber::ReadParm()
int Parm_Amber::ReadParm(std::string const& fname, Topology &TopIn ) {
  if (file_.OpenRead( fname )) return 1;
  int err = ReadParmAmber( TopIn );
  file_.CloseFile();
  return err;
}

// Parm_Amber::AmberIfbox()
/** Return Amber IFBOX type:
  *   0: No box
  *   1: Box
  *   2: Truncated octahedral box
  */
int Parm_Amber::AmberIfbox(const Box& boxIn) {
  if (boxIn.Type() == Box::NOBOX) return 0;
  else if (boxIn.Type() == Box::TRUNCOCT) return 2;
  return 1;
}

void Parm_Amber::CheckNameWidth(const char* typeIn, NameType const& nameIn) {
  if (nameIn[4] != '\0')
    mprintf("Warning: Parm_Amber: %s name (%s) is too large and will be truncated (4 chars max).\n",
            typeIn, *nameIn);
}

// Parm_Amber::WriteParm()
int Parm_Amber::WriteParm(std::string const& fname, Topology const& parmIn) {
  // For date and time
  time_t rawtime;
  struct tm *timeinfo;

  // Create arrays of atom info
  std::vector<NameType> names;
  std::vector<double> charge;
  std::vector<int> at_num;
  std::vector<double> mass;
  std::vector<int> atype_index;
  std::vector<NameType> types;
  std::vector<double> gb_radii;
  std::vector<double> gb_screen;
  std::vector<int> numex;
  std::vector<int> excluded;
  for (Topology::atom_iterator atom = parmIn.begin(); atom != parmIn.end(); atom++) 
  {
    names.push_back( (*atom).Name() );
    CheckNameWidth("Atom",names.back());
    charge.push_back( (*atom).Charge() * ELECTOAMBER );
    at_num.push_back( (*atom).AtomicNumber() );
    mass.push_back( (*atom).Mass() );
    atype_index.push_back( (*atom).TypeIndex() );
    types.push_back( (*atom).Type() );
    CheckNameWidth("Type",types.back());
    gb_radii.push_back( (*atom).Radius() );
    gb_screen.push_back( (*atom).Screen() );
    // Amber atom exclusion list prints a 0 placeholder for atoms with
    // no exclusions, so always print at least 1 for numex
    int nex = (*atom).Nexcluded();
    if (nex == 0) {
      numex.push_back( 1 );
      excluded.push_back( 0 );
    } else {
      numex.push_back( nex );
      for (Atom::excluded_iterator ex = (*atom).excludedbegin();
                                   ex != (*atom).excludedend(); ex++)
        // Amber atom #s start from 1
        excluded.push_back( (*ex) + 1 );
    }
  }

  // Create arrays of residue info
  std::vector<int> resnums;
  std::vector<NameType> resnames;
  for (Topology::res_iterator res = parmIn.ResStart(); res != parmIn.ResEnd(); res++)
  {
    resnames.push_back( (*res).Name() );
    // Amber atom #s start from 1
    resnums.push_back( (*res).FirstAtom()+1 );
  }

  // Get bond information. If atom bonding info is present but bond
  // parameters are not, create placeholder since other programs
  // may expect bond parameters if bonds/bondsh arrays defined.
  std::vector<int> Bonds = parmIn.Bonds();
  std::vector<int> BondsH = parmIn.BondsH();
  std::vector<double> BondRk = parmIn.BondRk();
  std::vector<double> BondReq = parmIn.BondReq();
  if ( (!Bonds.empty() || !BondsH.empty()) &&
       (BondRk.empty() && BondReq.empty())    )
  {
    mprintf("Warning: [%s] Bond information present but no bond parameters.\n",
            file_.Filename().base());
    mprintf("Warning: This can occur e.g. when bonds are determined from PDB.\n");
    mprintf("Warning: Dummy parameters will be created as placeholders.\n");
    mprintf("Warning: DO NOT USE THIS AMBER TOPOLOGY FOR SIMULATIONS!\n");
    BondRk.push_back(1.0);
    BondReq.push_back(1.0);
    for (unsigned int idx = 2; idx < Bonds.size(); idx += 3)
      Bonds[idx] = 1;
    for (unsigned int idx = 2; idx < BondsH.size(); idx += 3)
      BondsH[idx] = 1;
  }

  // Create pointer array
  std::vector<int> values(AMBERPOINTERS, 0);
  values[NATOM] = parmIn.Natom();
  values[NTYPES] = parmIn.Ntypes();
  values[NBONH] = (int)BondsH.size() / 3; // NOTE: Check divisible by 3?
  values[MBONA] = (int)Bonds.size() / 3; // NOTE: Check divisible by 3?
  values[NTHETH] = (int)parmIn.AnglesH().size() / 4;
  values[MTHETA] = (int)parmIn.Angles().size() / 4;
  values[NPHIH] = (int)parmIn.DihedralsH().size() / 5; 
  values[MPHIA] = (int)parmIn.Dihedrals().size() / 5;
  values[NNB] = (int)excluded.size();
  values[NRES] = parmIn.Nres();
  //   NOTE: Assuming NBONA == MBONA etc
  values[NBONA] = values[MBONA];
  values[NTHETA] = values[MTHETA];
  values[NPHIA] = values[MPHIA];
  values[NUMBND] = (int)BondRk.size();
  values[NUMANG] = (int)parmIn.AngleTk().size();
  values[NPTRA] = (int)parmIn.DihedralPk().size();
  values[NATYP] = (int)parmIn.Solty().size(); // Only for SOLTY
  values[NPHB] = (int)parmIn.Asol().size();
  values[IFBOX] = AmberIfbox( parmIn.ParmBox() );
  values[NMXRS] = parmIn.FindResidueMaxNatom();
    
  // Write parm
  if (file_.OpenWrite( fname )) return 1;
  // HEADER AND TITLE (4 lines, version, flag, format, title)
  time( &rawtime );
  timeinfo = localtime(&rawtime);
  file_.Printf("%-44s%02i/%02i/%02i  %02i:%02i:%02i                  \n",
               "%VERSION  VERSION_STAMP = V0001.000  DATE = ",
               timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year%100,
               timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  std::string title = parmIn.ParmName();
  // Resize title to max 80 char
  if (title.size() > 80)
    title.resize(80);
  file_.Printf("%-80s\n%-80s\n%-80s\n","%FLAG TITLE","%FORMAT(20a4)",title.c_str());
  // POINTERS
  WriteInteger(F_POINTERS, values);
  // NAMES, CHARGE, ATOMIC NUMBER, MASS, TYPE INDEX, EXCLUDED, NB INDEX
  WriteName(F_NAMES, names);
  WriteDouble(F_CHARGE, charge);
  WriteInteger(F_ATOMICNUM, at_num);
  WriteDouble(F_MASS, mass);
  WriteInteger(F_ATYPEIDX, atype_index);
  WriteInteger(F_NUMEX, numex);
  WriteInteger(F_NB_INDEX, parmIn.NB_index());
  // RESIDUE NAME, POINTER
  WriteName(F_RESNAMES, resnames);
  WriteInteger(F_RESNUMS, resnums);
  // BOND, ANGLE, and DIHEDRAL FORCE CONSTANT and EQUIL VALUES
  WriteDouble(F_BONDRK, BondRk);
  WriteDouble(F_BONDREQ, BondReq);
  WriteDouble(F_ANGLETK, parmIn.AngleTk());
  WriteDouble(F_ANGLETEQ, parmIn.AngleTeq());
  WriteDouble(F_DIHPK, parmIn.DihedralPk());
  WriteDouble(F_DIHPN, parmIn.DihedralPn());
  WriteDouble(F_DIHPHASE, parmIn.DihedralPhase());
  WriteDouble(F_SCEE, parmIn.SCEE());
  WriteDouble(F_SCNB, parmIn.SCNB());
  // SOLTY - Currently unused
  WriteDouble(F_SOLTY, parmIn.Solty());
  // LJ params
  WriteDouble(F_LJ_A, parmIn.LJA());
  WriteDouble(F_LJ_B, parmIn.LJB());
  // BONDS/ANGLES/DIHEDRAL INDICES WITH AND WITHOUT HYDROGEN
  WriteInteger(F_BONDSH, BondsH); 
  WriteInteger(F_BONDS, Bonds);
  WriteInteger(F_ANGLESH, parmIn.AnglesH());
  WriteInteger(F_ANGLES, parmIn.Angles());
  WriteInteger(F_DIHH, parmIn.DihedralsH());
  WriteInteger(F_DIH, parmIn.Dihedrals());
  // EXCLUDED ATOMS LIST
  WriteInteger(F_EXCLUDE, excluded);
  // HBOND
  WriteDouble(F_ASOL, parmIn.Asol());
  WriteDouble(F_BSOL, parmIn.Bsol());
  WriteDouble(F_HBCUT, parmIn.HBcut());
  // AMBER ATOM TYPE
  WriteName(F_TYPES, types);
  // TREE CHAIN CLASSIFICATION, JOIN, IROTAT
  // TODO: Generate automatically
  WriteName(F_ITREE, parmIn.Itree());
  WriteInteger(F_JOIN, parmIn.Join());
  WriteInteger(F_IROTAT, parmIn.Irotat());
  // Write solvent info if IFBOX>0
  if (values[IFBOX] > 0) {
    // Determine first solvent residue
    int firstSolventMol = -1;
    for (Topology::mol_iterator mol = parmIn.MolStart(); mol != parmIn.MolEnd(); ++mol) {
      if ( (*mol).IsSolvent() ) { 
        firstSolventMol = (int)(mol - parmIn.MolStart() + 1); 
        break;
      }
    }
    // If no solvent, just set to 1 beyond # of molecules
    if (firstSolventMol == -1)
      firstSolventMol = parmIn.Nmol() + 1;
    // Solvent Pointers
    std::vector<int> solvent_pointer(3);
    solvent_pointer[0] = parmIn.FinalSoluteRes(); // Already +1
    solvent_pointer[1] = parmIn.Nmol();
    solvent_pointer[2] = firstSolventMol;
    if (solvent_pointer[2] == 0)
      solvent_pointer[2] = solvent_pointer[1] + 1;
    WriteInteger(F_SOLVENT_POINTER, solvent_pointer);
    // ATOMS PER MOLECULE
    std::vector<int> APM;
    APM.reserve( solvent_pointer[1] );
    for (Topology::mol_iterator mol = parmIn.MolStart(); mol != parmIn.MolEnd(); mol++)
      APM.push_back( (*mol).NumAtoms() );
    WriteInteger(F_ATOMSPERMOL, APM);
    // BOX DIMENSIONS
    std::vector<double> betaLengths(4); // Beta X Y Z
    betaLengths[0] = parmIn.ParmBox().Beta();
    betaLengths[1] = parmIn.ParmBox().BoxX();
    betaLengths[2] = parmIn.ParmBox().BoxY();
    betaLengths[3] = parmIn.ParmBox().BoxZ();
    WriteDouble(F_PARMBOX, betaLengths);
  }
  // GB RADIUS SET, RADII, SCREENING PARAMETERS
  std::string radius_set = parmIn.GBradiiSet();
  if (!radius_set.empty()) {
    WriteSetup(F_RADSET, 1);
    if (radius_set.size()>80)
      radius_set.resize(80);
    file_.Printf("%-80s\n",radius_set.c_str());
  }
  WriteDouble(F_RADII, gb_radii);
  WriteDouble(F_SCREEN, gb_screen);
 
  file_.CloseFile();

  return 0;
}

// -----------------------------------------------------------------------------
// Parm_Amber::ReadParmAmber()
int Parm_Amber::ReadParmAmber( Topology &TopIn ) {
  std::vector<int> atomsPerMol;
  //int finalSoluteRes = -1;
  //int firstSolvMol = -1;
  Box parmbox;
  bool chamber = false; // True if this top is a chamber-created top file.
  std::string title;
  int Npointers = AMBERPOINTERS; 

  if (newParm_) {
    if (debug_>0) 
      mprintf("\tReading Amber Topology file %s\n",file_.Filename().base());
    // Title. If not found check for CTITLE (chamber)
    if (PositionFileAtFlag(F_TITLE)) {
      title = GetLine();
    } else {
      if (PositionFileAtFlag(F_CTITLE)) {
        title = GetLine();
        chamber = true;
      } else {
        // No TITLE or CTITLE, weird, but dont crash out yet.
        mprintf("Warning: [%s] No TITLE in Amber Parm.\n",file_.Filename().base());
      }
    }
  } else {
    mprintf("\tReading old (<v7) Amber Topology file.\n");
    title = GetLine();
    // One less POINTERS (no NEXTRA)
    Npointers = 30;
  }
  mprintf("\tAmberParm Title: [%s]\n",title.c_str());
  TopIn.SetParmName( title, file_.Filename().Base() );
  // POINTERS
  std::vector<int> values = GetFlagInteger(F_POINTERS, Npointers);
  if (values.empty()) {
    mprinterr("Error: [%s] Could not get POINTERS from Amber Topology.\n",file_.Filename().base());
    return 1;
  }
  // DEBUG
  //for (std::vector<int>::iterator val = values.begin(); val != values.end(); ++val)
  //  mprintf("\t%i\n", *val);

  // Read parm variables
  std::vector<NameType> names = GetFlagName(F_NAMES, values[NATOM]);
  std::vector<double> charge = GetFlagDouble(F_CHARGE,values[NATOM]);
  // New parm only: ATOMICNUM
  std::vector<int> at_num;
  if (newParm_)
    at_num = GetFlagInteger(F_ATOMICNUM,values[NATOM]);
  std::vector<double> mass = GetFlagDouble(F_MASS,values[NATOM]);
  std::vector<int> atype_index = GetFlagInteger(F_ATYPEIDX,values[NATOM]);
  // For old parm need to read past NUMEX
  if (!newParm_)
    GetFlagInteger(F_NUMEX,values[NATOM]);
  std::vector<int> NB_index = GetFlagInteger(F_NB_INDEX,values[NTYPES]*values[NTYPES]);
  std::vector<NameType> resnames = GetFlagName(F_RESNAMES,values[NRES]);
  std::vector<int> resnums = GetFlagInteger(F_RESNUMS,values[NRES]);
  std::vector<double> bond_rk = GetFlagDouble(F_BONDRK, values[NUMBND]);
  std::vector<double> bond_req = GetFlagDouble(F_BONDREQ, values[NUMBND]);
  std::vector<double> angle_tk = GetFlagDouble(F_ANGLETK, values[NUMANG]);
  std::vector<double> angle_teq = GetFlagDouble(F_ANGLETEQ, values[NUMANG]);
  std::vector<double> dihedral_pk = GetFlagDouble(F_DIHPK, values[NPTRA]);
  std::vector<double> dihedral_pn = GetFlagDouble(F_DIHPN, values[NPTRA]);
  std::vector<double> dihedral_phase = GetFlagDouble(F_DIHPHASE, values[NPTRA]);
  // New parm only: SCEE and SCNB
  std::vector<double> scee_scale;
  std::vector<double> scnb_scale;
  if (newParm_) {
    scee_scale = GetFlagDouble(F_SCEE,values[NPTRA]);
    scnb_scale = GetFlagDouble(F_SCNB,values[NPTRA]);
  }
  // SOLTY: currently unused
  std::vector<double> solty = GetFlagDouble(F_SOLTY, values[NATYP]);
  int nlj = values[NTYPES] * (values[NTYPES]+1) / 2;
  std::vector<double> LJ_A = GetFlagDouble(F_LJ_A,nlj);
  std::vector<double> LJ_B = GetFlagDouble(F_LJ_B,nlj);
  std::vector<int> bondsh = GetFlagInteger(F_BONDSH,values[NBONH]*3);
  std::vector<int> bonds = GetFlagInteger(F_BONDS,values[NBONA]*3);
  std::vector<int> anglesh = GetFlagInteger(F_ANGLESH, values[NTHETH]*4);
  std::vector<int> angles = GetFlagInteger(F_ANGLES, values[NTHETA]*4);
  std::vector<int> dihedralsh = GetFlagInteger(F_DIHH, values[NPHIH]*5);
  std::vector<int> dihedrals = GetFlagInteger(F_DIH, values[NPHIA]*5);
  // For old parm need to read past EXCLUDE
  if (!newParm_)
    GetFlagInteger(F_EXCLUDE, values[NNB]);
  std::vector<double> asol = GetFlagDouble(F_ASOL,values[NPHB]);
  std::vector<double> bsol = GetFlagDouble(F_BSOL,values[NPHB]);
  std::vector<double> hbcut = GetFlagDouble(F_HBCUT,values[NPHB]);
  std::vector<NameType> types = GetFlagName(F_TYPES,values[NATOM]);
  std::vector<NameType> itree = GetFlagName(F_ITREE,values[NATOM]);
  std::vector<int> join_array = GetFlagInteger(F_JOIN,values[NATOM]);
  std::vector<int> irotat = GetFlagInteger(F_IROTAT,values[NATOM]);
  // Get solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    std::vector<int> solvent_pointer = GetFlagInteger(F_SOLVENT_POINTER,3);
    if (solvent_pointer.empty()) {
      mprintf("Could not get %s from Amber Topology file.\n",AmberParmFlag[F_SOLVENT_POINTER]);
      return 1;
    }
    //finalSoluteRes = solvent_pointer[0] - 1;
    int molecules = solvent_pointer[1];
    //firstSolvMol = solvent_pointer[2] - 1;
    atomsPerMol = GetFlagInteger(F_ATOMSPERMOL,molecules);
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    std::vector<double> boxFromParm = GetFlagDouble(F_PARMBOX,4);
    // If no box information present in the parm (such as with Chamber prmtops)
    // set the box info if ifbox = 2, otherwise default is NOBOX; the box info 
    // will eventually be set by angles from the first trajectory associated 
    // with this parm.
    if (boxFromParm.empty()) {
      if (not chamber) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (values[IFBOX] == 2) 
        parmbox.SetTruncOct(); 
    // Determine box type, set Box angles and lengths from beta (boxFromParm[0])
    } else {
      parmbox.SetBetaLengths( boxFromParm[0], boxFromParm[1], boxFromParm[2], boxFromParm[3] );
    }
    // Check for IFBOX/BoxType mismatch
    if (values[IFBOX]==2 && parmbox.Type() != Box::TRUNCOCT) {
      mprintf("Warning: Amber Parm Box should be Truncated Octahedron (ifbox==2)\n");
      mprintf("         but BOX_DIMENSIONS indicate %s - may cause imaging problems.\n",
              parmbox.TypeName());
    }
  }
  // New Parm only: GB parameters; radius set, radii, and screening parameters
  std::vector<double> gb_radii;
  std::vector<double> gb_screen;
  if (newParm_) {
    if (PositionFileAtFlag(F_RADSET)) {
      std::string radius_set = GetLine();
      mprintf("\tRadius Set: %s\n",radius_set.c_str());
      TopIn.SetGBradiiSet( radius_set );
    }
    gb_radii = GetFlagDouble(F_RADII,values[NATOM]);
    gb_screen = GetFlagDouble(F_SCREEN,values[NATOM]); 
  }
 
  if (error_count_==0) {
    // Convert Amber Charge to elec
    for (std::vector<double>::iterator q = charge.begin(); q != charge.end(); q++)
      *q *= (AMBERTOELEC);
    // Shift atom #s in resnums by -1 so they start from 0
    for (std::vector<int>::iterator r = resnums.begin(); r != resnums.end(); r++)
      --(*r);
    // If ever used, shift atom #s in excludedAtoms by -1 so they start from 0 
    error_count_ += TopIn.CreateAtomArray( names, charge, at_num, mass, atype_index, 
                                           types, gb_radii, gb_screen,
                                           resnames, resnums );
    // Shift indices in parm arrays?
    error_count_ += TopIn.SetBondInfo(bonds, bondsh, bond_rk, bond_req);
    error_count_ += TopIn.SetAngleInfo(angles, anglesh, angle_tk, angle_teq);
    error_count_ += TopIn.SetDihedralInfo(dihedrals, dihedralsh, dihedral_pk,
                                          dihedral_pn, dihedral_phase,
                                          scee_scale, scnb_scale);
    error_count_ += TopIn.SetAmberHbond(asol, bsol, hbcut);
    error_count_ += TopIn.SetAmberExtra(solty, itree, join_array, irotat);
    error_count_ += TopIn.SetNonbondInfo(values[NTYPES], NB_index, LJ_A, LJ_B);
    if (values[IFBOX]>0) 
      TopIn.SetBox( parmbox );
  }

  return error_count_;
}

// -----------------------------------------------------------------------------
// Parm_Amber::GetLine()
std::string Parm_Amber::GetLine() {
  file_.Gets(lineBuffer_, BUF_SIZE);
  std::string line( lineBuffer_ );
  // Remove any trailing whitespace
  RemoveTrailingWhitespace( line );
  return line;
}

// Parm_Amber::GetInteger()
std::vector<int> Parm_Amber::GetInteger(int width, int ncols, int maxval) {
  std::vector<int> iarray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return iarray;
  else if (err == -1) {
    mprinterr("Error in read of integer values from %s\n",file_.Filename().base());
    ++error_count_;
    return iarray;
  }
  // Reserve variable memory
  iarray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    iarray.push_back( atoi(ptrbegin) );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return iarray;
}

// Parm_Amber::GetDouble()
std::vector<double> Parm_Amber::GetDouble(int width, int ncols, int maxval) {
  std::vector<double> darray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return darray;
  else if (err == -1) {
    mprinterr("Error in read of double values from %s\n",file_.Filename().base());
    ++error_count_;
    return darray;
  }
  // Reserve variable memory
  darray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    darray.push_back( atof(ptrbegin) );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return darray;
}

// Parm_Amber::GetName()
std::vector<NameType> Parm_Amber::GetName(int width, int ncols, int maxval) {
  std::vector<NameType> carray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return carray;
  else if (err == -1) {
    mprinterr("Error in read of Name values from %s\n",file_.Filename().base());
    ++error_count_;
    return carray;
  }
  // Reserve variable memory
  carray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    carray.push_back( ptrbegin );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return carray;
}

// Parm_Amber::GetFlagInteger()
std::vector<int> Parm_Amber::GetFlagInteger(AmberParmFlagType fflag, int maxval) {
  std::vector<int> iarray;
  if (newParm_) {
    // Seek to flag and set up fncols, fwidth
    if (!SeekToFlag(fflag)) return iarray;
    // NOTE: Check that type matches?
    iarray = GetInteger( fwidth_, fncols_, maxval );
  } else {
    iarray = GetInteger( 6, 12, maxval );
  }
  if ((int)iarray.size() != maxval)
    mprinterr("Error: Reading INTEGER section %s\n",AmberParmFlag[fflag]);
  return iarray;
}

// Parm_Amber::GetFlagDouble()
std::vector<double> Parm_Amber::GetFlagDouble(AmberParmFlagType fflag, int maxval) {
  std::vector<double> darray;
  if (newParm_) {
    // Seek to flag and set up fncols, fwidth
    if (!SeekToFlag(fflag)) return darray;
    // NOTE: Check that type matches?
    darray = GetDouble( fwidth_, fncols_, maxval );
  } else {
    darray = GetDouble( 16, 5, maxval );
  }
  if ((int)darray.size() != maxval)
    mprinterr("Error: Reading DOUBLE section %s\n",AmberParmFlag[fflag]); 
  return darray;
}

// Parm_Amber::GetFlagName()
std::vector<NameType> Parm_Amber::GetFlagName(AmberParmFlagType fflag, int maxval) {
  std::vector<NameType> carray;
  if (newParm_) {
    // Seek to flag and set up fncols, fwidth
    if (!SeekToFlag(fflag)) return carray;
    // NOTE: Check that type matches?
    carray = GetName( fwidth_, fncols_, maxval );
  } else {
    carray = GetName( 4, 20, maxval );
  }
  if ((int)carray.size() != maxval)
    mprinterr("Error: Reading CHAR section %s\n",AmberParmFlag[fflag]);
  return carray;
}

// Parm_Amber::SeekToFlag()
bool Parm_Amber::SeekToFlag(AmberParmFlagType fflag) {
  // Find flag, set format line
  if (!PositionFileAtFlag(fflag)) return false;
  // Set up cols, width, etc from format
  if (!SetFortranType()) return false;
  return true;
}  

// Parm_Amber::AllocateAndRead()
int Parm_Amber::AllocateAndRead(int width, int ncols, int maxval) {
  char temp[3]; // Only for when maxval is 0, space for \n, \r, null
  // Sanity Check
  if (maxval < 0) {
    mprinterr("Internal Error: AllocateAndRead called with maxval < 0 (%i)\n",maxval);
    return -1;
  }
  // If # expected values is 0 there will still be a newline placeholder
  // in the parmtop. Read past that and return
  if (maxval==0) {
    file_.Gets(temp,2);
    return 0;
  }
  // Allocate buffer to read in entire section
  size_t BufferSize = GetFortranBufferSize(width, ncols, maxval);
  if (buffer_!=0) delete[] buffer_;
  buffer_ = new char[ BufferSize ];
  // Read section from file
  int err = file_.Read(buffer_,BufferSize);
  return err;
}

// Parm_Amber::PositionFileAtFlag()
bool Parm_Amber::PositionFileAtFlag(AmberParmFlagType flag) {
  const char *Key = AmberParmFlag[flag];
  char value[BUF_SIZE];
  bool searchFile = true;
  bool hasLooped = false;

  if (debug_ > 0) mprintf("Reading %s\n",Key);
  // Search for '%FLAG <Key>'
  while ( searchFile ) {
    while ( file_.Gets(lineBuffer_, BUF_SIZE) == 0 ) {
      if ( strncmp(lineBuffer_,"%FLAG",5)==0 ) {
        // Check flag Key
        sscanf(lineBuffer_, "%*s %s",value);
        if (strcmp(value,Key)==0) {
          if (debug_>1) mprintf("DEBUG: Found Flag Key [%s]\n",value);
          // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
          // read past until you get to the FORMAT line
          file_.Gets(lineBuffer_, BUF_SIZE); 
          while (strncmp(lineBuffer_, "%FORMAT",7)!=0)
            file_.Gets(lineBuffer_, BUF_SIZE);
          if (debug_>1) mprintf("DEBUG: Format line [%s]\n",lineBuffer_);
          // Set format
          // NOTE: Check against stored formats?
          fformat_.assign(lineBuffer_);
          return true;
        } // END found Key
      } // END found FLAG line
    } // END scan through file
    // If we havent yet tried to search from the beginning, try it now.
    // Otherwise the Key has not been found.
    if (!hasLooped) {
      file_.Rewind();
      hasLooped = true;
    } else
      searchFile = false;
  }

  // If we reached here Key was not found.
  if (debug_>0)
    mprintf("Warning: [%s] Could not find Key %s in file.\n",file_.Filename().base(),Key);
  fformat_.clear();
  return false;
}

// -----------------------------------------------------------------------------
// Parm_Amber::WriteSetup()
int Parm_Amber::WriteSetup(AmberParmFlagType fflag, size_t Nelements) {
  // Assign format string
  fformat_.assign( AmberParmFmt[fflag] );
  //mprintf("DEBUG: FlagFormat[%s], Nelements=%zu\n",fformat_.c_str(),Nelements);
  // Set type, cols, width, and precision from format string
  if (!SetFortranType()) return 1;
  // Write FLAG and FORMAT lines
  file_.Printf("%%FLAG %-74s\n%-80s\n",AmberParmFlag[fflag],AmberParmFmt[fflag]);
  // If Nelements is 0 just print newline and exit
  if (Nelements == 0) {
    file_.Printf("\n");
    return 1;
  }
  // Allocate character buffer space for memory, +1 for null
  size_t bufsize = GetFortranBufferSize(fwidth_, fncols_, Nelements);
  //mprinterr("DEBUG: Current max=%zu  Current size=%zu  Requested size=%zu\n",
  //          buffer_max_size_, buffer_size_, bufsize);
  if (bufsize+1 > buffer_max_size_) {
    if (buffer_!=0) delete[] buffer_;
    buffer_ = new char[ bufsize+1 ];
    buffer_max_size_ = bufsize+1;
  }
  buffer_size_ = bufsize;
  return 0;
}

// Parm_Amber::WriteInteger()
int Parm_Amber::WriteInteger(AmberParmFlagType fflag, std::vector<int>const& iarray)
{
  std::string FS;
  if (WriteSetup(fflag, iarray.size())) return 0;
  // Set up printf format string; true == no leading space
  SetIntegerFormatString(FS, fwidth_, true);
  const char *FORMAT = FS.c_str();
  char *ptr = buffer_;
  int col = 0;
  for (std::vector<int>::const_iterator it = iarray.begin(); it != iarray.end(); it++) {
    sprintf(ptr, FORMAT, *it);
    // If # chars written > width, this silently truncates
    ptr += fwidth_;
    ++col;
    if (col == fncols_) {
      sprintf(ptr,"\n");
      ++ptr;
      col = 0;
    }
  }
  //mprinterr("INT: Last col written = %i (%i)\n",col,fncols_);
  if (col != 0) sprintf(ptr,"\n");
  file_.Write(buffer_, buffer_size_);

  return 0;
}

// Parm_Amber::WriteDouble()
int Parm_Amber::WriteDouble(AmberParmFlagType fflag, std::vector<double>const& darray)
{
  std::string FS;
  if (WriteSetup(fflag, darray.size())) return 0;
  // Set up printf format string; true == no leading space
  SetDoubleFormatString(FS, fwidth_, fprecision_, 2, true);
  const char *FORMAT = FS.c_str();
  char *ptr = buffer_;
  int col = 0;
  for (std::vector<double>::const_iterator it = darray.begin(); it != darray.end(); it++) {
    sprintf(ptr, FORMAT, *it);
    // If # chars written > width, this silently truncates
    ptr += fwidth_;
    ++col;
    if (col == fncols_) {
      sprintf(ptr,"\n");
      ++ptr;
      col = 0;
    }
  }
  //mprinterr("DOUBLE: Last col written = %i (%i)\n",col,fncols_);
  if (col != 0) sprintf(ptr,"\n");
  file_.Write(buffer_, buffer_size_);

  return 0;
}

// Parm_Amber::WriteName()
int Parm_Amber::WriteName(AmberParmFlagType fflag, std::vector<NameType>const& carray)
{
  std::string FS;
  if (WriteSetup(fflag, carray.size())) return 0;
  // Set up printf format string; true == no leading space
  SetStringFormatString(FS, fwidth_, true);
  const char *FORMAT = FS.c_str();
  char *ptr = buffer_;
  int col = 0;
  for (std::vector<NameType>::const_iterator it = carray.begin(); it != carray.end(); it++) {
    sprintf(ptr, FORMAT, *(*it));
    // If # chars written > width, this silently truncates
    ptr += fwidth_;
    ++col;
    if (col == fncols_) {
      sprintf(ptr,"\n");
      ++ptr;
      col = 0;
    }
  }
  //mprinterr("NAME: Last col written = %i (%i)\n",col,fncols_);
  if (col != 0) sprintf(ptr,"\n");
  file_.Write(buffer_, buffer_size_);

  return 0;
}

// -----------------------------------------------------------------------------
// Parm_Amber::GetFortranBufferSize()
/** Given number of columns and the width of each column, return the 
  * necessary char buffer size for N data elements.
  */
size_t Parm_Amber::GetFortranBufferSize(int width, int ncols, int N) {
  size_t bufferLines=0;
  size_t BufferSize=0;

  BufferSize = N * width;
  bufferLines = N / ncols;
  if ((N % ncols)!=0) ++bufferLines;
  // If DOS file there are CRs before Newlines
  if (file_.IsDos()) bufferLines *= 2;
  BufferSize += bufferLines;
  //if (debug_>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
  return BufferSize;
}

// Parm_Amber::SetFortranType()
/** Given a fortran-type format string, set the corresponding fortran
  * type. Set fncols (if present), fwidth, and fprecision (if present).
  */
// 01234567
// %FORMAT([<cols>][(]<type><width>[<precision>][)])
bool Parm_Amber::SetFortranType() {
  std::locale loc;
  std::string arg;

  if ( fformat_.empty() ) return false;
  // Make sure characters are upper case.
  for (std::string::iterator p = fformat_.begin(); p != fformat_.end(); p++)
    toupper(*p, loc);
  // Advance past left parentheses
  std::string::iterator ptr = fformat_.begin() + 7;
  while (*ptr=='(') ++ptr;
  // If digit, have number of data columns
  fncols_ = 0;
  if (isdigit(*ptr, loc)) {
    while (ptr!=fformat_.end() && isdigit(*ptr, loc)) {
      arg += *ptr;
      ++ptr;
    }
    fncols_ = atoi( arg.c_str() );
  }
  // Advance past any more left parentheses
  while (ptr!=fformat_.end() && *ptr=='(') ++ptr;
  // Type
  if (ptr==fformat_.end()) {
    mprinterr("Error: Malformed fortran format string (%s) in Amber Topology %s\n",
              fformat_.c_str(), file_.Filename().base());
    return false;
  }
  switch (*ptr) {
    case 'I' : ftype_ = FINT;    break;
    case 'E' : ftype_ = FDOUBLE; break;
    case 'A' : ftype_ = FCHAR;   break;
    case 'F' : ftype_ = FFLOAT;  break;
    default  : ftype_ = UNKNOWN_FTYPE;
  }
  ++ptr;
  // Width
  fwidth_ = 0;
  arg.clear(); 
  while (isdigit(*ptr,loc)) {
    arg += *ptr;
    ++ptr;
  }
  fwidth_ = atoi( arg.c_str() );
  // Precision
  fprecision_ = 0;
  if (*ptr == '.') {
    ++ptr;
    arg.clear();
    while (isdigit(*ptr,loc)) {
      arg += *ptr;
      ++ptr;
    }
    fprecision_ = atoi( arg.c_str() );
  }
  if (debug_ > 2)
    mprintf("[%s]: cols=%i type=%i width=%i precision=%i\n",fformat_.c_str(),
            fncols_,(int)ftype_,fwidth_,fprecision_);

  return true;
}
