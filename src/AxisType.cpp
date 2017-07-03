#include <map>
#include "AxisType.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "TorsionRoutines.h" // pucker calc
#include "Constants.h" // pucker calc
#include "PDBfile.h" // load base reference
// ---------- NA_Reference -----------------------------------------------------
/** Add all common variations of NA base name given single letter X:
  * DX/DX3/DX5 (not for U), RX/RX3/RX5 (not for T), X3/X5, X
  */
static inline void AddBaseNames(std::string const& rname, RefBase& base) {
  if (rname != "U") {
    base.AddName("D" + rname );
    base.AddName("D" + rname + "3" );
    base.AddName("D" + rname + "5" );
  }
  if (rname != "T") {
    base.AddName("R" + rname );
    base.AddName("R" + rname + "3" );
    base.AddName("R" + rname + "5" );
  }
  base.AddName(rname + "3");
  base.AddName(rname + "5");
  base.AddName(rname);
}

/// CONSTRUCTOR - load default base references.
NA_Reference::NA_Reference() {
  // NOTE: This is slightly more wasteful of memory since the reference
  //       instance is created every time NA_Reference is constructed
  //       but allows for the addition of custom residues.
  // ADE
  bases_.push_back( RefBase('A', "ADE", NA_Base::ADE) );
  RefBase& ADE = bases_.back();
  AddBaseNames("A", ADE);
  ADE.AddAtom(NA_Atom( -2.479000,  5.346000,  0.000000, NA_Base::NONE,     0, "C1' "));
  ADE.AddAtom(NA_Atom( -1.291000,  4.498000,  0.000000, NA_Base::NONE,     1, "N9  "));
  ADE.AddAtom(NA_Atom(  0.024000,  4.897000,  0.000000, NA_Base::NONE,     1, "C8  "));
  ADE.AddAtom(NA_Atom(  0.877000,  3.902000,  0.000000, NA_Base::ACCEPTOR, 1, "N7  "));
  ADE.AddAtom(NA_Atom(  0.071000,  2.771000,  0.000000, NA_Base::NONE,     1, "C5  "));
  ADE.AddAtom(NA_Atom(  0.369000,  1.398000,  0.000000, NA_Base::NONE,     1, "C6  "));
  ADE.AddAtom(NA_Atom(  1.611000,  0.909000,  0.000000, NA_Base::DONOR,    0, "N6  "));
  ADE.AddAtom(NA_Atom( -0.668000,  0.532000,  0.000000, NA_Base::ACCEPTOR, 1, "N1  "));
  ADE.AddAtom(NA_Atom( -1.912000,  1.023000,  0.000000, NA_Base::NONE,     1, "C2  "));
  ADE.AddAtom(NA_Atom( -2.320000,  2.290000,  0.000000, NA_Base::ACCEPTOR, 1, "N3  "));
  ADE.AddAtom(NA_Atom( -1.267000,  3.124000,  0.000000, NA_Base::NONE,     1, "C4  "));
  // CYT
  bases_.push_back( RefBase('C', "CYT", NA_Base::CYT) );
  RefBase& CYT = bases_.back();
  AddBaseNames("C", CYT);
  CYT.AddAtom(NA_Atom( -2.477000,  5.402000,  0.000000, NA_Base::NONE,     0, "C1' "));
  CYT.AddAtom(NA_Atom( -1.285000,  4.542000,  0.000000, NA_Base::NONE,     1, "N1  "));
  CYT.AddAtom(NA_Atom( -1.472000,  3.158000,  0.000000, NA_Base::NONE,     1, "C2  "));
  CYT.AddAtom(NA_Atom( -2.628000,  2.709000,  0.000000, NA_Base::ACCEPTOR, 0, "O2  "));
  CYT.AddAtom(NA_Atom( -0.391000,  2.344000,  0.000000, NA_Base::ACCEPTOR, 1, "N3  "));
  CYT.AddAtom(NA_Atom(  0.837000,  2.868000,  0.000000, NA_Base::NONE,     1, "C4  "));
  CYT.AddAtom(NA_Atom(  1.875000,  2.027000,  0.000000, NA_Base::DONOR,    0, "N4  "));
  CYT.AddAtom(NA_Atom(  1.056000,  4.275000,  0.000000, NA_Base::NONE,     1, "C5  "));
  CYT.AddAtom(NA_Atom( -0.023000,  5.068000,  0.000000, NA_Base::NONE,     1, "C6  "));
  // GUA
  bases_.push_back( RefBase('G', "GUA", NA_Base::GUA) );
  RefBase& GUA = bases_.back();
  AddBaseNames("G", GUA);
  GUA.AddAtom(NA_Atom( -2.477000,  5.399000,  0.000000, NA_Base::NONE,     0, "C1' "));
  GUA.AddAtom(NA_Atom( -1.289000,  4.551000,  0.000000, NA_Base::NONE,     1, "N9  "));
  GUA.AddAtom(NA_Atom(  0.023000,  4.962000,  0.000000, NA_Base::NONE,     1, "C8  "));
  GUA.AddAtom(NA_Atom(  0.870000,  3.969000,  0.000000, NA_Base::ACCEPTOR, 1, "N7  "));
  GUA.AddAtom(NA_Atom(  0.071000,  2.833000,  0.000000, NA_Base::NONE,     1, "C5  "));
  GUA.AddAtom(NA_Atom(  0.424000,  1.460000,  0.000000, NA_Base::NONE,     1, "C6  "));
  GUA.AddAtom(NA_Atom(  1.554000,  0.955000,  0.000000, NA_Base::ACCEPTOR, 0, "O6  "));
  GUA.AddAtom(NA_Atom( -0.700000,  0.641000,  0.000000, NA_Base::DONOR,    1, "N1  "));
  GUA.AddAtom(NA_Atom( -1.999000,  1.087000,  0.000000, NA_Base::NONE,     1, "C2  "));
  GUA.AddAtom(NA_Atom( -2.949000,  0.139000, -0.001000, NA_Base::DONOR,    0, "N2  "));
  GUA.AddAtom(NA_Atom( -2.342000,  2.364000,  0.001000, NA_Base::ACCEPTOR, 1, "N3  "));
  GUA.AddAtom(NA_Atom( -1.265000,  3.177000,  0.000000, NA_Base::NONE,     1, "C4  "));
  // THY
  bases_.push_back( RefBase('T', "THY", NA_Base::THY) );
  RefBase& THY = bases_.back();
  AddBaseNames("T", THY);
  THY.AddAtom(NA_Atom( -2.481000,  5.354000,  0.000000, NA_Base::NONE,     0, "C1' "));
  THY.AddAtom(NA_Atom( -1.284000,  4.500000,  0.000000, NA_Base::NONE,     1, "N1  "));
  THY.AddAtom(NA_Atom( -1.462000,  3.135000,  0.000000, NA_Base::NONE,     1, "C2  "));
  THY.AddAtom(NA_Atom( -2.562000,  2.608000,  0.000000, NA_Base::ACCEPTOR, 0, "O2  "));
  THY.AddAtom(NA_Atom( -0.298000,  2.407000,  0.000000, NA_Base::DONOR,    1, "N3  "));
  THY.AddAtom(NA_Atom(  0.994000,  2.897000,  0.000000, NA_Base::NONE,     1, "C4  "));
  THY.AddAtom(NA_Atom(  1.944000,  2.119000,  0.000000, NA_Base::ACCEPTOR, 0, "O4  "));
  THY.AddAtom(NA_Atom(  1.106000,  4.338000,  0.000000, NA_Base::NONE,     1, "C5  "));
  THY.AddAtom(NA_Atom(  2.466000,  4.961000,  0.001000, NA_Base::NONE,     0, "C7  "));
  THY.AddAtom(NA_Atom( -0.024000,  5.057000,  0.000000, NA_Base::NONE,     1, "C6  "));
  // URA
  bases_.push_back( RefBase('U', "URA", NA_Base::URA) );
  RefBase& URA = bases_.back();
  AddBaseNames("U", URA);
  URA.AddAtom(NA_Atom( -2.481000,  5.354000,  0.000000, NA_Base::NONE,     0, "C1' "));
  URA.AddAtom(NA_Atom( -1.284000,  4.500000,  0.000000, NA_Base::NONE,     1, "N1  "));
  URA.AddAtom(NA_Atom( -1.462000,  3.131000,  0.000000, NA_Base::NONE,     1, "C2  "));
  URA.AddAtom(NA_Atom( -2.563000,  2.608000,  0.000000, NA_Base::ACCEPTOR, 0, "O2  "));
  URA.AddAtom(NA_Atom( -0.302000,  2.397000,  0.000000, NA_Base::DONOR,    1, "N3  "));
  URA.AddAtom(NA_Atom(  0.989000,  2.884000,  0.000000, NA_Base::NONE,     1, "C4  "));
  URA.AddAtom(NA_Atom(  1.935000,  2.094000, -0.001000, NA_Base::ACCEPTOR, 0, "O4  "));
  URA.AddAtom(NA_Atom(  1.089000,  4.311000,  0.000000, NA_Base::NONE,     1, "C5  "));
  URA.AddAtom(NA_Atom( -0.024000,  5.053000,  0.000000, NA_Base::NONE,     1, "C6  "));
  // DEBUG
  //for (BaseArray::const_iterator b = bases_.begin(); b != bases_.end(); ++b)
  //  b->PrintInfo();
}

static inline void CheckHbondValue(int& hbond) {
  if (hbond < 0 || hbond > 2) {
    mprintf("Warning: Invalid value for hbond column (%i): setting to 0 (none).\n", hbond);
    hbond = 0;
  }
}

static inline void CheckRmsValue(int& rms) {
  if (rms < 0 || rms > 1) {
    mprintf("Warning: Invalid value for rms column (%i): setting to 0 (not used for fit).\n",
            rms);
    rms = 0;
  }
}

/** Load a NA base reference from a file. Eventually want to be able to
  * read 3DNA files, but right now the way nastruct works it needs more
  * info than the 3DNA reference PDBs provide, namely each atoms
  * hydrogen bonding status and whether it should be used for rms-fitting.
  * To make it as simple as possible for now, accept the following formats.
  * 1 - Simple (all whitespace-delimited, ignore '#' lines):
  *    Header:
  *      NASTRUCT REFERENCE
  *      <1 char name> <res name 0> [ <res name 1> ...]
  * The 1 char name will be used to ID the underlying base type: A G C T U
  *    Atom lines:
  *      <atom name> <X> <Y> <Z> <hbond> <rms>
  * 2 - Modified PDB:
  *    Atom lines:
  *      ATOM <#> <atom name> <res name> <res num> <X> <Y> <Z> <hbond/occ> <rms/bfac>
  * For hbond, 0 is none, 1 is donor, 2 is acceptor.
  */
int NA_Reference::LoadFromFile(FileName const& fname) {
  CpptrajFile infile;
  if (infile.SetupRead(fname, 0)) return 1;
  RefBase baseIn;
  // Determine if the file is PDB format (like 3DNA references)
  if (PDBfile::ID_PDB(infile)) {
    // PDB format.
    PDBfile pdbIn;
    if (pdbIn.OpenRead(fname)) return 1;
    double XYZ[3];
    float occ, bfac;
    while (pdbIn.NextRecord() != PDBfile::END_OF_FILE) {
      if (pdbIn.RecType() == PDBfile::ATOM) {
        pdbIn.pdb_XYZ( XYZ );
        Atom pAtm = pdbIn.pdb_Atom();
        Residue pRes = pdbIn.pdb_Residue();
        //pdbIn.pdb_OccupancyAndBfactor(occ, bfac);
        // Use ChargeAndRadius, less strict about spacing
        pdbIn.pdb_ChargeAndRadius(occ, bfac);
        int hbond = (int)occ;
        int rms   = (int)bfac;
        CheckHbondValue( hbond );
        CheckRmsValue( rms );
        if (baseIn.empty())
          baseIn = RefBase( pRes.c_str()[0], pRes.Name(), NA_Base::UNKNOWN_BASE );
        // TODO determine base type, hbond, rms fit
        baseIn.AddAtom(NA_Atom(XYZ[0], XYZ[1], XYZ[2], NA_Base::HBType(hbond), rms, pAtm.c_str()));
      } else if (pdbIn.RecType() == PDBfile::END)
        break;
    }
    pdbIn.CloseFile();
  } else {
    if (infile.OpenFile()) return 1;
    std::string line = infile.GetLine();
    if (line.empty() || line.compare(0,18,"NASTRUCT REFERENCE") != 0) {
      mprinterr("Error: Unrecognized reference base format.\n");
      return 1;
    }
    // Next line will contain the 1 char name and all residue names
    const char* ptr = infile.NextLine();
    while (ptr != 0 && ptr[0] == '#') ptr = infile.NextLine();
    ArgList argline( ptr );
    if (argline.Nargs() < 2) {
      mprinterr("Error: Could not read 1 char base name and residue name.\n");
      return 1;
    }
    char baseChar = argline[0][0];
    NA_Base::NAType baseType;
    switch (baseChar) {
      case 'A' : baseType = NA_Base::ADE; break;
      case 'G' : baseType = NA_Base::GUA; break;
      case 'C' : baseType = NA_Base::CYT; break;
      case 'T' : baseType = NA_Base::THY; break;
      case 'U' : baseType = NA_Base::URA; break;
      default  : baseType = NA_Base::UNKNOWN_BASE;
    }
    argline.GetStringNext(); // Essentially "pops" the 1 char name
    line = argline.GetStringNext(); // First residue name
    RefBase::NameArray resNames;
    while (!line.empty()) {
      resNames.push_back( line );
      line = argline.GetStringNext();
    }
    baseIn = RefBase( baseChar, resNames, baseType );
    // Remaining lines contain '<atom name> <X> <Y> <Z> <hbond> <rms>'
    ptr = infile.NextLine();
    while (ptr != 0) {
      argline.SetList( ptr, " " );
      // Skip empty lines and comments
      if (argline.Nargs() > 0 && argline[0][0] != '#') {
        if (argline.Nargs() != 6) {
          mprinterr("Error: Expected 6 columns, got %i\n", argline.Nargs());
          return 1;
        }
        NA_Base::HBType hbtype = NA_Base::NONE;
        if (validInteger(argline[4])) {
          int hbond = convertToInteger(argline[4]);
          CheckHbondValue( hbond );
          hbtype = (NA_Base::HBType)hbond;
        } else {
          switch (argline[4][0]) {
            case 'a' : hbtype = NA_Base::ACCEPTOR; break;
            case 'd' : hbtype = NA_Base::DONOR; break;
            case 'n' : hbtype = NA_Base::DONOR; break;
            default  :
              mprintf("Warning: Unrecognized HB type '%s'. Using NONE.\n", argline[4].c_str());
              hbtype = NA_Base::NONE;
          }
        }
        int rms = convertToInteger(argline[5]);
        CheckRmsValue( rms );
        baseIn.AddAtom(NA_Atom(convertToDouble(argline[1]),
                               convertToDouble(argline[2]),
                               convertToDouble(argline[3]),
                               hbtype, rms, argline[0].c_str()));
      }
      ptr = infile.NextLine();
    }
    infile.CloseFile();
  }
  mprintf("Loaded reference from file '%s' for residues named:", fname.base());
  for (RefBase::name_iterator it = baseIn.nameBegin(); it != baseIn.nameEnd(); ++it)
    mprintf(" %s", *(*it));
  mprintf("\n");
  //baseIn.PrintInfo();
  if (AddBase( baseIn )) return 1;
  return 0;
}

/** Add given reference base. Determine if any of the residue names for the
  * input reference match any existing references. If so, they will be
  * overridden. Easiest way to do this is to insert the given base
  * before any existing references.
  */
int NA_Reference::AddBase(RefBase const& refIn) {
  if (refIn.empty()) {
    mprinterr("Internal Error: Attempting to add an empty reference base.\n");
    return 1;
  }
  BaseArray newBases;
  newBases.reserve( bases_.size() + 1 );
  newBases.push_back( refIn );
  for (BaseArray::const_iterator b = bases_.begin(); b != bases_.end(); ++b) {
    // Check if any refIn names match this base.
    for (RefBase::name_iterator name = refIn.nameBegin(); name != refIn.nameEnd(); ++name)
      if ( b->NameMatches( *name ) )
        mprintf("Warning: New reference residue '%s' will override existing reference.\n",
                *(*name));
    newBases.push_back( *b );
  }
  bases_ = newBases;
  return 0;
}

/** Attempt to find a reference for the given NA residue. */
NA_Reference::RetType
  NA_Reference::SetupBaseRef(NA_Base& baseIn, Topology const& currentParm, int resnum,
                             DataSetList& masterDSL, std::string const& dataname)
{
  // TODO clear baseIn?
  // Determine base from residue name.
  Residue const& RES = currentParm.Res(resnum);
  BaseArray::const_iterator REF = bases_.begin();
  for (; REF != bases_.end(); ++REF)
    if (REF->NameMatches( RES.Name() )) break;
  if (REF == bases_.end()) return NOT_FOUND;

  if (baseIn.Setup_Base(*REF, RES, resnum, currentParm.Atoms(), masterDSL, dataname))
    return BASE_ERROR;
  return BASE_OK;
}

// NA_Reference::AddNameToBaseType()
void NA_Reference::AddNameToBaseType(NameType const& nameIn, NA_Base::NAType typeIn) {
  for (BaseArray::iterator b = bases_.begin(); b != bases_.end(); ++b) {
    if (b->Type() == typeIn) {
      mprintf("\tAdding name '%s' to base '%c'\n", *nameIn, b->BaseChar());
      b->AddName( nameIn );
      break;
    }
  }
}

// ---------- RefBase ----------------------------------------------------------
/** \return true if any of this reference bases names matches given name. */
bool RefBase::NameMatches(NameType const& nameIn) const {
  for (NameArray::const_iterator n = names_.begin(); n != names_.end(); ++n)
    if ( *n == nameIn ) return true;
  return false;
}

/** Print reference base info to STDOUT. */
void RefBase::PrintInfo() const {
  mprintf("Base '%c':", baseChar_);
  for (NameArray::const_iterator n = names_.begin(); n != names_.end(); ++n)
    mprintf(" %s", *(*n));
  mprintf("\n");
  for (NA_Array::const_iterator at = atoms_.begin(); at != atoms_.end(); ++at)
    mprintf("\t%s %6.3f %6.3f %6.3f %i %i\n", at->name(), at->X(), at->Y(), at->Z(),
            (int)at->HB_type(), at->RmsFit());
}

// ---------- NA_Atom ----------------------------------------------------------
/// CONSTRUCTOR - replace any asterisks in name with quote
NA_Atom::NA_Atom(double x, double y, double z, NA_Base::HBType t, int r, const char* n) :
      x_(x), y_(y), z_(z), hb_type_(t), rms_fit_(r), aname_(n)
{
  aname_.ReplaceAsterisk();
}

// ---------- NA_Base ----------------------------------------------------------
NA_Base::NA_Base() :
  pucker_(0),
  rnum_(0),
  c3idx_(-1),
  c5idx_(-1),
  strandNum_(-1),
  bchar_('?'),
  type_(UNKNOWN_BASE)
{
  std::fill( atomIdx_, atomIdx_+6, -1 );
}

// COPY CONSTRUCTOR
NA_Base::NA_Base(const NA_Base& rhs) :
  pucker_(rhs.pucker_),
  rnum_(rhs.rnum_),
  c3idx_(rhs.c3idx_),
  c5idx_(rhs.c5idx_),
  strandNum_(rhs.strandNum_),
  bchar_(rhs.bchar_),
  type_(rhs.type_),
  Ref_(rhs.Ref_),
  anames_(rhs.anames_),
  basename_(rhs.basename_),
# ifdef NASTRUCTDEBUG
  rname_(rhs.rname_),
  refnames_(rhs.refnames_),
# endif
  Inp_(rhs.Inp_),
  hb_(rhs.hb_),
  parmMask_(rhs.parmMask_),
  inpFitMask_(rhs.inpFitMask_),
  refFitMask_(rhs.refFitMask_)
{
  std::copy( rhs.atomIdx_, rhs.atomIdx_+6, atomIdx_ );
}

// ASSIGNMENT
NA_Base& NA_Base::operator=(const NA_Base& rhs) {
  if (this != &rhs) {
    pucker_ = rhs.pucker_;
    rnum_ = rhs.rnum_;
    c3idx_ = rhs.c3idx_;
    c5idx_ = rhs.c5idx_;
    strandNum_ = rhs.strandNum_;
    bchar_ = rhs.bchar_;
    type_ = rhs.type_;
    Ref_ = rhs.Ref_;
    anames_ = rhs.anames_;
    basename_ = rhs.basename_;
#   ifdef NASTRUCTDEBUG
    rname_ = rhs.rname_;
    refnames_ = rhs.refnames_;
#   endif
    Inp_ = rhs.Inp_;
    hb_ = rhs.hb_;
    std::copy( rhs.atomIdx_, rhs.atomIdx_+6, atomIdx_ );
    parmMask_ = rhs.parmMask_;
    inpFitMask_ = rhs.inpFitMask_;
    refFitMask_ = rhs.refFitMask_;
  }
  return *this;
}

/** \return index of given input atom name. */
int NA_Base::FindAtom(NameType const& atname) const {
  int atom = 0;
  for (std::vector<NameType>::const_iterator Name = anames_.begin();
                                             Name != anames_.end(); ++Name)
  {
    if (*Name == atname) return atom;
    ++atom;
  }
  return -1;
}

#ifdef NASTRUCTDEBUG
static const char* HBSTRING[] = {" 0 ", "HBD", "HBA"};
static const char* NAbaseName[] = { "UNK", "ADE", "CYT", "GUA", "THY", "URA" };
#endif
/** Set NA residue reference coordinates for given NA base. Ensure that
  * the atom ordering in the reference matches that in the given parm.
  * If an error occurs the type will be set to UNKNOWN_BASE. 
  */
int NA_Base::Setup_Base(RefBase const& REF, Residue const& RES, int resnum,
                        std::vector<Atom> const& Atoms, 
                        DataSetList& masterDSL, std::string const& dataname) 
{
  type_ = REF.Type();
  int resstart = RES.FirstAtom();
  int resstop  = RES.LastAtom();
  // Create mask for all input coords for this residue
  parmMask_.AddAtomRange(resstart, resstop);
  // Allocate space to hold input coords
  Inp_.SetupFrame( parmMask_.Nselected() );
  hb_.resize( parmMask_.Nselected(), NONE );
  // Save atom names for input coords. Look for specific atom names for
  // calculating things like groove width and pucker.
  int inpatom = 0;
  std::fill( atomIdx_, atomIdx_+6, -1 );
  for (int atom = resstart; atom < resstop; ++atom) {
    anames_.push_back( Atoms[atom].Name() );
    // Is this atom P?
    if (anames_.back() == "P   ")
      atomIdx_[PHOS] = inpatom;
    // Is this atom O4'/O4*?
    else if (anames_.back() == "O4' " || anames_.back() == "O4* ")
      atomIdx_[O4p] = inpatom;
    else if (anames_.back() == "C1' " || anames_.back() == "C1* ")
      atomIdx_[C1p] = inpatom;
    else if (anames_.back() == "C2' " || anames_.back() == "C2* ")
      atomIdx_[C2p] = inpatom;
    else if (anames_.back() == "C3' " || anames_.back() == "C3* ")
      atomIdx_[C3p] = inpatom;
    else if (anames_.back() == "C4' " || anames_.back() == "C4* ")
      atomIdx_[C4p] = inpatom;
    inpatom++;
  }
  // Determine whether sugar atoms are all present.
  bool hasSugarAtoms = (atomIdx_[O4p] != -1 && atomIdx_[C1p] != -1 && 
                        atomIdx_[C2p] != -1 && atomIdx_[C3p] != -1 &&
                        atomIdx_[C4p] != -1);
  // For each atom defined as a reference atom for this base, find the
  // corresponding atom in the parm.
  std::map<int,int> BaseMap;
  int refatom = 0;
  
  for (RefBase::const_iterator ref = REF.begin(); ref != REF.end(); ++ref) {
    inpatom = FindAtom( ref->Name() );
    // Sometimes C1' is listed as C1*; if search for C1' fails look for C1*.
    if (inpatom < 0 && ref->Name() == "C1' ")
      inpatom = FindAtom("C1* ");
    if (inpatom < 0) {
      mprinterr("Error: Ref Atom [%s] not found in NA base [%s].\n",
                ref->name(), RES.c_str());
      return 1;
    } else {
      BaseMap.insert( std::pair<int,int>(inpatom, refatom) );
#     ifdef NASTRUCTDEBUG
      mprintf("Ref atom %i:%s found in parm (%i:%s)\n",refatom+1,ref->name(),
              inpatom+1,*anames_[inpatom]);
#     endif
    }
    ++refatom;
  }
  if (!BaseMap.empty()) {
    // Now create reference frame with same order as parm. Create RMS fit
    // masks for Ref and Inp frames. If atom is indicated as H-bonding store
    // its index in Inp.
    int refidx = 0; // Index in this NA_Base Ref
    for (std::map<int,int>::iterator atom = BaseMap.begin(); 
                                     atom != BaseMap.end(); atom++, refidx++) {
      inpatom = atom->first;  // Index in Inp 
      refatom = atom->second; // Index in NA_RefAtom array 
      // Store type of hbonding atom.
      hb_[inpatom] = REF[refatom].HB_type();
      // Store coords
      Ref_.AddVec3( Vec3(REF[refatom].X(), REF[refatom].Y(), REF[refatom].Z()) );
#     ifdef NASTRUCTDEBUG
      // Store reference atom names
      refnames_.push_back( REF[refatom].Name() );
#     endif
      // Will this atom be used for RMS fitting?
      if (REF[refatom].RmsFit() == 1) {
        inpFitMask_.AddAtom( inpatom );
        refFitMask_.AddAtom( refidx );
#       ifdef NASTRUCTDEBUG
        mprintf("\tFit atom: input parm atom= %i, ref atom= %i\n", inpatom+1, refidx+1);
#       endif
      }
    }
    // Make sure all masks have atoms
    if (parmMask_.None() || inpFitMask_.None() || refFitMask_.None()) {
      mprinterr("Error: One or more masks for NA residue %i has no atoms.\n", resnum+1);
      return 1;
    } else {
      rnum_ = resnum;
      c3idx_ = -1;
      c5idx_ = -1;
      strandNum_ = -1;
      bchar_ = REF.BaseChar();
#     ifdef NASTRUCTDEBUG
      rname_ = RES.Name();
      mprintf("\tSet up residue %i:%s as %s (%c)\n", rnum_+1, *rname_, NAbaseName[type_], bchar_);
      mprintf("\tReference Atoms:\n");
      for (int atom = 0; atom < Ref_.Natom(); ++atom) {
        mprintf("\t\t%s: ", *(refnames_[atom]));
        Ref_.printAtomCoord(atom);
      }
      mprintf("\tResidue is %i atoms:\n", Inp_.Natom());
      for (int atom = 0; atom < (int)anames_.size(); ++atom)
        mprintf("\t\t%s: %i (%s)\n", *(anames_[atom]), atom+1, HBSTRING[hb_[atom]]);
      mprintf("\tP=%i  O4'=%i\n", atomIdx_[PHOS]+1, atomIdx_[O4p]+1);
      parmMask_.PrintMaskAtoms("ParmMask");
      inpFitMask_.PrintMaskAtoms("InputFitMask");
      refFitMask_.PrintMaskAtoms("RefFitMask");
#     endif
    }
  } else {
    mprinterr("Error: Could not set up reference for residue %i\n", resnum+1);
    return 1;
  }
  basename_ = integerToString( rnum_+1 ) + BaseChar();
  // Add any base-related DataSets
  if (hasSugarAtoms) {
    MetaData md(dataname, "pucker", rnum_+1);
    md.SetLegend( basename_ );
    md.SetScalarType( MetaData::PUCKER );
    md.SetScalarMode( MetaData::M_PUCKER );
    // Check if pucker data set is already present
    DataSet* ds = masterDSL.CheckForSet( md );
    if (ds == 0) {
      pucker_ = (DataSet_1D*)masterDSL.AddSet(DataSet::FLOAT, md);
      if (pucker_ == 0) return 1;
    } else {
      // Check that it is the correct type.
      if (ds->Type() != DataSet::FLOAT) {
        mprinterr("Error: Set '%s' already present but is not FLOAT.\n", ds->legend());
        return 1;
      }
      pucker_ = (DataSet_1D*)ds;
    }
  } else
    pucker_ = 0;
  return 0;
}

void NA_Base::SetInputFrame(Frame const& inputFrame) {
  Inp_.SetCoordinates( inputFrame, parmMask_ );
}

void NA_Base::PrintAtomNames() const {
  mprintf("\tInp Atoms:");
  for (std::vector<NameType>::const_iterator aname = anames_.begin();
                                             aname != anames_.end(); ++aname)
    mprintf(" %s", *(*aname));
  mprintf("\n");
}

void NA_Base::CalcPucker(int framenum, PmethodType puckerMethod) {
  if (pucker_ != 0) {
    double pval=0.0, aval, tval;
    switch (puckerMethod) {
      case ALTONA:
        pval = Pucker_AS( C1xyz(), C2xyz(), C3xyz(),
                          C4xyz(), O4xyz(), aval );
      break;
      case CREMER:
        pval = Pucker_CP( C1xyz(), C2xyz(), C3xyz(),
                          C4xyz(), O4xyz(), 0,
                          5, aval, tval );
        break;
    }
    float fval = (float)(pval * Constants::RADDEG);
    pucker_->Add(framenum, &fval);
  }
}

// ---------- NA_Axis ----------------------------------------------------------
// CONSTRUCTOR
NA_Axis::NA_Axis() {}

void NA_Axis::StoreRotMatrix(Matrix_3x3 const& Rin, Vec3 const& vIn) {
  R_ = Rin;
  RX_ = R_.Col1();
  RY_ = R_.Col2();
  RZ_ = R_.Col3();
  origin_ = vIn;
}
 
void NA_Axis::PrintAxisInfo(const char *title) const {
  mprintf("         %s origin: %8.4f %8.4f %8.4f\n",title,origin_[0],origin_[1],origin_[2]);
  mprintf("         %s R_x vec: %8.4f %8.4f %8.4f\n",title,R_[0],R_[3],R_[6]);
  mprintf("         %s R_y vec: %8.4f %8.4f %8.4f\n",title,R_[1],R_[4],R_[7]);
  mprintf("         %s R_z vec: %8.4f %8.4f %8.4f\n",title,R_[2],R_[5],R_[8]);
}

/** Flip the Z and Y axes. Equivalent to rotation around the X axis.
  * Done for antiparallel stranded DNA.
  */
void NA_Axis::FlipYZ() {
  R_[1] = -R_[1]; // -Yx
  R_[4] = -R_[4]; // -Yy
  R_[7] = -R_[7]; // -Yz
  R_[2] = -R_[2]; // -Zx
  R_[5] = -R_[5]; // -Zy
  R_[8] = -R_[8]; // -Zz
  RY_.Neg();
  RZ_.Neg();
}

/** Flip the X and Y axes. Equivalent to rotation around the Z axis.
  * Done for parallel stranded DNA.
  */
void NA_Axis::FlipXY() {
  R_[0] = -R_[0]; // -Xx
  R_[3] = -R_[3]; // -Xy
  R_[6] = -R_[6]; // -Xz
  R_[1] = -R_[1]; // -Yx
  R_[4] = -R_[4]; // -Yy
  R_[7] = -R_[7]; // -Yz
  RX_.Neg();
  RY_.Neg();
}
