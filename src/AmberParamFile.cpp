#include "AmberParamFile.h"
#include "ArgList.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "DataSet_LeapOpts.h"
#include "StringRoutines.h"
#include "Parm/DihedralParmSet.h"
#include "Parm/ParameterSet.h"
#include "Parm/ParmHolder.h"
#include "Constants.h"
#include <cstdio> // sscanf
#include <cctype> // isspace

using namespace Cpptraj::Parm;

const int AmberParamFile::MAXSYMLEN = 16;

/// CONSTRUCTOR
AmberParamFile::AmberParamFile() :
  debug_(0),
  default_scee_(1.2), // AMBER DEFAULT
  default_scnb_(2.0)  // AMBER DEFAULT
{}

/** Set debug level */
void AmberParamFile::SetAmberParamDebug(int d) {
  debug_ = d;
}

/** Set defaults from leap options */
void AmberParamFile::SetDefaults(DataSet_LeapOpts* leapopts) {
  if (leapopts == 0) return;
  default_scee_ = leapopts->SCEE();
  default_scnb_ = leapopts->SCNB();
}

/** Read symbols delimited by - and space. */
int AmberParamFile::read_symbols(const char* ptrIn, std::vector<std::string>& symbols, int nsymbols)
{
  int isymbol = 0;
  bool char_has_been_read = false;
  for (const char* ptr = ptrIn; *ptr != '\0'; ++ptr)
  {
    if (*ptr == '-') {
      isymbol++;
      char_has_been_read = false;
    } else if ((*ptr == ' ' || *ptr == '\t') && isymbol + 1 == nsymbols && char_has_been_read) {
      return (ptr - ptrIn);
    } else {
      symbols[isymbol] += *ptr;
      if (*ptr != ' ') char_has_been_read = true;
    }
  }
  return -1;
}

/// Hold a set of nonbonded parameters
class AmberParamFile::NonbondSet {
  public:
    NonbondSet(std::string const& n) : name_(n) {}

    std::string name_;          ///< Name of set parameters
    ParmHolder<LJparmType> LJ_; ///< Hold LJ 6-12 parameters
};

/// Hold an off-diagonal NB modification
class AmberParamFile::OffdiagNB {
  public:
    OffdiagNB(NameType const& AT1, NameType const& AT2, double sig1, double eps1, double sig2, double eps2) :
      types_(2), LJ1_(sig1, eps1), LJ2_(sig2, eps2)
    {
      types_.AddName(AT1);
      types_.AddName(AT2);
    }

    TypeNameHolder types_;
    LJparmType LJ1_;
    LJparmType LJ2_;
};

/** Read input for atom symbols and masses. */
int AmberParamFile::read_atype(ParameterSet& prm, const char* ptr)
const
{
  // Format (A2,2X,F10.2x,f10.2)
  if (debug_ > 1) mprintf("DEBUG: Atype: %s\n", ptr);
  char kndsym[MAXSYMLEN];
  double amass = 0;
  double atpol = 0;
  // FIXME Atom class needs to be redone to make reading in elements easier.
  int nscan = sscanf(ptr, "%s %lf %lf", kndsym, &amass, &atpol);
  TypeNameHolder types(kndsym);
  Cpptraj::Parm::RetType ret;
  if (nscan > 3)
    mprintf("Warning: CPPTRAJ currently only reads atomic symbol, mass, and polarizability.\n"
            "Warning: Ignoring information after the third column: %s\n", ptr);
  if (nscan >= 3) {
    // Mass and polarization. Check polarization.
    if (atpol < 0.0) {
      bool pol_is_neg_one = (atpol < -0.99 && atpol > -1.01); // NOTE: To be consistent with LEAP behavior
      if (!pol_is_neg_one) {
        mprintf("Warning: %s negative polarization %g - omitting.\n", *(types[0]), atpol);
        return 0;
      }
    } else if (atpol > 15.0) {
      mprintf("Warning: %s polarization is large (%g).\n", *(types[0]), atpol);
    }
    ret = prm.AT().AddParm( types,
                            AtomType(amass, atpol),
                            true );
  } else if (nscan == 2) {
    // Only mass
    ret = prm.AT().AddParm( types,
                            AtomType(amass),
                            true );
  } else {
    mprinterr("Error: Expected atom type, mass, polarizability, got only %i columns.\n", nscan);
    return 1;
  }
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated %s\n", types.TypeNameStr("atom").c_str());
  else if (ret == Cpptraj::Parm::UPDATED) {
    AtomType const& prev = prm.AT().PreviousParm();
    bool diff_mass = (prev.Mass() != amass);
    bool diff_pol  = (prev.Polarizability() != atpol);
    mprintf("Warning: Redefining %s", types.TypeNameStr("atom").c_str());
             //prm.AT().PreviousParm().Mass(), prm.AT().PreviousParm().Polarizability(), amass, atpol);
    if (diff_mass) mprintf(" from mass %g to %g", prev.Mass(), amass);
    if (diff_pol)  mprintf(" from pol %g to %g", prev.Polarizability(), atpol);
    mprintf("\n");
  } else if (ret == Cpptraj::Parm::ERR) {
    mprinterr("Error: Reading %s\n", types.TypeNameStr("atom").c_str());
    return 1;
  }

  //if (ret == Cpptraj::Parm::UPDATED)
  //  mprintf("Warning: Redefining atom type %s\n", kndsym);
  return 0;
}

/** Read input for bond. */
int AmberParamFile::read_bond(ParameterSet& prm, const char* ptr)
const
{
  // Bond parameters
  // IBT , JBT , RK , REQ
  // FORMAT(A2,1X,A2,2F10.2)
  if (debug_ > 1) mprintf("DEBUG: Bond: %s\n", ptr);
  std::vector<std::string> symbols(2);
  int pos = read_symbols(ptr, symbols, 2);
  if (pos < 0) {
    mprinterr("Error: Could not read symbols for bond from %s\n", ptr);
    return 1;
  }
  //mprintf("DEBUG: %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), ptr+pos);
  double RK, REQ;
  int nscan = sscanf(ptr+pos, "%lf %lf", &RK, &REQ);
  if (nscan != 2) {
    mprinterr("Error: Expected RK, REQ, got only %i elements\n", nscan);
    return 1;
  }
  TypeNameHolder types(2);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  Cpptraj::Parm::RetType ret = prm.BP().AddParm(types, BondParmType(RK, REQ), true);

  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated %s\n", types.TypeNameStr("bond").c_str());
  else if (ret == Cpptraj::Parm::UPDATED)
    mprintf("Warning: Redefining %s from RK= %g REQ= %g to RK= %g REQ= %g\n",
            types.TypeNameStr("bond").c_str(),
            prm.BP().PreviousParm().Rk(), prm.BP().PreviousParm().Req(), RK, REQ);
  else if (ret == Cpptraj::Parm::ERR) {
    mprinterr("Error: Reading %s\n", types.TypeNameStr("bond").c_str());
    return 1;
  }

  return 0;
}

/** Read input for angle. */
int AmberParamFile::read_angle(ParameterSet& prm, const char* ptr)
const
{
  // Angle parameters
  // ITT , JTT , KTT , TK , TEQ
  // FORMAT(A2,1X,A2,1X,A2,2F10.2)
  if (debug_ > 1) mprintf("DEBUG: Angle: %s\n", ptr);
  std::vector<std::string> symbols(3);
  int pos = read_symbols(ptr, symbols, 3);
  if (pos < 0) {
    mprinterr("Error: Could not read symbols for angle from %s\n", ptr);
    return 1;
  }
  //mprintf("DEBUG: %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), ptr+pos);
  double TK, TEQ;
  int nscan = sscanf(ptr+pos, "%lf %lf", &TK, &TEQ);
  if (nscan != 2) {
    mprinterr("Error: Expected TK, TEQ, got only %i elements\n", nscan);
    return 1;
  }
  TypeNameHolder types(3);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  types.AddName( symbols[2] );
  Cpptraj::Parm::RetType ret = prm.AP().AddParm(types, AngleParmType(TK, TEQ*Constants::DEGRAD), true);
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated %s\n", types.TypeNameStr("angle").c_str());
  else if (ret == Cpptraj::Parm::UPDATED)
    mprintf("Warning: Redefining %s from TK= %g TEQ= %g to RK= %g REQ= %g\n",
             types.TypeNameStr("angle").c_str(),
             prm.AP().PreviousParm().Tk(), prm.AP().PreviousParm().Teq()*Constants::RADDEG, TK, TEQ);
  else if (ret == Cpptraj::Parm::ERR) {
    mprinterr("Error: Reading %s\n", types.TypeNameStr("angle").c_str());
    return 1;
  }
  //if (ret == Cpptraj::Parm::UPDATED)
  //  mprintf("Warning: Redefining angle type %s - %s - %s\n", *(types[0]), *(types[1]), *(types[2]));

  return 0;
}

/** Read input for dihedral.
  * IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN
  * FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
  * If IPT .eq. 'X ' .and. LPT .eq. 'X ' then any dihedrals in the
  * system involving the atoms "JPT" and and "KPT" are assigned 
  * the same parameters.  This is called the general dihedral type
  * and is of the form "X "-"JPT"-"KPT"-"X ".
  * IDIVF is the factor by which the torsional barrier is divided.
  * Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
  * details. Basically, the actual torsional potential is
  *   (PK/IDIVF) * (1 + cos(PN*phi - PHASE))
  * If PN .lt. 0.0 then the torsional potential is assumed to have more
  * than one term, and the values of the rest of the terms are read from
  * the next cards until a positive PN is encountered.  The negative value
  * of pn is used only for identifying the existence of the next term and 
  * only the absolute value of PN is kept.
  */
int AmberParamFile::read_dihedral(DihedralParmSet& prm, const char* ptr,
                                  std::vector<std::string>& last_symbols,
                                  bool first_char_is_space)
const
{
  // Dihedral parameters
  if (debug_ > 1)
    mprintf("DEBUG: Dihedral: %s\n", ptr);
  std::vector<std::string> symbols(4);
  int pos = 0;
  if (first_char_is_space) {
    // Assume no symbols on this line. Use previous symbols.
    if (last_symbols.empty()) {
      mprinterr("Error: No symbols in dihedral line: %s\n", ptr);
      return 1;
    }
    symbols = last_symbols;
    // Advance past the whitespace
    while (ptr[pos] == ' ') pos++;
  } else {
    pos = read_symbols(ptr, symbols, 4);
    if (pos < 0) {
      mprinterr("Error: Could not read symbols for dihedral from %s\n", ptr);
      return 1;
    }
  }
  //mprintf("DEBUG: %s %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), ptr+pos);
  int IDIVF;
  double PK, PHASE, PN;
  char sSCEE[128];
  char sSCNB[128];
  // TODO note when PN is negative and expect more terms?
  int nscan = sscanf(ptr+pos, "%i %lf %lf %lf %s %s", &IDIVF, &PK, &PHASE, &PN, sSCEE, sSCNB);
  if (nscan < 4) {
    mprinterr("Error: Expected IDIVF, PK, PHASE, PN, got only %i elements\n", nscan);
    return 1;
  }
  if (PN < 0.0) PN = -PN;
  double scee = default_scee_; // AMBER DEFAULT
  double scnb = default_scnb_; // AMBER DEFAULT
  if (nscan == 6) {
    // Check for SCEE/SCNB (GLYCAM)
    if (sSCEE[0] == 'S' && sSCEE[1] == 'C' && sSCEE[2] == 'E' && sSCEE[3] == 'E')
     sscanf( sSCEE, "SCEE=%lf", &scee);
    if (sSCNB[0] == 'S' && sSCNB[1] == 'C' && sSCNB[2] == 'N' && sSCNB[3] == 'B')
     sscanf( sSCNB, "SCNB=%lf", &scnb);
  }
  TypeNameHolder types(4);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  types.AddName( symbols[2] );
  types.AddName( symbols[3] );
  Cpptraj::Parm::RetType ret =
    prm.AddDihParm(types, DihedralParmType(PK / (double)IDIVF, PN, PHASE*Constants::DEGRAD, scee, scnb), true);
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated %s\n", types.TypeNameStr("dihedral").c_str());
  else if (ret == Cpptraj::Parm::UPDATED) // TODO SCEE/SCNB?
    mprintf("Warning: Redefining %s (PN=%g) from PK= %g Phase= %g to PK= %g Phase= %g\n",
             types.TypeNameStr("dihedral").c_str(), PN,
             prm.PreviousParm().Pk(), prm.PreviousParm().Phase()*Constants::RADDEG, PK/(double)IDIVF, PHASE);
  else if (ret == Cpptraj::Parm::ERR) {
    mprinterr("Error: Reading %s\n", types.TypeNameStr("dihedral").c_str());
    return 1;
  }
  //if (ret == Cpptraj::Parm::UPDATED) {
  //  mprintf("Warning: Redefining dihedral type %s - %s - %s - %s (PN=%g)\n",
  //          *(types[0]), *(types[1]), *(types[2]), *(types[3]), PN);
    //mprintf("DEBUG: %s\n", ptr);
  //}
  last_symbols = symbols;
  return 0;
}

/** Read input for improper. */
int AmberParamFile::read_improper(DihedralParmSet& prm, const char* ptr)
const
{
  // Improper parameters
  // IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN
  // FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)
  if (debug_ > 1) mprintf("DEBUG: Improper: %s\n", ptr);
  std::vector<std::string> symbols(4);
  int pos = read_symbols(ptr, symbols, 4);
  if (pos < 0) {
    mprinterr("Error: Could not read symbols for improper from %s\n", ptr);
    return 1;
  }
  //mprintf("DEBUG: %s %s %s %s '%s'\n", symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), ptr+pos);
  double PK, PHASE, PN;
  int nscan = sscanf(ptr+pos, "%lf %lf %lf", &PK, &PHASE, &PN);
  if (nscan != 3) {
    mprinterr("Error: Expected PK, PHASE, PN, got only %i elements\n", nscan);
    return 1;
  }
  if (PN < 0.0) {
    mprintf("Warning: Improper for %s-%s-%s-%s has negative phase (%g)\n",
            symbols[0].c_str(), symbols[1].c_str(), symbols[2].c_str(), symbols[3].c_str(), PN);
    PN = -PN;
  }
  TypeNameHolder types(4);
  types.AddName( symbols[0] );
  types.AddName( symbols[1] );
  types.AddName( symbols[2] );
  types.AddName( symbols[3] );
  //types.SortImproperByAlpha("X"); // FIXME wildcard should be a static var
  Cpptraj::Parm::RetType ret =
    prm.AddDihParm(types, DihedralParmType(PK, PN, PHASE*Constants::DEGRAD), true);
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated %s\n", types.TypeNameStr("improper").c_str());
  else if (ret == Cpptraj::Parm::UPDATED)
    mprintf("Warning: Redefining %s from PK= %g Phase= %g to PK= %g Phase= %g\n",
             types.TypeNameStr("improper").c_str(),
             prm.PreviousParm().Pk(), prm.PreviousParm().Phase()*Constants::RADDEG, PK, PHASE);
  else if (ret == Cpptraj::Parm::ERR) {
    mprinterr("Error: Reading %s\n", types.TypeNameStr("improper").c_str());
    return 1;
  }

  //if (ret == Cpptraj::Parm::UPDATED)
  //  mprintf("Warning: Redefining improper type %s - %s - %s - %s\n",
  //          *(types[0]), *(types[1]), *(types[2]), *(types[3]));

  return 0;
}

/** Read input for LJ 10-12 hbond. */
int AmberParamFile::read_lj1012(ParameterSet& prm, const char* ptr)
const
{
  // Lennard-Jones 10-12 hydrogen bonding term.
  // According to the docs, the format should be:
  // KT1 , KT2 , A , B , ASOLN , BSOLN , HCUT , IC
  // FORMAT(2X,A2,2X,A2,2x,5F10.2,I2)
  // In practive (in e.g. parm91.dat), the actual format appears to be:
  // KT1, KT2, A, B, HCUT
  char KT1[MAXSYMLEN], KT2[MAXSYMLEN];
  double A, B, HCUT;
  //double ASOLN, BSOLN, HCUT;
  //int IC;
  if (debug_ > 1) mprintf("DEBUG: LJ 10-12: %s\n", ptr);
  //int nscan = sscanf(ptr, "%s %s %lf %lf %lf %lf %lf %i\n",
  int nscan = sscanf(ptr, "%s %s %lf %lf %lf", KT1, KT2, &A, &B, &HCUT);
  if (nscan < 5) {
    mprinterr("Error: Expected at least KT1, KT2, A, B, HCUT, got only %i elements\n", nscan);
    return 1;
  }
  TypeNameHolder types(2);
  types.AddName( KT1 );
  types.AddName( KT2 );
  Cpptraj::Parm::RetType ret = prm.HB().AddParm(types, HB_ParmType(A, B, HCUT), true);
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated %s\n", types.TypeNameStr("LJ 10-12 hbond").c_str());
  else if (ret == Cpptraj::Parm::UPDATED) {
    mprintf("Warning: Redefining %s from Asol= %g Bsol= %g Cut= %g to Asol= %g Bsol= %g Cut=%g\n",
            types.TypeNameStr("LJ 10-12 hbond").c_str(),
            prm.HB().PreviousParm().Asol(), prm.HB().PreviousParm().Bsol(), prm.HB().PreviousParm().HBcut(),
            A, B, HCUT);
  } else if (ret == Cpptraj::Parm::ERR) {
    mprinterr("Error: Reading %s\n", types.TypeNameStr("LJ 10-12 hbond").c_str());
    return 1;
  }

  //if (ret == Cpptraj::Parm::UPDATED)
  //  mprintf("Warning: Redefining LJ 10-12 hbond type %s %s\n", *(types[0]), *(types[1]));
  return 0;
}

/** Read IPOL section */
int AmberParamFile::read_ipol(ParameterSet& prm, const char* ptr)
const
{
  //mprintf("DEBUG: Read IPOL.\n");
  int ipolIn = -1;
  if (sscanf(ptr, "%i", &ipolIn) != 1) {
    mprinterr("Error: Could not read IPOL value from line: %s\n", ptr);
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: IPOL value is %i\n", ipolIn);
  return 0;
}

/** Read radius/depth input for LJ 6-12 nonbond. */
int AmberParamFile::read_nb_RE(NonbondSet& nbset, const char* ptr)
const
{
  // ***** ONLY IF KINDNB .EQ. 'RE' *****
  // LTYNB , R , EDEP
  if (debug_ > 1) mprintf("DEBUG: Nonbond: %s\n", ptr);
  //char LTYNB[MAXSYMLEN];
  //double R, EDEP;
  //double R14, E14;
  // This section is a little tricky. Apparently CHARMM-style Amber FF
  // files can have 14 LJ params here. Try to detect this.
  ArgList nbargs( ptr, " " );
  if (nbargs.Nargs() < 3) {
    mprinterr("Error: Expected at least TYPE, R, DEPTH, got %i elements.\n", nbargs.Nargs());
    return 1;
  }
  bool has_14 = false;
  if (nbargs.Nargs() >= 5) {
    if (validDouble( nbargs[3] ) && validDouble( nbargs[4] )) {
      has_14 = true;
    }
  }
  if (has_14) mprintf("DEBUG: NB HAS 1-4 LJ PARAMS.\n"); // FIXME save these
  double R = convertToDouble( nbargs[1] );
  double EDEP = convertToDouble( nbargs[2] );
  TypeNameHolder types( nbargs[0] );
  Cpptraj::Parm::RetType ret = nbset.LJ_.AddParm( types, LJparmType(R, EDEP), true );
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated LJ 6-12 type %s\n", *(types[0]));
  else if (ret == Cpptraj::Parm::UPDATED)
    mprintf("Warning: Redefining LJ 6-12 type %s from R= %g depth= %g to R= %g depth= %g\n", *(types[0]),
            nbset.LJ_.PreviousParm().Radius(), nbset.LJ_.PreviousParm().Depth(), R, EDEP);
  else if (ret == Cpptraj::Parm::ERR)
    mprinterr("Error: Adding LJ 6-12 type %s\n", *(types[0]));
  return 0;
}

/** Read A/C coefficient input for LJ 6-12 nonbond. */
int AmberParamFile::read_nb_AC(NonbondSet& nbset, const char* ptr)
const
{
  // ***** ONLY IF KINDNB .EQ. 'AC' *****
  // LTYNB , A , C
  if (debug_ > 1) mprintf("DEBUG: Nonbond: %s\n", ptr);

  ArgList nbargs( ptr, " " );
  if (nbargs.Nargs() < 3) {
    mprinterr("Error: Expected at least TYPE, A, C, got %i elements.\n", nbargs.Nargs());
    return 1;
  }

  double A = convertToDouble( nbargs[1] );
  double C = convertToDouble( nbargs[2] );
  TypeNameHolder types( nbargs[0] );

  NonbondType nbt(A, C);

  Cpptraj::Parm::RetType ret = nbset.LJ_.AddParm( types, LJparmType(nbt.Radius(), nbt.Depth()), true );
  if (ret == Cpptraj::Parm::SAME)
    mprintf("Warning: Duplicated AC LJ 6-12 type %s\n", *(types[0]));
  else if (ret == Cpptraj::Parm::UPDATED)
    mprintf("Warning: Redefining AC LJ 6-12 type %s from R= %g depth= %g to R= %g depth= %g\n", *(types[0]),
            nbset.LJ_.PreviousParm().Radius(), nbset.LJ_.PreviousParm().Depth(), nbt.Radius(), nbt.Depth());
  else if (ret == Cpptraj::Parm::ERR)
    mprinterr("Error: Adding AC LJ 6-12 type %s\n", *(types[0]));
  return 0;
}

/** Read LJ off-diagonal modifications */
int AmberParamFile::read_ljedit(Oarray& Offdiag, const char* ptr)
const
{
  if (debug_ > 1) mprintf("DEBUG: LJedit: %s\n", ptr);
  // Lennard-Jones sigma and epsilon of the first atom type when it
  // interacts with anything under the normal rules, then the sigma
  // and epsilon of the second atom type when it interacts with the first.
  char AT1[MAXSYMLEN], AT2[MAXSYMLEN];
  double sig1, eps1, sig2, eps2;
  int nscan = sscanf(ptr, "%s %s %lf %lf %lf %lf", AT1, AT2, &sig1, &eps1, &sig2, &eps2);
  if (nscan != 6) {
    //mprinterr("Error: Expected AT1, AT2, SIG1, EPS1, SIG2, EPS2, got %i elements.\n", nscan);
    //return 1;
    // NOTE: Change to a warning to match LEAP behavior
    mprintf("Warning: Expected AT1, AT2, SIG1, EPS1, SIG2, EPS2, got %i elements. Skipping.\n", nscan);
    mprintf("Warning: Incomplete LJEDIT line: %s\n", ptr);
    return 0;
  }
  Offdiag.push_back( OffdiagNB(AT1, AT2, sig1, eps1, sig2, eps2) );
  return 0;
}

/// \return 1 if problem detected with CMAP term
int AmberParamFile::check_cmap(int currentCmapIdx, CmapGridType const& cmap) const {
  if (debug_ > 0) {
    mprintf("DEBUG: Cmap %i '%s' (gridSize= %i) %zu res:",
            currentCmapIdx, cmap.Title().c_str(), cmap.Size(), cmap.ResNames().size());
    for (std::vector<std::string>::const_iterator it = cmap.ResNames().begin();
                                                it != cmap.ResNames().end(); ++it)
      mprintf(" %s", it->c_str());
    mprintf("\n");
  }
  if (!cmap.CmapIsValid()) {
    mprinterr("Error: CMAP term %i '%s' is not valid.\n", currentCmapIdx, cmap.Title().c_str());
    // Why not?
    if (cmap.Resolution() < 1) mprinterr("Error: Resolution is 0.\n");
    if ((int)(cmap.Resolution() * cmap.Resolution()) != cmap.Size()) mprinterr("Error: Resolution %u*%u != grid size %i\n",cmap.Resolution(), cmap.Resolution(), cmap.Size());
    if (cmap.AtomNames().size() != 5) mprinterr("Error: Expected 5 atom names (got %zu)\n", cmap.AtomNames().size());
    if (cmap.NcmapResNames() < 1) mprinterr("Error: No residue names set up.\n");
    return 1;
  }
  if (!cmap.CmapNresIsValid()) {
    mprintf("Warning: The number of expected residue names for CMAP %i '%s' (%i) is not equal to the actual number %zu\n",
            currentCmapIdx,
            cmap.Title().c_str(),
            cmap.NcmapResNames(),
            cmap.ResNames().size());
  }
  return 0;
}

/// Add default atom names to CMAP term
static inline void add_cmap_default_atoms(CmapGridType& cmap) {
  cmap.AddAtomName("C");
  cmap.AddAtomName("N");
  cmap.AddAtomName("CA");
  cmap.AddAtomName("C");
  cmap.AddAtomName("N");
  cmap.AddResOffset( -1 );
  cmap.AddResOffset(  0 );
  cmap.AddResOffset(  0 );
  cmap.AddResOffset(  0 );
  cmap.AddResOffset(  1 );
}

/** Read CMAP section
  * Default CMAP atoms are Amber phi/psi:
  * Res  : 0   1   1    1   2
  * Atom : C - N - CA - C - N
  */ 
int AmberParamFile::read_cmap(CmapGridType& currentCmap, ParameterSet& prm, CmapType& currentCmapFlag,
                              std::string const& line, int& cmap_count_is_index)
const
{
  if (line.empty()) return 0;
  // At this point, line should have no leading whitespace

  if (line[0] == '%') {
    ArgList argline(line);
    if (argline.Nargs() > 1) {
      // NOTE: CMAP_COUNT is used differently in different files. For example,
      //       in ff19SB it is used as an index of the current CMAP parameter
      //       being read. However, in ff12polL it is used to give a total
      //       count of CMAP parameters
      if (argline[0] == "%FLAG") {
        if (argline[1] == "CMAP_COUNT") {
          // New CMAP term. Ignore the index for now. If a previous CMAP
          // was read make sure its OK.
          //if (!currentCmap.empty()) {
          //  if (currentCmap.AtomNames().empty()) add_cmap_default_atoms( currentCmap );
          //  if (check_cmap(prm.CMAP().size()+1, currentCmap)) return 1;
          //  prm.CMAP().AddParm( currentCmap, true, debug_ );
          //  currentCmap = CmapGridType();
          //}
          int cmapcount = argline.getKeyInt("CMAP_COUNT", -1);
          if (debug_ > 0)
            mprintf("DEBUG: Cmap count: %i\n", cmapcount);
          //if ( cmapcount != (int)prm.CMAP().size()+1 )
          //  mprintf("Warning: CMAP term is not in numerical order. CMAP_COUNT %i, expected %zu\n",
          //          cmapcount, prm.CMAP().size());
          if (cmap_count_is_index == -1) {
            if  (cmapcount == (int)prm.CMAP().size()) {
              if (debug_ > 0)
                mprintf("DEBUG: Assuming CMAP count is based on index.\n");
              cmap_count_is_index = 1;
            } else {
              if (debug_ > 0)
                mprintf("DEBUG: Assuming CMAP count is the total number of CMAPS\n");
              cmap_count_is_index = 0;
            }
          } else if (cmap_count_is_index == 1 && cmapcount != (int)prm.CMAP().size()) {
            mprintf("Warning: CMAP term is not in numerical order. CMAP_COUNT %i, expected %zu\n",
                    cmapcount, prm.CMAP().size());
  
          }
          return 0;
        } else if (argline[1] == "CMAP_TITLE") {
          currentCmapFlag = CMAP_TITLE;
          // New CMAP term. If a previous CMAP was read make sure it is OK.
          if (!currentCmap.empty()) {
            if (currentCmap.AtomNames().empty()) add_cmap_default_atoms( currentCmap );
            if (check_cmap(prm.CMAP().size()+1, currentCmap)) return 1;
            prm.CMAP().AddParm( currentCmap, true, debug_ );
            currentCmap = CmapGridType();
          }

          return 0;
        } else if (argline[1] == "CMAP_RESLIST") {
          currentCmapFlag = CMAP_RESLIST;
          int nres = argline.getKeyInt("CMAP_RESLIST", -1);
          if (nres < 1) {
            mprinterr("Error: Bad CMAP # residues: %s\n", line.c_str());
            return 1;
          }
          if (debug_ > 0)
            mprintf("DEBUG: CMAP_RESLIST expect %i residues.\n", nres);
          currentCmap.SetNumCmapRes( nres );
          return 0;
        } else if (argline[1] == "CMAP_RESOLUTION") {
          int cmapres = argline.getKeyInt("CMAP_RESOLUTION", -1);
          if (debug_ > 0)
            mprintf("DEBUG: Cmap resolution: %i\n", cmapres);
          if (cmapres < 1) {
            mprinterr("Error: Bad CMAP resolution: %s\n", line.c_str());
            return 1;
          }
          currentCmap.SetResolution( cmapres );
          return 0;
        } else if (argline[1] == "CMAP_PARAMETER") {
          currentCmapFlag = CMAP_PARAMETER;
          return 0;
        } else if (argline[1] == "CMAP_ATMLIST") {
          currentCmapFlag = CMAP_ATMLIST;
          // Remove any dashes as well
          ArgList atmArgs(line, " -");
          if (atmArgs.Nargs() < 7) {
            mprinterr("Error: Malformed CMAP_ATMLIST line; expected 5 elements, got %i\n", atmArgs.Nargs() - 2);
            mprinterr("Error: %s\n", line.c_str());
            return 1;
          }
          if (debug_ > 0)
            mprintf("DEBUG: Atom names: %s %s %s %s %s\n",
                    atmArgs[2].c_str(),
                    atmArgs[3].c_str(),
                    atmArgs[4].c_str(),
                    atmArgs[5].c_str(),
                    atmArgs[6].c_str());
          currentCmap.AddAtomName( atmArgs[2] );
          currentCmap.AddAtomName( atmArgs[3] );
          currentCmap.AddAtomName( atmArgs[4] );
          currentCmap.AddAtomName( atmArgs[5] );
          currentCmap.AddAtomName( atmArgs[6] );
        } else if (argline[1] == "CMAP_RESIDX") {
          currentCmapFlag = CMAP_RESIDX;
          if (argline.Nargs() < 7) {
            mprinterr("Error: Malformed CMAP_RESIDX line; expected 5 elements, got %i\n", argline.Nargs() - 2);
            mprinterr("Error: %s\n", line.c_str());
            return 1;
          }
          int residxs[5];
          for (int i = 0; i < 5; i++)
            residxs[i] = argline.getNextInteger(0);
          if (debug_ > 0)
            mprintf("DEBUG: Residue offsets: %i %i %i %i %i\n",
                    residxs[0], residxs[1], residxs[2], residxs[3], residxs[4]);
          for (int i = 0; i < 5; i++)
            currentCmap.AddResOffset( residxs[i] );
        } else {
          mprintf("Warning: Unhandled CMAP FLAG %s\n", argline[1].c_str());
        }
      } else if (argline[0][0] == '%' && argline[0][1] == 'C' && argline[0][2] == 'O') {
        // Assume comment
        currentCmapFlag = CMAP_COMMENT;
        if (debug_ > 0) mprintf("DEBUG: CMAP comment: %s\n", line.c_str());
        return 0;
      }
    } // END line args > 1
  } // END line starts with %

  if (currentCmapFlag == CMAP_PARAMETER) {
    // Another annoyance; CMAP parameters can be in regular formatted 
    // 8 column F9.5, or comma-separated list
    size_t found = line.find_first_of(",");
    if (found != std::string::npos) {
      // Comma-separated
      ArgList termArgs(line, " ,\t\r");
      for (int i = 0; i != termArgs.Nargs(); i++) {
        double term = termArgs.getNextDouble(0.0);
        if (debug_ > 1) mprintf("DEBUG: cmap term %f\n", term);
        currentCmap.AddToGrid( term );
      }
    } else {
      // Assume 8F9.5
      double terms[8];
      int nterms = sscanf(line.c_str(), "%9lf%9lf%9lf%9lf%9lf%9lf%9lf%9lf",
                          terms, terms+1, terms+2, terms+3,
                          terms+4, terms+5, terms+6, terms+7);
      if (debug_ > 1) {
        for (int i = 0; i != nterms; i++)
          mprintf("DEBUG: cmap term %f\n", terms[i] );
      }
      for (int i = 0; i != nterms; i++) {
        currentCmap.AddToGrid( terms[i] );
      }
    }
  } else if (currentCmapFlag == CMAP_RESLIST) {
    ArgList resnames( line );
    std::string rn = resnames.GetStringNext();
    while (!rn.empty()) {
      currentCmap.AddResName( rn );
      rn = resnames.GetStringNext();
    }
  } else if (currentCmapFlag == CMAP_TITLE) {
    if (debug_ > 0)
      mprintf("DEBUG: CMAP TITLE: %s\n", line.c_str());
    currentCmap.SetTitle( line );
  }
  return 0;
}

/** Assign nonbond parameters from a NonbondSet to ParameterSet */
int AmberParamFile::assign_nb(ParameterSet& prm, NonbondSet const& nbset) const {
  // Nonbonds
  for (ParmHolder<LJparmType>::const_iterator it = nbset.LJ_.begin();
                                              it != nbset.LJ_.end(); ++it)
  {
    ParmHolder<AtomType>::iterator at = prm.AT().GetParam( it->first );
    if (at == prm.AT().end()) {
      mprinterr("Error: Nonbond parameters defined for previously undefined type '%s'.\n",
                *(it->first[0]));
      return 1;
    }
    //mprintf("DEBUG: searched type %s hasLJ=%i\n", *(at->first[0]), (int)at->second.HasLJ());
    if (at->second.HasLJ() && at->second.LJ() != it->second)
      mprintf("Warning: Changing type %s radius from %g to %g, depth from %g to %g\n", *(it->first[0]),
              at->second.LJ().Radius(), it->second.Radius(),
              at->second.LJ().Depth(), it->second.Depth());
    at->second.SetLJ( it->second );
  }
  return 0;
}

/** Assign off-diagonal nonbond parameters from an OffdiagNB array to ParameterSet */
int AmberParamFile::assign_offdiag(ParameterSet& prm, Oarray const& Offdiag) const {
  // Do off diagonal NB mods
  if (!Offdiag.empty()) {
    for (Oarray::const_iterator od = Offdiag.begin(); od != Offdiag.end(); ++od)
    {
      if (debug_ > 1) mprintf("DEBUG: Off diag %s %s\n", *(od->types_[0]), *(od->types_[1]));
      // Set type 1
      ParmHolder<AtomType>::iterator it = prm.AT().GetParam( od->types_[0] );
      if (it == prm.AT().end()) {
        mprinterr("Error: Off-diagonal nonbond parameters defined for previously undefined type '%s'.\n",
                  *(od->types_[0]));
        return 1;
      }
      it->second.SetLJ( od->LJ1_ );
      // Set off-diagonal for type1-type2
      // FIXME different combine rules?
      prm.NB().AddParm(od->types_, od->LJ1_.Combine_LB(od->LJ2_), false);
    }
  } // END off-diagonal NB mods
  return 0;
}

/** Read parameters from Amber frcmod file. */
int AmberParamFile::ReadFrcmod(ParameterSet& prm, FileName const& fname) const
{
  int cmap_count_is_index = -1;
  // Set wildcard character for dihedrals and impropers
  prm.DP().SetWildcard('X');
  prm.IP().SetWildcard('X');
  // Read title
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open file '%s' as Amber FF.\n", fname.full());
    return 1;
  }
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read anything from Amber FF file %s\n", fname.full());
    return 1;
  }
  std::vector<std::string> last_symbols;
  last_symbols.reserve(4);
  std::string title(ptr);
  mprintf("\tTitle: %s\n", title.c_str());
  prm.SetParamSetName( title );
  prm.SetParamSetFile( fname );
  NonbondSet nbset(title);
  Oarray Offdiag;
  DihedralParmSet dihPrm( debug_ );
  DihedralParmSet impPrm( debug_ );
  // Read file
  SectionType section = UNKNOWN;
  CmapType currentCmapFlag = CMAP_INITIAL;
  CmapGridType currentCmap;
  ptr = infile.Line();
  while (ptr != 0) {
    bool first_char_is_space = (*ptr == ' ');
    // Advance to first non-space char
    while (isspace(*ptr) && *ptr != '\0') ++ptr;
    // Is this a recognized section keyword?
    if (*ptr != '\0') {
      std::string line(ptr);
      if      (line.compare(0, 4, "MASS") == 0) section = ATYPE;
      else if (line.compare(0, 4, "BOND") == 0) section = BOND;
      else if (line.compare(0, 4, "ANGL") == 0) section = ANGLE;
      else if (line.compare(0, 4, "DIHE") == 0) section = DIHEDRAL;
      else if (line.compare(0, 4, "IMPR") == 0) section = IMPROPER;
      else if (line.compare(0, 4, "HBON") == 0) section = LJ1012;
      else if (line.compare(0, 4, "CMAP") == 0) section = CMAP;
      else if (line.compare(0, 4, "IPOL") == 0) section = IPOL;
      else if (line.compare(0, 6, "LJEDIT") == 0) {
        section = LJEDIT;
        prm.SetHasLJparams( true );
      } else if (line.compare(0, 4, "NONB") == 0) {
        section = NONBOND;
        prm.SetHasLJparams( true );
        // TODO check RE
      } else {
        //mprintf("DEBUG: Section %i: %s\n", (int)section, ptr);
        int err = 0;
        if (section == ATYPE)
          err = read_atype(prm, ptr);
        else if (section == BOND)
          err = read_bond(prm, ptr);
        else if (section == ANGLE)
          err = read_angle(prm, ptr);
        else if (section == DIHEDRAL)
          err = read_dihedral(dihPrm, ptr, last_symbols, first_char_is_space);
        else if (section == IMPROPER)
          err = read_improper(impPrm, ptr);
        else if (section == LJ1012)
          err = read_lj1012(prm, ptr);
        else if (section == NONBOND)
          err = read_nb_RE(nbset, ptr);
        else if (section == LJEDIT)
          err = read_ljedit(Offdiag, ptr);
        else if (section == CMAP)
          err = read_cmap(currentCmap, prm, currentCmapFlag, line, cmap_count_is_index);
        else if (section == IPOL)
          err = read_ipol(prm, ptr);
        if (err != 0) {
          mprinterr("Error: Reading line: %s\n", ptr);
          return 1;
        }
      }
    }
    ptr = infile.Line();
  }
  // Check last cmap
  if (!currentCmap.empty()) {
    if (currentCmap.AtomNames().empty()) add_cmap_default_atoms( currentCmap );
    if (check_cmap(prm.CMAP().size(), currentCmap)) return 1;
    prm.CMAP().AddParm( currentCmap, true, debug_ );
  }
  // Add dihedrals
  if (dihPrm.ToDihParm(prm.DP())) return 1;
  if (impPrm.ToImpParm(prm.IP())) return 1;
  // Nonbonds
  if (assign_nb(prm, nbset)) return 1;
  // Off-diagonal NB modifications
  if (assign_offdiag(prm, Offdiag)) return 1;

  if (debug_ > 0) prm.Debug();
  infile.CloseFile();

  return 0;
}

/** Read parameters from Amber main FF parameter file. */
int AmberParamFile::ReadParams(ParameterSet& prm, FileName const& fname,
                               std::string const& nbsetnameIn) const
{
  // Set wildcard character for dihedrals and impropers
  prm.DP().SetWildcard('X');
  prm.IP().SetWildcard('X');
  // For files with > 1 set of NB params
  typedef std::vector<NonbondSet> NbSetArrayType;
  NbSetArrayType NBsets;
  // For holding equivalent NB type names
  typedef std::vector<NameType> Narray;
  typedef std::vector<Narray> XNarray;
  XNarray EquivalentNames;
  // For holding off-diagonal mods
  Oarray Offdiag;

  // Read title
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open file '%s' as Amber FF.\n", fname.full());
    return 1;
  }
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Could not read anything from Amber FF file %s\n", fname.full());
    return 1;
  }
  std::string title(ptr);
  mprintf("\tTitle: %s\n", title.c_str());
  prm.SetParamSetName( title );
  prm.SetParamSetFile( fname );
  std::vector<std::string> last_symbols;
  last_symbols.reserve(4);
  DihedralParmSet dihPrm( debug_ );
  DihedralParmSet impPrm( debug_ );
  // Read file
  int readNbType = 0;
  SectionType section = ATYPE;
  ptr = infile.Line();
  while (ptr != 0) {
    bool first_char_is_space = (*ptr == ' ');
    // Advance to first non-space char
    while (*ptr == ' ' && *ptr != '\0') ++ptr;
    //mprintf("DEBUG: First char: %c (%i)\n", *ptr, (int)*ptr);
    int read_err = 0;
    if (*ptr == '\0') {
      // Section Change
      if (section != UNKNOWN) {
        if (section == NONBOND) {
          // Do a lookahead to see if there are multiple NB sets.
          // It will either be another set, END, or LJEDIT
          ptr = infile.Line();
          while (*ptr == ' ' && *ptr != '\0') ++ptr;
          //mprintf("DEBUG: NONBOND First char: %c (%i)\n", *ptr, (int)*ptr);
          std::string nbline(ptr);
          if (nbline == "END") {
            if (debug_ > 0) mprintf("END\n");
            section = UNKNOWN;
          } else if (nbline == "LJEDIT") {
            section = LJEDIT;
            prm.SetHasLJparams( true );
          } // Otherwise assume another nonbond section
            //else mprintf("DEBUG: Assuming another nonbond section.\n");
        } else {
          if (debug_ > 0) mprintf("SECTION %i change to %i\n", (int)section, (int)section + 1);
          section = (SectionType)((int)section + 1);
        }
      }
      // Special cases
      if (section == HYDROPHILIC) {
        // Special case: hydrophilic atom types. Could be multiple lines.
        // FORMAT(20(A2,2X))
        bool read_hydrophilic = true;
        while (read_hydrophilic) {
          ptr = infile.Line();
          if (debug_ > 1) mprintf("DEBUG: Hydrophilic: %s\n", ptr);
          // Take advantage of the fact that we expect whitespace-delimiters
          ArgList hsymbols( ptr, " " );
          for (int iarg = 0; iarg != hsymbols.Nargs(); ++iarg)
            prm.AddHydrophilicAtomType( NameType( hsymbols[iarg] ) );
          read_hydrophilic = (hsymbols.Nargs() > 19);
        }
        if (debug_ > 0) mprintf("DEBUG: Read %u hydrophilic atom types.\n", prm.NhydrophilicAtomTypes());
        section = (SectionType)((int)section + 1);
        // Look ahead. Older parm files have the hydrophilic section delimited
        // with a newline.
        ptr = infile.Line();
        while (*ptr == ' ' && *ptr != '\0') ++ptr;
        if (*ptr != '\0') continue;
      } else if (section == NONBOND) {
        prm.SetHasLJparams( true );
        if (debug_ > 1) mprintf("DEBUG: Begin nonbond section.\n");
        // Special case: first read the line.
        // LABEL , KINDNB
        // FORMAT(A4,6X,A2)
        if (NBsets.empty())
          ptr = infile.Line();
        ArgList nb_args(ptr);
        //nb_args.PrintDebug();
        if (nb_args.Nargs() != 2) {
          mprintf("Warning: Line %i: Expected 2 elements, nonbond label and kind, got %i elements: %s\n", infile.LineNumber(), nb_args.Nargs(), ptr);
          section = UNKNOWN;
        } else {
          if (debug_ > 0) mprintf("DEBUG: NB label= %s  NB kind = %s\n", nb_args[0].c_str(), nb_args[1].c_str());
          std::string const& nb_kind = nb_args[1];
          if (nb_kind[0] == 'R' && nb_kind[1] == 'E')
            readNbType = 0;
          else if (nb_kind[0] == 'A' && nb_kind[1] == 'C')
            readNbType = 1;
          else {
            mprinterr("Error: Nonbond parameters are not of type 'RE' (Rmin, Epsilon) or 'AC' (A, C coefficients).\n");
            return 1;
          }
          NBsets.push_back( NonbondSet( nb_args[0] ) );
        }
      }
    } else if (section == ATYPE) {
      read_err = read_atype(prm, ptr);
    } else if (section == BOND) {
      read_err = read_bond(prm, ptr);
    } else if (section == ANGLE) {
      read_err = read_angle(prm, ptr);
    } else if (section == DIHEDRAL) {
      read_err = read_dihedral(dihPrm, ptr, last_symbols, first_char_is_space);
    } else if (section == IMPROPER) {
      read_err = read_improper(impPrm, ptr);
    } else if (section == LJ1012) {
      read_err = read_lj1012(prm, ptr);
    } else if (section == NB_EQUIV) {
      // EQUIVALENCING ATOM SYMBOLS FOR THE NON-BONDED 6-12 POTENTIAL PARAMETERS
      // IORG , IEQV(I) , I = 1 , 19
      // FORMAT(20(A2,2X))
      if (debug_ > 1) mprintf("DEBUG: Nonbond equiv: %s\n", ptr);
      EquivalentNames.push_back( Narray() );
      ArgList equiv_line( ptr, " " );
      for (int iarg = 0; iarg != equiv_line.Nargs(); iarg++)
        EquivalentNames.back().push_back( equiv_line[iarg] );
    } else if (section == NONBOND) {
      /// Check for a badly formatted END
      if (ptr[0] != '\0' && ptr[1] != '\0' && ptr[2] != '\0' &&
          ptr[0] == 'E'  && ptr[1] == 'N'  && ptr[2] == 'D')
      {
        mprintf("Warning: 'END' encountered before blank line in NONBOND section.\n");
        section = UNKNOWN;
      } else if (readNbType == 0)
        // ***** ONLY IF KINDNB .EQ. 'RE' *****
        read_err = read_nb_RE(NBsets.back(), ptr);
      else if (readNbType == 1)
        read_err = read_nb_AC(NBsets.back(), ptr);
    } else if (section == LJEDIT) {
      read_err = read_ljedit(Offdiag, ptr);
    }
    if (read_err != 0) {
      mprinterr("Error: Reading line: %s\n", ptr);
      return 1;
    } 
    ptr = infile.Line();
  } // END loop over file.

  // Add dihedrals
  if (dihPrm.ToDihParm(prm.DP())) return 1;
  if (impPrm.ToImpParm(prm.IP())) return 1;
  // Deal with nonbond and equivalence
  if (!NBsets.empty()) {
    int nbsetidx = 0;
    if (NBsets.size() > 1 || !nbsetnameIn.empty()) {
      // Choose from multiple nbsets
      if (nbsetnameIn.empty()) {
        mprintf("Warning: Parm set in file '%s' contains %zu nonbonded parameter sets\n"
                "Warning:  but no specific set has been specified. First set will be used.\n", fname.full(), NBsets.size());
        mprintf("Warning: Nonbonded sets:");
        for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it)
          mprintf(" %s", it->name_.c_str());
        mprintf("\n");
      } else {
        nbsetidx = -1;
        for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it) {
          if (it->name_ == nbsetnameIn) {
            nbsetidx = it - NBsets.begin();
            break;
          }
        }
        if (nbsetidx < 0) {
          mprinterr("Error: Nonbonded set '%s' not found.\n", nbsetnameIn.c_str());
          mprinterr("Error: Need to specify one of");
          for (NbSetArrayType::const_iterator it = NBsets.begin(); it != NBsets.end(); ++it)
            mprinterr(" %s", it->name_.c_str());
          mprinterr("\n");
          return 1;
        }
      }
    }
    mprintf("\tUsing nonbonded parm set: %s\n", NBsets[nbsetidx].name_.c_str());
    prm.SetNbParamName( NBsets[nbsetidx].name_ );
    if (assign_nb(prm, NBsets[nbsetidx])) return 1;
    
    // Do equivalent atoms.
    for (XNarray::const_iterator equivAts = EquivalentNames.begin();
                                 equivAts != EquivalentNames.end(); ++equivAts)
    {
      // First name is the type to copy TODO check size?
      Narray::const_iterator typeName = equivAts->begin();
      ParmHolder<AtomType>::const_iterator at0 = prm.AT().GetParam( *typeName );
      if (at0 == prm.AT().end()) {
        mprinterr("Error: Equivalent atom type '%s' not found.\n", *(*typeName) );
        return 1;
      }
      ++typeName;
      for (; typeName != equivAts->end(); ++typeName) {
        ParmHolder<AtomType>::iterator at1 = prm.AT().GetParam( *typeName );
        if (at1 == prm.AT().end()) {
          mprintf("Warning: Equivalent atom type '%s' (base type '%s') not found.\n", *(*typeName), *(at0->first[0]));
          //return 1;
        } else {
          if (debug_ > 1) mprintf("DEBUG: Equiv '%s' => '%s'\n", *(at0->first[0]), *(*typeName));
          at1->second.SetLJ( at0->second.LJ() );
        }
      }
    } // END loop over EquivalentNames
  } // END nonbond parameters

  // Do off diagonal NB mods
  if (assign_offdiag(prm, Offdiag)) return 1;

  if (debug_ > 0) prm.Debug();
  infile.CloseFile();
  return 0;
}

// DataIO_AmberFF::WriteData()
int AmberParamFile::WriteParams(ParameterSet& prm, FileName const& fname) const
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) return 1;
  // 1 - Title
  std::string title;
  if (prm.ParamSetName().empty())
    title.assign("Parameters written from CPPTRAJ");
  else
    title = prm.ParamSetName();
  outfile.Printf("%s\n", title.c_str());
  // 2 - Atom symbols and masses
  for (ParmHolder<AtomType>::const_iterator at = prm.AT().begin(); at != prm.AT().end(); ++at)
  {
    std::string asym = at->first[0].Truncated();
    if (asym.size() > 2)
      mprintf("Warning: Atom symbol %s is larger than 2 characters, which breaks Amber FF format.\n", asym.c_str());
    outfile.Printf("%-2s", asym.c_str());
    if (at->second.Mass() < 10.0)
      outfile.Printf(" %-10.3f", at->second.Mass());
    else if (at->second.Mass() < 100.0)
      outfile.Printf(" %-10.2f", at->second.Mass());
    else
      outfile.Printf(" %-10.1f", at->second.Mass());
    outfile.Printf(" %10.3f\n", at->second.Polarizability());

    //outfile.Printf("%-2s  %-10.2f %-10.2f\n", asym.c_str(), at->second.Mass(), at->second.Polarizability());
  }

  return 0;
}
