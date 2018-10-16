#include "CharmmParamFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // RemoveTrailingWhitespace()
#include "ArgList.h"
#include "Constants.h"

static inline std::string Input(const char* line) {
  std::string input;
  const char* ptr = line;
  while (*ptr != '\0') {
    if (*ptr == '!') break;
    input += *ptr;
    ++ptr;
  }
  return input;
}

/** Read a single line of CHARMM input from file, respecting any line
  * continuation and ignoring any comments. All letters are converted to
  * upper case.
  * \return 1 if there is more input, 0 if EOF.
  */
int CharmmParamFile::ReadInput( std::string& input, BufferedLine& infile ) const {
  input.clear();
  const char* line = infile.Line();
  if (line == 0) {
    return 0; 
  }
  while (line != 0) {
    // Strip away leading whitespace.
    while (*line != '\0' && isspace(*line)) ++line;
    // If title or comment, return now.
    if (line[0] == '*') {
      input.assign(line);
      return 1;
    } else if (line[0] == '!') {
      input.assign(line);
      return 1; 
    }
    // Read up to any comment character.
    const char* ptr = line;
    while (*ptr != '\0') {
      if (*ptr == '!') break;
      input += toupper(*ptr);
      ++ptr;
    }
    // Remove trailing whitespace
    RemoveTrailingWhitespace( input );
    // See if we need more input
    if (input.size() > 0 && input[input.size()-1] == '-')
      line = infile.Line();
    else
      line = 0;
  }
  return 1;
}

static inline bool ChmCmd(std::string const& cmd, const char* key) {
  return (cmd.compare(0, 4, key)==0);
}

/** Read CHARMM parameters from specified file into given parameter set. */
int CharmmParamFile::ReadParams(ParameterSet& prm, FileName const& nameIn, int debugIn) const {
  BufferedLine infile;

  mprintf("\tReading CHARMM parameters from '%s'\n", nameIn.full());

  if (infile.OpenFileRead( nameIn )) {
    mprinterr("Error: Could not open file '%s'\n", nameIn.full());
    return 1;
  }

  std::string input;

  enum ModeType { NONE = 0, PARAM, TOP };
  ModeType mode = NONE;
  enum SectionType { UNKNOWN, ATOMS, BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED, RES, IGNORE };
  SectionType currentSection = UNKNOWN;

  std::string currentResName;
  double currentResQ = 0.0;

  while (ReadInput(input, infile)) {
    if (input.empty()) continue;
    if (input[0] == '*') {
      if (debugIn > 0) mprintf("DEBUG: Title: %s\n", input.c_str());
    } else if (input[0] == '!') {
      if (debugIn > 0) mprintf("DEBUG: Comment: %s\n", input.c_str());
    } else {
      ArgList args( input, " \t" );
      if (args.Nargs() > 0) {
        const char* line = input.c_str(); // DEBUG
        // Line not blank.
        //while (args[args.Nargs()-1] == "-") {
        //  args.MarkArg(args.Nargs()-1);
        //  line = infile.Line();
        //  args.Append( ArgList(Input(line)) );
        //}
        if (debugIn > 1)
          mprintf("DBG: %s\n", args.ArgLine());
        // Process Command
        // ----- Common to TOP and PARAM --------- 
        if (ChmCmd(args[0],"END")) {
          // END read FIXME is it ok for this to end everything?
          currentSection = UNKNOWN;
          mode = NONE;
        } else if (ChmCmd(args[0],"READ")) {
          currentSection = UNKNOWN;
          if (args.hasKey("RTF"))
            mode = TOP;
          else if (args.hasKey("PARAM")) // FIXME really only need para
            mode = PARAM;
          else
            mode = NONE;
        } else if (ChmCmd(args[0],"MASS") && args.Nargs() == 4) {
          args.MarkArg(0);
          args.MarkArg(1);
          args.MarkArg(2);
          prm.AT().AddAtomType(args[2], AtomType(args.getNextDouble(0)));
        } else if (ChmCmd(args[0], "ATOM")) {
          if (mode != TOP) {
            currentSection = ATOMS; //mode = PARAM;
            if (debugIn > 0)
              mprintf("DEBUG: Section ATOMS, line %i\n", infile.LineNumber());
          }
        } else if (ChmCmd(args[0], "BOND")) {
          if (mode != TOP) {
            currentSection = BONDS; //mode = PARAM;
            if (debugIn > 0)
              mprintf("DEBUG: Section BONDS, line %i\n", infile.LineNumber());
          }
        } else if (ChmCmd(args[0],"ANGL") || ChmCmd(args[0],"THET")) {
          if (mode != TOP) {
            currentSection = ANGLES; //mode = PARAM;
            if (debugIn > 0)
              mprintf("DEBUG: Section ANGLES, line %i\n", infile.LineNumber());
          }
        } else if (ChmCmd(args[0],"DIHE") || ChmCmd(args[0],"PHI")) {
          if (mode != TOP) {
            currentSection = DIHEDRALS; //mode = PARAM;
            if (debugIn > 0)
              mprintf("DEBUG: Section DIHEDRALS, line %i\n", infile.LineNumber());
          }
        } else if (ChmCmd(args[0],"IMPR") || ChmCmd(args[0],"IMPH")) {
          if (mode != TOP) {
            currentSection = IMPROPERS; //mode = PARAM;
            if (debugIn > 0)
              mprintf("DEBUG: Section IMPROPERS, line %i\n", infile.LineNumber());
          }
        } else if (ChmCmd(args[0], "CMAP")) {
          currentSection = IGNORE; //mode = NONE;
          mprintf("Warning: Ignoring CMAP section.\n");
        // ----- Below are PARAM only ------------
        } else if (ChmCmd(args[0], "NONB")) {
          currentSection = NONBONDED; //mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section NONBONDED, line %i\n", infile.LineNumber());
        } else if (ChmCmd(args[0],"HBON")) {
          currentSection = IGNORE; //mode = NONE;
            mprintf("Warning: Ignoring HBOND section.\n");
        } else if (ChmCmd(args[0],"NBFI")) {
          currentSection = IGNORE;
          mprintf("Warning: Ignoring NBFIX section.\n");
        // ----- Below are TOP only --------------
        } else if (ChmCmd(args[0], "RESI") || ChmCmd(args[0], "PRES")) {
          if (args.Nargs() < 3)
            mprintf("Warning: Malformed residue command: %s\n", args.ArgLine());
          args.MarkArg(0);
          currentResName = args.GetStringNext();
          currentResQ = args.getNextDouble(0.0);
          mode = TOP;
          mprintf("Residue: %s %g\n", currentResName.c_str(), currentResQ);
        } else if (ChmCmd(args[0],"DECL")) {
          mprintf("Warning: Skipping DECL %s\n", args[1].c_str());
          mode = TOP;
        } else {
          // Parm section-specific commands
          if (currentSection == BONDS) {
            // BOND PARAMETERS
            if (args.Nargs() < 4)
              mprintf("Warning: Bad syntax for bond parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              AtomTypeHolder types(2);
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              prm.AT().CheckForAtomType( types[0] );
              prm.AT().CheckForAtomType( types[1] );
              double rk = args.getNextDouble(0);
              double req = args.getNextDouble(0);
              prm.BP().AddParm(types, BondParmType(rk, req), false);
            }
          } else if (currentSection == ANGLES) {
            // ANGLE PARAMETERS
            if (args.Nargs() < 5)
              mprintf("Warning: Bad syntax for angle parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              AtomTypeHolder types(3);
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              prm.AT().CheckForAtomType( types[0] );
              prm.AT().CheckForAtomType( types[1] );
              prm.AT().CheckForAtomType( types[2] );
              double tk = args.getNextDouble(0);
              double teq = args.getNextDouble(0);
              prm.AP().AddParm(types, AngleParmType(tk, teq*Constants::DEGRAD), false);
              if (args.Nargs() > 5) {
                // UREY-BRADLEY
                AtomTypeHolder utypes(2);
                utypes.AddName(types[0]);
                utypes.AddName(types[2]);
                tk = args.getNextDouble(0);
                teq = args.getNextDouble(0);
                prm.UB().AddParm(utypes, BondParmType(tk, teq), false);
              }
            }
          } else if (currentSection == DIHEDRALS || currentSection == IMPROPERS) {
            // DIHEDRAL/IMPROPER PARAMETERS
            if (args.Nargs() < 7)
              mprintf("Warning: Bad syntax for dihedral parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              AtomTypeHolder types(4, "X"); // X is wildcard character
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              prm.AT().CheckForAtomType( types[0] );
              prm.AT().CheckForAtomType( types[1] );
              prm.AT().CheckForAtomType( types[2] );
              prm.AT().CheckForAtomType( types[3] );
              double pk = args.getNextDouble(0);
              double pn = args.getNextDouble(0);
              double phase = args.getNextDouble(0) * Constants::DEGRAD;
              if (currentSection == DIHEDRALS)
                prm.DP().AddParm(types, DihedralParmType(pk, pn, phase, 1.0, 1.0), false);
              else
                prm.IP().AddParm(types, DihedralParmType(pk, pn, phase), false);
            }
          } else if (currentSection == NONBONDED) {
            // NONBONDED PARAMETERS TODO do not add if not already present
            if (args.Nargs() < 4)
              mprintf("Warning: Bad syntax for nonbond parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              prm.SetHasLJparams( true );
              // TODO handle 1-4 stuff
              NameType at = args.GetStringNext();
              double epsilon = args.getNextDouble(0.0); // skip
              epsilon = args.getNextDouble(0.0); // negative by convention
              double radius = args.getNextDouble(0.0);
              // Determine if there are wildcard chars.
              bool hasWC = false;
              for (char* atp = at.C_Array(); *atp != '\0'; ++atp) {
                if (*atp == '*') hasWC = true;
                // Replace the CHARMM single wildcard with Amber/UNIX
                if (*atp == '%') {
                  *atp = '?';
                  hasWC = true;
                }
              }
              if (hasWC) {
                // Check against all current atom types.
                for (AtomTypeArray::const_iterator it = prm.AT().begin();
                                                   it != prm.AT().end(); ++it)
                {
                  if (it->first.Match( at )) {
                    int idx = it->second;
                    prm.AT().UpdateType(idx).SetRadius( radius );
                    prm.AT().UpdateType(idx).SetDepth( -epsilon );
                  }
                }
              } else {
                // Single type
                int idx = prm.AT().AtomTypeIndex( at );
                if (idx == -1) {
                  mprintf("Warning: Nonbond parameters defined for type '%s' without MASS card."
                          " Skipping.\n", *at);
                } else {
                  prm.AT().UpdateType(idx).SetRadius( radius );
                  prm.AT().UpdateType(idx).SetDepth( -epsilon );
                }
              }
            }
          } 
        } // END input determination 
      }  // END line not blank
    } // END if not title or comment 
  } // END loop over file read
  if (debugIn > 0) 
    prm.Debug();

  return 0;
}

/** Write CHARMM parameters from specified set into given file. */
int CharmmParamFile::WriteParams(ParameterSet& prm, FileName const& nameIn, int debugIn) const {
  CpptrajFile outfile;
  if (outfile.OpenWrite(nameIn)) return 1;
  // Title
  outfile.Printf("* CHARMM parameters stream file generated by cpptraj.\n"
                 "*\n");
  outfile.Printf("\nread param card flex append\n"
                 "* Parameters written from cpptraj.\n"
                 "*\n");
  // ATOMS
  outfile.Printf("\nATOMS\n");
  for (AtomTypeArray::const_iterator it = prm.AT().begin();
                                     it != prm.AT().end(); ++it)
    outfile.Printf("%-4s %3i  %-8s%10.5f\n", "MASS", -1, *(it->first), prm.AT()[it->second].Mass());

  // BONDS
  if (!prm.BP().empty()) {
    outfile.Printf("\nBONDS\n");
    for (ParmHolder<BondParmType>::const_iterator it = prm.BP().begin();
                                                  it != prm.BP().end(); ++it)
      outfile.Printf("%-8s %-8s %8.3f %10.4f\n",
                     *(it->first[0]), *(it->first[1]),
                     it->second.Rk(), it->second.Req());
  }

  // ANGLES
  if (!prm.AP().empty()) {
    outfile.Printf("\nANGLES\n");
    for (ParmHolder<AngleParmType>::const_iterator it = prm.AP().begin();
                                                   it != prm.AP().end(); ++it)
      outfile.Printf("%-8s %-8s %-8s %8.3f %10.4f\n",
                     *(it->first[0]), *(it->first[1]), *(it->first[2]),
                     it->second.Tk(), it->second.Teq());
  }

  // DIHEDRALS
  if (!prm.DP().empty()) {
    outfile.Printf("\nDIHEDRALS\n");
    for (ParmHolder<DihedralParmType>::const_iterator it = prm.DP().begin();
                                                      it != prm.DP().end(); ++it)
      outfile.Printf("%-8s %-8s %-8s %-8s %10.4f %2i %8.2f\n",
                     *(it->first[0]), *(it->first[1]), *(it->first[2]), *(it->first[3]),
                     it->second.Pk(), (int)it->second.Pn(),
                     it->second.Phase()*Constants::RADDEG);
  }

  // END
  outfile.Printf("\nEND\n");
  outfile.CloseFile();
  return 0;
}
