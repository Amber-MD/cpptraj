#include "CharmmParamFile.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
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

/** Read CHARMM parameters from specified file into given parameter set. */
int CharmmParamFile::ReadParams(ParameterSet& prm, FileName const& nameIn, int debugIn) {
  BufferedLine infile;

  mprintf("\tReading CHARMM parameters from '%s'\n", nameIn.full());

  if (infile.OpenFileRead( nameIn )) {
    mprinterr("Error: Could not open file '%s'\n", nameIn.full());
    return 1;
  }

  const char* line = infile.Line();
  if (line == 0) return 1;

  enum ModeType { NONE = 0, PARAM, TOP };
  ModeType mode = NONE;
  enum SectionType { UNKNOWN, ATOMS, BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED, IGNORE };
  SectionType currentSection = UNKNOWN;
  while (line != 0) {
    // Strip away leading whitespace.
    while (*line != '\0' && isspace(*line)) ++line;
    if (line[0] == '*') {
      if (debugIn > 0)
        mprintf("DEBUG: Title: %s\n", line);
    } else if (line[0] == '!') {
      if (debugIn > 0)
        mprintf("DEBUG: Comment: %s\n", line);
    } else {
      // Input() will read everything up to comment character
      ArgList args( Input(line), " \t" );
      if (args.Nargs() > 0) {
        // Line not blank. Check for continuation
        while (args[args.Nargs()-1] == "-") {
          args.MarkArg(args.Nargs()-1);
          line = infile.Line();
          args.Append( ArgList(Input(line)) );
        }
        if (debugIn > 1)
          mprintf("DBG: %s\n", args.ArgLine());
        // Determine input
        if (args.hasKey("END")) {
          // END read
          currentSection = UNKNOWN;
          mode = NONE;
        } else if (args.Nargs() >= 2 && args[0] == "read" && args.hasKey("param")) {
          // 'read param' command
          mode = PARAM; 
        } else if (args.hasKey("ATOMS"))     {
          currentSection = ATOMS; mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section ATOMS, line %i\n", infile.LineNumber());
        } else if (args.hasKey("BONDS"))     {
          currentSection = BONDS; mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section BONDS, line %i\n", infile.LineNumber());
        } else if (args.hasKey("ANGLES"))    {
          currentSection = ANGLES; mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section ANGLES, line %i\n", infile.LineNumber());
        } else if (args.hasKey("DIHEDRALS")) {
          currentSection = DIHEDRALS; mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section DIHEDRALS, line %i\n", infile.LineNumber());
        } else if (args.hasKey("IMPROPER") || args.hasKey("IMPROPERS")) {
          currentSection = IMPROPERS; mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section IMPROPERS, line %i\n", infile.LineNumber());
        } else if (args.hasKey("NONBONDED")) {
          currentSection = NONBONDED; mode = PARAM;
          if (debugIn > 0)
            mprintf("DEBUG: Section NONBONDED, line %i\n", infile.LineNumber());
        } else if (args.hasKey("HBOND")) {
          currentSection = IGNORE; mode = NONE;
          mprintf("Warning: Ignoring HBOND section.\n");
        } else if (args.hasKey("CMAP")) {
          currentSection = IGNORE; mode = NONE;
          mprintf("Warning: Ignoring CMAP section.\n");
        } else if (mode == PARAM) {
          // ----- Reading parameters ------------
          if (currentSection == ATOMS) {
            // ATOM TYPES (masses really)
            if (args.hasKey("MASS") && args.Nargs() == 4) {
              args.MarkArg(0);
              args.MarkArg(1);
              args.MarkArg(2);
              prm.AT().AddAtomType(args[2], AtomType(args.getNextDouble(0)));
            }
          } else if (currentSection == BONDS) {
            if (args.Nargs() < 4)
              mprintf("Warning: Bad syntax for bond parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              // BOND PARAMETERS
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
            if (args.Nargs() < 5)
              mprintf("Warning: Bad syntax for angle parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              // ANGLE PARAMETERS
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
            if (args.Nargs() < 7)
              mprintf("Warning: Bad syntax for dihedral parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              // DIHEDRAL/IMPROPER PARAMETERS
              AtomTypeHolder types(4, "X"); // X is wildcard character
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              types.AddName( args.GetStringNext() );
              prm.AT().CheckForAtomType( types[0] );
              prm.AT().CheckForAtomType( types[1] );
              prm.AT().CheckForAtomType( types[2] );
              double pk = args.getNextDouble(0);
              double pn = args.getNextDouble(0);
              double phase = args.getNextDouble(0) * Constants::DEGRAD;
              if (currentSection == DIHEDRALS)
                prm.DP().AddParm(types, DihedralParmType(pk, pn, phase, 1.0, 1.0), false);
              else
                prm.IP().AddParm(types, DihedralParmType(pk, pn, phase), false);
            }
          } else if (currentSection == NONBONDED) {
            if (args.Nargs() < 4)
              mprintf("Warning: Bad syntax for nonbond parameter on line %i: %s\n", infile.LineNumber(), line);
            else {
              prm.SetHasLJparams( true );
              // NONBONDED PARAMETERS TODO do not add if not already present
              // TODO handle 1-4 stuff
              NameType at = args.GetStringNext();
              int idx = prm.AT().AtomTypeIndex( at );
              if (idx == -1) {
                mprinterr("Error: Nonbond parameters defined for type '%s' without MASS card.\n",
                          *at);
                return 1;
              }
              double epsilon = args.getNextDouble(0.0); // skip
              epsilon = args.getNextDouble(0.0); // negative by convention
              double radius = args.getNextDouble(0.0);
              prm.AT().UpdateType(idx).SetRadius( radius );
              prm.AT().UpdateType(idx).SetDepth( -epsilon );
            }
          }
          // -------------------------------------
        } // END input determination 
      }
    }
    line = infile.Line();
  }
  if (debugIn > 0) 
    prm.Debug();

  return 0;
}

/** Read CHARMM parameters from specified file into given parameter set. */
int CharmmParamFile::WriteParams(ParameterSet& prm, FileName const& nameIn, int debugIn) {
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
