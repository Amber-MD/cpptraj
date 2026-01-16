#include "GB_Params.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

static const char* GB_RadiiTypeStr_[] = {
  "Bondi radii", // 0
  "Amber 6 modified Bondi radii", // 1
  "Modified Bondi radii", // 2
  "Radii optimized for Amber charges by Huo and Kollman", // 3
  "H(N)-modified Bondi radii", // 6
  "PARSE radii", // 7
  "ArgH and AspGluO modified Bondi2 radii", // 8
  "Unknown GB radii set"
};

static const char* GB_RadiiTypeKey_[] = {
  "bondi",   // 0
  "amber6",  // 1
  "mbondi",  // 2
  "pbamber", // 3
  "mbondi2", // 6
  "parse",   // 7
  "mbondi3", // 8
  0
};

static const char* GB_RadiiAmberFlag_[] = {
  "Bondi radii (bondi)",                              // 0
  "amber6 modified Bondi radii (amber6)",             // 1
  "modified Bondi radii (mbondi)",                    // 2
  "Huo and Kollman optimized radii (pbamber)",        // 3
  "H(N)-modified Bondi radii (mbondi2)",              // 6
  "Parse radii (parse)",                              // 7
  "ArgH and AspGluO modified Bondi2 radii (mbondi3)", // 8
  "Unknown radius set (leap needs to be modified!)"
};

/** \return key corresponding to type. */
const char* Cpptraj::Parm::GbTypeKey(GB_RadiiType t) {
  return GB_RadiiTypeKey_[t];
}

/** \return char string corresponding to type. */
std::string Cpptraj::Parm::GbTypeStr(GB_RadiiType t) {
  return std::string(GB_RadiiTypeStr_[t]);
}

/** \return Amber topology FLAG */
std::string Cpptraj::Parm::GbAmberFlag(GB_RadiiType t) {
  return std::string(GB_RadiiAmberFlag_[t]);
}

/** \return GB_RadiiType corresponding to string. */
Cpptraj::Parm::GB_RadiiType Cpptraj::Parm::GbTypeFromKey(std::string const& key) {
  if (!key.empty()) {
    for (int i = 0; i < (int)UNKNOWN_GB; i++) {
      if ( key == std::string(GB_RadiiTypeKey_[i]) ) return (GB_RadiiType)i;
    }
  }
  return UNKNOWN_GB;
}

/// These are the numeric iGBparm values used in LEaP
static const int GB_RadiiTypeIGB_[] = {
  0, // bondi
  1, // amber6
  2, // mbondi
  3, // pbamber
  6, // mbondi2
  7, // parse
  8, // mbondi3
  999999
};

// =============================================================================
/** CONSTRUCTOR */
Cpptraj::Parm::GB_Params::GB_Params() : gbradii_(Cpptraj::Parm::UNKNOWN_GB) {}

/** Print keywords to stdout. */
std::string Cpptraj::Parm::GB_Params::HelpText() {
  std::string out("gb {");
  for (int ig = 0; ig != (int)UNKNOWN_GB; ig++) {
    if (ig > 0) out.append("|");
    out.append(std::string( GbTypeKey((GB_RadiiType)ig) ));
  }
  out.append("}");
  return out;
}

/** Init radii set */
int Cpptraj::Parm::GB_Params::Init_GB_Radii(ArgList& argIn, Cpptraj::Parm::GB_RadiiType gbRadIn)
{
  if (gbRadIn == Cpptraj::Parm::UNKNOWN_GB) {
    // No radii set specified. Check for keyword.
    gbradii_ = Cpptraj::Parm::MBONDI; // Default
    std::string gbset = argIn.GetStringKey("gb");
    if (!gbset.empty()) {
      gbradii_ = Cpptraj::Parm::GbTypeFromKey( gbset );
      if (gbradii_ == Cpptraj::Parm::UNKNOWN_GB) {
        mprinterr("Error: Unknown GB radii set: %s\n", gbset.c_str());
        return 1;
      }
    }
  } else {
    // Use passed-in GB radii set
    gbradii_ = gbRadIn;
  }
  mprintf("\tGB radii set: %s\n", Cpptraj::Parm::GbTypeStr(gbradii_).c_str());
  return 0;
}

/** Assign GB radii. */
void Cpptraj::Parm::GB_Params::assign_gb_radii(Topology& top, Cpptraj::Parm::GB_RadiiType radiiSet)
{
  int iGBparm = GB_RadiiTypeIGB_[radiiSet];
  for (int iat = 0; iat != top.Natom(); iat++)
  {
    Atom const& currentAtom = top[iat];
    Residue const& currentRes = top.Res(currentAtom.ResNum());
    double dGBrad = 0.0;
//    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
//    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
//    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar,
//                &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
//                &iElement, &iHybridization, sDesc);
    if (iGBparm < 3 || iGBparm == 6 || iGBparm == 8) {
      // Bondi or modified Bondi radii
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBrad = 1.2;

          // Make the modifications that hydrogen radii
          // depend upon the atoms they are bonded to.  
          // iGBparm=1 corresponds to Amber 6, JACS 122:2489 (2000);
          // iGBparm=2 adds the update of Biopolymers 56: 275 (2001)   
          if (currentAtom.Nbonds() > 0) {

	    // For multiply bonded Hydrogen atoms use the first bond for
            // determining modified GB radii.  WAT contains multiply bonded
	    // Hydrogen atoms so do not emit a warning.
            Atom const& aAtomA = top[currentAtom.Bond(0)];
            if (iGBparm == 1 || iGBparm == 2) {
              switch (aAtomA.Element()) {
                case Atom::CARBON: dGBrad = 1.3; break;      // Carbon
                case Atom::OXYGEN: dGBrad = 0.8; break;      // Oxygen
                case Atom::SULFUR: dGBrad = 0.8; break;      // Sulfur
                case Atom::NITROGEN:
                  if (iGBparm == 2) {
                    dGBrad = 1.3;
                  }
                  break;      // Nitrogen, mbondi
                case Atom::HYDROGEN:        // Special case: water hydrogen
                  if (aAtomA.Type() == "HW" || aAtomA.Type() == "hw") {
                                dGBrad = 0.8;
		  }
                  break;
                default : break;
              }
            }
	    else if (iGBparm == 6 || iGBparm == 8) {

              // Try Alexey's scheme
              if (aAtomA.Element() == Atom::NITROGEN) {
                dGBrad = 1.3;
                if (iGBparm == 8) {

                  // Update residue as appropriate
                  //if (saPAtom->iResidueIndex != iResidueIndex) {
                  //  iResidueIndex = saPAtom->iResidueIndex;
                  //  cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                  //                iResidueIndex - 1)->sName;
                  //  if (strlen(cPTemp) > 3) {
                  //    cPTemp += (strlen(cPTemp) - 3);
		  //  }
                  //}

		  // Adjust Arg HH and HE
                  if (currentRes.Name() == "ARG" && (currentAtom.Name().Match("HH*") || currentAtom.Name() == "HE"))
                  {
                  //if (!strcmp(cPTemp, "ARG") &&
		  //    !(strncmp(sAtomName(saPAtom->aAtom), "HH", 2) &&
                  //      strcmp(sAtomName(saPAtom->aAtom), "HE"))) {
                    dGBrad = 1.17;
                  }
                }
              }
            }
          }
          else {
            mprintf("Warning: Unbonded Hydrogen atom %s in %s.\n"
                    " Cannot determine the requested GB radius for this atom.\n"
                    " Writing the unmodified Bondi GB radius.\n",
                    *(currentAtom.Name()),
                    top.TruncResNameNum(currentAtom.ResNum()).c_str());
          }
          break;
        case Atom::CARBON:

	  // Use the mass of the carbon atom. We are testing for
          // carbons here. C1 == CH, C2 == CH2, C3 == CH3. UA carbons
          // have a larger radius (2.2), so we want to make sure that
          // the C1, C2, and C3 atom types _really_ correspond to UA
          // UA carbons. C1 atoms should have a mass of 12.01 + 1.01,
          // C2 should be 12.01 + 2.02, and C3 should be 12.01 + 3.03.
          // This mneumonic will not work for 13C named "C1". This is
          // a (hopefully) temporary kludge. FIXME not sure this first line is right
          //if (strncmp(sType, "C1", 2) && strncmp(sType, "C2", 2) && strncmp(sType, "C3", 2)) {
	  if (currentAtom.Type() != "C1" && currentAtom.Type() != "C2" && currentAtom.Type() != "C3") {
            dGBrad = 1.7;
	  }
          else if (currentAtom.Type() == "C1" && currentAtom.Mass() < 13.0) {
            dGBrad = 1.7;
	  }
	  else if (currentAtom.Type() == "C2" && currentAtom.Mass() < 14.0) {
	    dGBrad = 1.7;
	  }
	  else if (currentAtom.Type() == "C3" && currentAtom.Mass() < 15.0) {
	    dGBrad = 1.7;
	  }
          else {
            dGBrad = 2.2;
	  }
          break;
        case Atom::NITROGEN:
          dGBrad = 1.55;
          break;
        case Atom::OXYGEN:
          dGBrad = 1.5;
          if (iGBparm == 8) {

            // Update residue as appropriate
            //if (saPAtom->iResidueIndex != iResidueIndex) {
            //  iResidueIndex = saPAtom->iResidueIndex;
            //  cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, iResidueIndex - 1)->sName;
            //  if (strlen(cPTemp) > 3) {
            //    cPTemp += (strlen(cPTemp) - 3);
	    //  }
            //}
	    
            // Adjust Asp OD and Glu OE, and terminal OXT
            if ( (currentRes.Name() == "ASP" && currentAtom.Name().Match("OD*")) ||
                 (currentRes.Name() == "AS4" && currentAtom.Name().Match("OD*")) ||
                 (currentRes.Name() == "GLU" && currentAtom.Name().Match("OE*")) ||
                 (currentRes.Name() == "GL4" && currentAtom.Name().Match("OE*")) ||
                 (currentAtom.Name() == "OXT") ||
                 (iat+1 < top.Natom() && top[iat+1].Name() == "OXT") // FIXME O next to OXT in topology seems very hacky
               )
            {
            //if (!(strcmp(cPTemp, "ASP") || strncmp(sAtomName(saPAtom->aAtom), "OD", 2)) ||
            //      !(strcmp(cPTemp, "AS4") || strncmp(sAtomName(saPAtom->aAtom), "OD", 2)) ||
            //      !(strcmp(cPTemp, "GLU") || strncmp(sAtomName(saPAtom->aAtom), "OE", 2)) ||
            //      !(strcmp(cPTemp, "GL4") || strncmp(sAtomName(saPAtom->aAtom), "OE", 2)) ||
            //      (!strcmp(sAtomName(saPAtom->aAtom), "OXT") ||
            //      (i + 1 < iVarArrayElementCount(uUnit->vaAtoms) &&
            //       !strcmp(sAtomName(PVAI(uUnit->vaAtoms, SAVEATOMt, i + 1)->aAtom),
            //      "OXT")))) {
	      dGBrad = 1.4;
	    }
          }
          break;
        case Atom::FLUORINE:
          dGBrad = 1.5;
          break;
        case Atom::SILICON:
          dGBrad = 2.1;
          break;
        case Atom::PHOSPHORUS:
          dGBrad = 1.85;
          break;
        case Atom::SULFUR:
          dGBrad = 1.8;
          break;
        case Atom::CHLORINE:
          dGBrad = 1.7;
          break;
        default:
          dGBrad = 1.5;
          break;
      }
    }
    else if (iGBparm == 3) {

      // Radii from Huo & Kollman
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBrad = 1.15;
          break;
        case Atom::CARBON:
          dGBrad = 1.85;
          break;
        case Atom::NITROGEN:
          dGBrad = 1.64;
          break;
        case Atom::OXYGEN:
          dGBrad = 1.53;
          break;
        case Atom::FLUORINE:
          dGBrad = 1.53;
          break;
        case Atom::PHOSPHORUS:
          dGBrad = 2.02;
          break;
        case Atom::SULFUR:
          dGBrad = 2.00;
          break;
        case Atom::CHLORINE:
          dGBrad = 1.97;
          break;
        case Atom::BROMINE:
          dGBrad = 2.03;
          break;
        case Atom::IODINE:
          dGBrad = 2.10;
          break;
        default:
          dGBrad = 1.5;
          break;
      }
    }
    else if (iGBparm == 7) {

      // Parse radii
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBrad = 1.00;
          break;
        case Atom::CARBON:
          dGBrad = 1.70;
          break;
        case Atom::NITROGEN:
          dGBrad = 1.50;
          break;
        case Atom::OXYGEN:
          dGBrad = 1.40;
          break;
        case Atom::SULFUR:
          dGBrad = 1.85;
          break;
        default:
          dGBrad = 1.50;
          break;
          // Radii from J. Phys. Chem. 1994, 98, 1978-1988
      }
    }
    top.SetAtom(iat).SetGBradius( dGBrad );
  } // END loop over atoms
}

/** Assign GB screening parameters. */
void Cpptraj::Parm::GB_Params::assign_gb_screen(Topology& top, Cpptraj::Parm::GB_RadiiType radiiSet)
{
  int iGBparm = GB_RadiiTypeIGB_[radiiSet];
  for (int iat = 0; iat != top.Natom(); iat++)
  {
    Atom const& currentAtom = top[iat];
    double dGBscreen = 0;
    if (iGBparm < 4 || iGBparm == 6 || iGBparm == 8) {

      // For now, hardwire the Bondi radii
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBscreen = 0.85;
          break;
        case Atom::CARBON:
          dGBscreen = 0.72;
          break;
        case Atom::NITROGEN:
          dGBscreen = 0.79;
          break;
        case Atom::OXYGEN:
          dGBscreen = 0.85;
          break;
        case Atom::FLUORINE:
          dGBscreen = 0.88;
          break;
        case Atom::PHOSPHORUS:
          dGBscreen = 0.86;
          break;
        case Atom::SULFUR:
          dGBscreen = 0.96;
          break;
        default:
          dGBscreen = 0.8;
          break;          // or should fail??
      }
    }
    else if (iGBparm == 4) {    // param for Jayaram et al. 'GB'
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBscreen = 0.8461;
          break;
        case Atom::CARBON:
          dGBscreen = 0.9615;
          break;
        case Atom::NITROGEN:
          dGBscreen = 0.9343;
          break;
        case Atom::OXYGEN:
          dGBscreen = 1.0088;
          break;
        case Atom::SODIUM:
          dGBscreen = 1.0000;
          break;
        case Atom::MAGNESIUM:
          dGBscreen = 1.0000;
          break;          // Set by HG
        case Atom::PHOSPHORUS:
          dGBscreen = 1.0700;
          break;
        case Atom::SULFUR:
          dGBscreen = 1.1733;
          break;
        default:
          dGBscreen = 0.8000;
          break;          // Set by HG
      }
    }
    else if (iGBparm == 5) {

      // Param for Jayaram et al. 'MGB'
      switch (currentAtom.Element()) {
        case Atom::HYDROGEN:
          dGBscreen = 0.8846;
          break;
        case Atom::CARBON:
          dGBscreen = 0.9186;
          break;
        case Atom::NITROGEN:
          dGBscreen = 0.8733;
          break;
        case Atom::OXYGEN:
          dGBscreen = 0.8836;
          break;
        case Atom::SODIUM:
          dGBscreen = 1.0000;
          break;
        case Atom::MAGNESIUM:
          dGBscreen = 1.0000;
          break;          // Set by HG
        case Atom::PHOSPHORUS:
          dGBscreen = 0.9604;
          break;
        case Atom::SULFUR:
          dGBscreen = 0.9323;
          break;
        default:
          dGBscreen = 0.8000;
          break;          // Set by HG
      }
    }
    top.SetAtom(iat).SetGBscreen( dGBscreen );
  } // END loop over atoms
}

/** Assign GB radii and screening parameters based on the given radius set. */
int Cpptraj::Parm::GB_Params::Assign_GB_Radii(Topology& top)
const
{
  if (gbradii_ == UNKNOWN_GB) {
    mprinterr("Error: Unknown GB radii set.\n");
    return 1;
  }
  mprintf("\tUsing GB radii set: %s\n", GB_RadiiTypeStr_[gbradii_]);
  top.SetGBradiiSet( std::string(GB_RadiiAmberFlag_[gbradii_]) );
  assign_gb_radii( top, gbradii_ );
  assign_gb_screen( top, gbradii_ );

  return 0;
}
