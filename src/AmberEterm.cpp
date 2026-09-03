#include "AmberEterm.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NoTrailingWhitespace, validDouble

using namespace Cpptraj;

AmberEterm::AmberEterm() { }

/** Names corresponding to FieldType. */
const char* AmberEterm::Enames_[] = {
  "Etot",   "EPtot",  "GMAX",  "BOND",
  "ANGLE",  "DIHED",  "VDW",   "EELEC",      "EGB",     "EPB", "ECAVITY", "EDISPER",
  "VDW1-4", "EEL1-4", "RST",   "EAMBER",     "Density",
  "RMS",    "EKtot",  "ESURF", "EAMD_BOOST", "VOLUME",  "TEMP",
  "PRESS",  "DVDL",   "CMAP",  0
};

/** \return FieldType corresponding to given term name, or N_FIELDTYPES if
  *         not recognized.
  */
AmberEterm::FieldType AmberEterm::getTermIdx(std::string const& name) const {
  NameIdxMap::const_iterator it = termIdxMap_.find( name );
  if (it == termIdxMap_.end()) {
    return (FieldType)N_FIELDTYPES;
  } else {
    return (FieldType)it->second;
  }
}


/** Parse the given line for energy terms of format <name>=<value>. */
int AmberEterm::GetAmberEterms(const char* ptr, Darray& Energy, std::vector<bool>& EnergyExists) const {
  //mprintf("DBG: [%s]\n", ptr);
  if (ptr == 0 || ptr[0] == '|') return 0;
  const char* beg = ptr;
  //          111111111122222222223
  //0123456789012345678901234567890
  // NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =   435.99  PRESS =-10207.6
  bool eol = false;
  while (!eol) {
    // Skip leading whitespace
    while (*beg == ' ' && *beg != '\0') ++beg;
    if (*beg == '\0') {
      // Line is blank or no more terms. Bail out.
      break;
    }
    //mprintf("DBG: beg= %c\n", *beg);
    // Search for next '='
    const char* eq = beg + 1;
    while (*eq != '=' && *eq != '\0') ++eq;
    if (*eq == '\0')
      eol = true;
    else {
      // Search for end token. Start just after '='.
      const char* val = eq + 1;
      // Skip leading whitespace
      while (*val == ' ' && *val != '\0') ++val;
      if (*val == '\0') {
        eol = true;
        mprintf("Warning: EOL encountered before energy term could be read.\n");
        return 1;
      } else {
        //mprintf("DBG: val= %c\n", *val);
        // Search for next whitespace or line end.
        const char* end = val + 1;
        while (*end != ' ' && *end != '\0' && *end != '\n' && *end != '\r') ++end;
        // Term is now complete. Convert.
        std::string valstr(val, end);
        //mprintf("DBG: valstr= '%s'\n", valstr.c_str());
        std::string termName = NoTrailingWhitespace(std::string(beg,eq));
        FieldType Eindex = getTermIdx(termName);
        if (Eindex != N_FIELDTYPES) {
          if (!validDouble(valstr)) {
            mprintf("Warning: Invalid number detected: %s = %s\n", termName.c_str(), valstr.c_str());
          } else {
            //mprintf("DBG: %s = %s\n", termName.c_str(), valstr.c_str());
            Energy[Eindex] = atof( valstr.c_str() );
            EnergyExists[Eindex] = true;
          }
        }
        beg = end;
      }
    }
  } // END loop over line
  
  return 0;
}

