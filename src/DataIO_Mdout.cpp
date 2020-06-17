#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include "DataIO_Mdout.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToDouble
#include "DataSet_double.h"

DataIO_Mdout::DataIO_Mdout() {
  // Populate the term name to index map. In some cases, multiple term names
  // map to the same index.
  termIdxMap_.insert(NameIdxPair("Etot", ETOT));
  termIdxMap_.insert(NameIdxPair("EPtot", EPTOT));
  termIdxMap_.insert(NameIdxPair("GMAX", GMAX)); // Not necessary?
  termIdxMap_.insert(NameIdxPair("BOND", BOND));
  termIdxMap_.insert(NameIdxPair("ANGLE", ANGLE));
  termIdxMap_.insert(NameIdxPair("DIHED", DIHED));
  termIdxMap_.insert(NameIdxPair("VDWAALS", VDWAALS));
  termIdxMap_.insert(NameIdxPair("EEL", EEL));
  termIdxMap_.insert(NameIdxPair("EELEC", EEL));
  termIdxMap_.insert(NameIdxPair("EGB", EGB));
  termIdxMap_.insert(NameIdxPair("EPB", EPB));
  termIdxMap_.insert(NameIdxPair("ECAVITY", ECAVITY));
  termIdxMap_.insert(NameIdxPair("EDISPER", EDISPER));
  termIdxMap_.insert(NameIdxPair("1-4 VDW", VDW14));
  termIdxMap_.insert(NameIdxPair("1-4 NB", VDW14));
  termIdxMap_.insert(NameIdxPair("1-4 EEL", EEL14));
  termIdxMap_.insert(NameIdxPair("RESTRAINT", RESTRAINT));
  termIdxMap_.insert(NameIdxPair("EAMBER", EAMBER));
  termIdxMap_.insert(NameIdxPair("Density", DENSITY));
  termIdxMap_.insert(NameIdxPair("RMS", RMS)); // Not necessary?
  termIdxMap_.insert(NameIdxPair("EKtot", EKTOT));
  termIdxMap_.insert(NameIdxPair("ESURF", ESURF));
  termIdxMap_.insert(NameIdxPair("EAMD_BOOST", EAMD_BOOST));
  termIdxMap_.insert(NameIdxPair("VOLUME", VOLUME));
  termIdxMap_.insert(NameIdxPair("TEMP(K)", TEMP));
  termIdxMap_.insert(NameIdxPair("PRESS", PRESS));
  termIdxMap_.insert(NameIdxPair("DV/DL", DVDL));
}

// DataIO_Mdout::ID_DataFormat()
bool DataIO_Mdout::ID_DataFormat(CpptrajFile& infile) {
  if (infile.OpenFile()) return false;
  bool isMdout = false;
  std::string line = infile.GetLine();
  if (line[0] == '\n') {
    line = infile.GetLine();
    if (line.compare(0, 15, "          -----") == 0) {
      line = infile.GetLine();
      if (line.compare(0, 15, "          Amber") == 0)
        isMdout = true;
    }
  }
  infile.CloseFile();
  return isMdout;
}

static inline int EOF_ERROR() {
  mprinterr("Error: Unexpected EOF in MDOUT file.\n");
  return 1;
}

/** Names corresponding to FieldType. */
const char* DataIO_Mdout::Enames_[] = {
  "Etot",   "EPtot",  "GMAX",  "BOND",
  "ANGLE",  "DIHED",  "VDW",   "EELEC",      "EGB",     "EPB", "ECAVITY", "EDISPER",
  "VDW1-4", "EEL1-4", "RST",   "EAMBER",     "Density",
  "RMS",    "EKtot",  "ESURF", "EAMD_BOOST", "VOLUME",  "TEMP",
  "PRESS",  "DVDL",   0
};

/** \return FieldType corresponding to given term name, or N_FIELDTYPES if
  *         not recognized.
  */
DataIO_Mdout::FieldType DataIO_Mdout::getTermIdx(std::string const& name) const {
  NameIdxMap::const_iterator it = termIdxMap_.find( name );
  if (it == termIdxMap_.end()) {
    return (FieldType)N_FIELDTYPES;
  } else {
    return (FieldType)it->second;
  }
}

/** Parse the given line for energy terms of format <name>=<value>. */
int DataIO_Mdout::GetAmberEterms(const char* ptr, Darray& Energy, std::vector<bool>& EnergyExists) {
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

// DataIO_Mdout::ReadData()
int DataIO_Mdout::ReadData(FileName const& fname,
                            DataSetList& datasetlist, std::string const& dsname)
{
  mprintf("\tReading from mdout file: %s\n", fname.full());
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  const char* ptr = buffer.Line();
  if (ptr == 0) {
    mprinterr("Error: Nothing in MDOUT file: %s\n", fname.full());
    return 1;
  }
  // ----- PARSE THE INPUT SECTION ----- 
  int imin = -1;           // imin for this file
  const char* Trigger = 0; // Trigger for storing energies, must be 8 chars long.
  int frame = 0;           // Frame counter for this file
  double dt = 1.0;         // Timestep for this file (MD)
  double t0 = 0.0;         // Initial time for this file (MD)
  int ntpr = 1;            // Value of ntpr
  int irest = 0;           // Value of irest
  while ( ptr != 0 && strncmp(ptr, "   2.  CONTROL  DATA", 20) != 0 )
    ptr = buffer.Line();
  if (ptr == 0) return EOF_ERROR();
  // Determine whether this is dynamics or minimization, get dt
  ptr = buffer.Line(); // Dashes 
  ptr = buffer.Line(); // Blank 
  ptr = buffer.Line(); // title line
  while ( strncmp(ptr, "   3.  ATOMIC", 13) != 0 ) 
  {
    ArgList mdin_args( ptr, " ,=" ); // Remove commas, equal signs
    // Scan for stuff we want
    //mprintf("DEBUG:\tInput[%i] %s", mdin_args.Nargs(), mdin_args.ArgLine());
    for (int col=0; col < mdin_args.Nargs(); col += 2) {
      int col1 = col + 1;
      if (mdin_args[col] == "imin") {
        imin = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: imin is %i\n", imin);
        // Set a trigger for printing. For imin5 this is the word minimization.
        // For imin0 or imin1 this is NSTEP.
        if      (imin==0) Trigger = " NSTEP =";
        else if (imin==1) Trigger = "   NSTEP";
        else if (imin==5) Trigger = "minimiza";
        // Since imin0 and imin1 first trigger has no data, set frame 1 lower.
        if (imin==1 || imin==0) frame = -1;
      } else if (mdin_args[col] == "dt") {
        dt = convertToDouble( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: dt is %f\n", dt);
      } else if (mdin_args[col] == "t") {
        t0 = convertToDouble( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: t is %f\n", t0);
      } else if (mdin_args[col] == "ntpr") {
        ntpr = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: ntpr is %i\n", ntpr);
      } else if (mdin_args[col] == "irest") {
        irest = convertToInteger( mdin_args[ col1 ] );
        if (debug_ > 0) mprintf("\t\tMDIN: irest is %i\n", irest);
      }
    }
    ptr = buffer.Line();
    if (ptr == 0) return EOF_ERROR();
  }
  if (Trigger == 0) {
    mprinterr("Error: Could not determine whether MDOUT is md, min, or post-process.\n");
    return 1;
  }
  // ----- PARSE THE ATOMIC ... SECTION -----
  while ( ptr != 0 && strncmp(ptr, "   4.  RESULTS", 14) != 0 )
  {
    ptr = buffer.Line();
    // If run is a restart, set the initial time value.
    if (irest == 1) {
      if (strncmp(ptr, " begin time", 11) == 0) {
        sscanf(ptr, " begin time read from input coords = %lf", &t0);
        if (debug_ > 0) mprintf("\t\tMD restart initial time= %f\n", t0);
      }
    }
  }
  if (ptr == 0) return EOF_ERROR();
  // ----- PARSE THE RESULTS SECTION -----
  bool finalE = false;
  int nstep;
  int minStep = 0; // For imin=1 only
  if (irest == 0)
    nstep = 0;
  else
    nstep = ntpr;
  Darray Energy(N_FIELDTYPES, 0);
  std::vector<bool> EnergyExists(N_FIELDTYPES, false);
  DataSetList::Darray TimeVals;
  DataSetList::DataListType inputSets(N_FIELDTYPES, 0);
  double time = 0.0;
  while (ptr != 0) {
    // Check for end of imin 0 or 1 run; do not record Average and Stdevs
    if ( (imin == 1 && (strncmp(ptr, "                    FINAL", 25) == 0 ||
                        strncmp(ptr, "   5.  TIMINGS",            14) == 0   )) ||
         (imin == 0 && strncmp(ptr, "      A V", 9) == 0))
      finalE = true;
    // Check for '| TI region  2' to prevent reading duplicate energies
    if ( strncmp(ptr, "| TI region  2", 14) == 0 ) {
      while (ptr != 0 && !(ptr[0] == ' ' && ptr[1] == '-'))
        ptr = buffer.Line();
      if (ptr == 0) return EOF_ERROR();
    }
    // Record set for energy post-processing
    if (imin == 5 && strncmp(ptr, "minimizing", 10) == 0)
      nstep = atoi( ptr + 22 );
    // MAIN OUTPUT ROUTINE
    // If the trigger has been reached print output.
    // For imin0 and imin1 the first trigger will have no data.
    // If the end of the file has been reached print then exit.
    if ( strncmp(ptr, Trigger, 8) == 0 || finalE ) {
      if (frame > -1) {
        // Store all energies present.
        for (int i = 0; i < (int)N_FIELDTYPES; i++) {
          if (EnergyExists[i]) {
            if (inputSets[i] == 0) {
              MetaData md( dsname, Enames_[i] );
              md.SetLegend( dsname + "_" + Enames_[i] );
              inputSets[i] = new DataSet_double();
              inputSets[i]->SetMeta( md );
            }
            // Since energy terms can appear and vanish over the course of the
            // mdout file, resize if necessary.
            if (frame > (int)inputSets[i]->Size())
              ((DataSet_double*)inputSets[i])->Resize( frame );
            ((DataSet_double*)inputSets[i])->AddElement( Energy[i] );
          }
        }
        TimeVals.push_back( time );
        nstep += ntpr;
      }
      frame++;
      if (finalE) break;
    }
    // Check for NSTEP in minimization or post-processing. Values will be
    // on the next line. NOTE: NSTEP means something different for imin=5.
    if ((imin == 1 || imin == 5) && strncmp(ptr, "   NSTEP", 8) == 0) {
      ptr = buffer.Line(); // Get next line
      //sscanf(ptr, " %6lf    %13lE  %13lE  %13lE", Energy+NSTEP, Energy+EPtot, Energy+RMS, Energy+GMAX);
      double* Eptr = &(Energy[0]);
      sscanf(ptr, " %i %lE %lE %lE", &minStep, Eptr+EPTOT, Eptr+RMS, Eptr+GMAX);
      EnergyExists[EPTOT] = true;
      EnergyExists[RMS] = true;
      EnergyExists[GMAX] = true;
      ptr = buffer.Line();
    }
    // Tokenize line, scan through until '=' is reached; value after is target.
    if (GetAmberEterms(buffer.CurrentLine(), Energy, EnergyExists))
      mprintf("Warning: Issue parsing line %i\n", buffer.LineNumber());
    // Set time
    switch (imin) {
      case 5: time = (double)nstep + t0; break;
      case 1: time = (double)minStep + t0; break;
      case 0: time = ((double)nstep * dt) + t0; break;
    }
    // Read in next line
    ptr = buffer.Line();
  }
  mprintf("\t%i frames\n", frame);
  buffer.CloseFile();
  std::string Xlabel;
  if      (imin == 5) Xlabel.assign("Set");
  else if (imin == 1) Xlabel.assign("Nstep");
  else                Xlabel.assign("Time"); // imin == 0
  if (datasetlist.AddOrAppendSets( Xlabel, TimeVals, inputSets )) return 1;
  return 0;
}
