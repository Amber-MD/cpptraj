#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include "DataIO_Mdout.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataIO_Mdout::DataIO_Mdout() {}

static inline int EOF_ERROR() {
  mprinterr("Error: Unexpected EOF in MDOUT file.\n");
  return 1;
}

const char* DataIO_Mdout::Enames[] = {
  "NSTEP",  "Etot",   "EPtot", "GMAX",   "BOND", 
  "ANGLE",  "DIHED",  "VDW",   "EELEC",  "EGB",
  "VDW1-4", "EEL1-4", "RST",   "EAMBER", "Density",
  "RMS",    "EKtot",  "ESURF", "EAMD_BOOST", 0
};

/// \return index of name in Energy[] array, N_FIELDTYPES if not recognized.
DataIO_Mdout::FieldType DataIO_Mdout::getEindex(std::vector<std::string> const& Name) {
  //mprintf("DEBUG:\tgetEindex(%s,%s)\n", Name[0].c_str(), Name[1].c_str());
  if (Name[0]=="NSTEP") return NSTEP;
  if (Name[0]=="Etot")  return Etot;
  if (Name[0]=="EPtot") return EPtot;
  if (Name[0]=="GMAX") return GMAX; // Not necessary?
  if (Name[0]=="BOND") return BOND;
  if (Name[0]=="ANGLE") return ANGLE;
  if (Name[0]=="DIHED") return DIHED;
  if (Name[0]=="VDWAALS") return VDWAALS;
  if (Name[0]=="EEL" || Name[0]=="EELEC") return EEL;
  if (Name[0]=="EGB") return EGB;
  if ((Name[0]=="1-4" && Name[1]=="VDW") || (Name[0]=="1-4" && Name[1]=="NB")) return VDW14;
  if  (Name[0]=="1-4" && Name[1]=="EEL") return EEL14;
  if (Name[0]=="RESTRAINT") return RESTRAINT;
  if (Name[0]=="EAMBER") return EAMBER;
  if (Name[0]=="Density") return Density;
  if (Name[0]=="RMS") return RMS; // Not necessary?
  if (Name[0]=="EKtot") return EKtot;
  if (Name[0]=="ESURF") return ESURF;
  if (Name[0]=="EAMD_BOOST") return EAMD_BOOST;
  return N_FIELDTYPES;
}

// DataIO_Mdout::ReadData()
int DataIO_Mdout::ReadData(std::string const& fname, ArgList& argIn,
                            DataSetList& datasetlist, std::string const& dsname)
{
  std::vector<std::string> mdoutFilenames;
  mdoutFilenames.push_back( fname );
  // Check if more than one mdout name was specified.
  ArgList mdoutnames(argIn.GetStringKey("mdoutfiles"), ",");
  if (!mdoutnames.empty()) {
    for (int i = 0; i < mdoutnames.Nargs(); i++)
      mdoutFilenames.push_back( mdoutnames[i] );
  }
  mprintf("\tReading from mdout files:");
  for (std::vector<std::string>::const_iterator it = mdoutFilenames.begin();
                                          it != mdoutFilenames.end(); ++it)
    mprintf(" %s", (*it).c_str());
  mprintf("\n");

  // ----- CREATE DATASETS FOR ENERGIES -----
  // TODO: Set it up so blank data sets can be REMOVED.
  std::vector<DataSet*> Esets( N_FIELDTYPES, 0 );
  for (int i = 1; i < (int)N_FIELDTYPES; i++) { // Do not store NSTEP
    Esets[i] = datasetlist.AddSetAspect( DataSet::DOUBLE, dsname, Enames[i] );
    // Make legend same as aspect.
    Esets[i]->SetLegend( Enames[i] );
  }

  // LOOP OVER ALL MDOUT FILES
  BufferedLine buffer;
  double lastx = 0.0;
  int count = 0;
  for (std::vector<std::string>::const_iterator mdoutname = mdoutFilenames.begin();
                                                mdoutname != mdoutFilenames.end();
                                                ++mdoutname)
  {
    mprintf("\t%s\n", (*mdoutname).c_str());
    if (buffer.OpenFileRead( *mdoutname )) return 1;
    // Read first line
    const char* ptr = buffer.Line();
    if (ptr == 0) {
      mprinterr("Error: Nothing in MDOUT file: %s\n", (*mdoutname).c_str());
      return 1;
    }
    int imin = -1;           // imin for this file
    const char* Trigger = 0; // Trigger must be 8 chars long.
    int frame = 0;           // Frame counter for this file
    double dt = 1.0;         // Timestep for this file
    // ----- PARSE THE INPUT SECTION ----- 
    while ( ptr != 0 && strncmp(ptr, " Here is the input", 18) != 0 )
      ptr = buffer.Line();
    if (ptr == 0) return EOF_ERROR();
    // Determine whether this is dynamics or minimization, get dt
    ptr = buffer.Line(); // Blank line
    ptr = buffer.Line(); // title line
    mprintf("DEBUG:\tProcessing MD input: %s", ptr);
    while ( strncmp(ptr, "   1.  RESOURCE", 15) != 0 ) 
    {
      ArgList mdin_args( ptr, " ,=" ); // Remove commas, equal signs
      // Scan for stuff we want
      //mprintf("DEBUG:\tInput[%i] %s", mdin_args.Nargs(), mdin_args.ArgLine());
      for (int col=0; col < mdin_args.Nargs(); col += 2) {
        int col1 = col + 1;
        if (mdin_args[col] == "imin") {
          imin = mdin_args.IntegerAt( col1 );
          mprintf("\tMDIN: imin is %i\n", imin);
          // Set a trigger for printing. For imin5 this is the word minimization.
          // For imin0 or imin1 this is NSTEP.
          if      (imin==0) Trigger = " NSTEP =";
          else if (imin==1) Trigger = "   NSTEP";
          else if (imin==5) Trigger = "minimiza";
          // Since imin0 and imin1 first trigger has no data, set frame 1 lower.
          if (imin==1 || imin==0) frame = -1;
        } else if (mdin_args[col] == "dt") {
          dt = mdin_args.DoubleAt( col1 );
          mprintf("\tMDIN: dt is %f\n", dt);
        }
      }
      ptr = buffer.Line();
    }
    if (Trigger == 0) {
      mprinterr("Error: Could not determine whether MDOUT is md, min, or post-process.\n");
      return 1;
    }

    // ----- PARSE THE RESULTS SECTION -----
    while ( ptr != 0 && strncmp(ptr, "   4.  RESULTS", 14) != 0 )
      ptr = buffer.Line();
    if (ptr == 0) return EOF_ERROR();
/*    CpptrajFile TestOut; // DEBUG
    TestOut.OpenWrite("TestOut.dat"); // DEBUG
    bool printHeader = true; // DEBUG*/
    bool finalE = false;
    int set = 0;
    double Energy[N_FIELDTYPES];
    std::fill( Energy, Energy+N_FIELDTYPES, 0.0 );
    std::vector<bool> EnergyExists(N_FIELDTYPES, false);
    std::vector<std::string> Name(2);
    double time = 0.0;
    while (ptr != 0) {
      // Check for end of imin 0 or 1 run; do not record Average and Stdevs
      if ( (imin == 1 && (strncmp(ptr, "                    FINAL", 25) == 0 ||
                          strncmp(ptr, "   5.  TIMINGS",            14) == 0   )) ||
           (imin == 0 && strncmp(ptr, "      A V", 9) == 0))
        finalE = true;
      // Record set for energy post-processing
      if (imin == 5 && strncmp(ptr, "minimizing", 10) == 0)
        set = atoi( ptr + 22 );
      // MAIN OUTPUT ROUTINE
      // If the trigger has been reached print output.
      // For imin0 and imin1 the first trigger will have no data.
      // If the end of the file has been reached print then exit.
      if ( strncmp(ptr, Trigger, 8) == 0 || finalE ) {
        if (frame > -1) {
          // Data storage should go here
          for (int i = 1; i < (int)N_FIELDTYPES; i++) // skip NSTEP
            if (EnergyExists[i]) Esets[i]->Add( count, Energy + i );
/*          // DEBUG
          if (printHeader) {
            TestOut.Printf("%-14s", "#Time");
            for (int i = 1; i < (int)N_FIELDTYPES; i++) // skip NSTEP
              if (EnergyExists[i]) TestOut.Printf(" %14s", Enames[i]);
            TestOut.Printf("\n");
            printHeader = false;
          }
          TestOut.Printf(" %14.4f", time);
          for (int i = 1; i < (int)N_FIELDTYPES; i++) // skip NSTEP
            if (EnergyExists[i]) TestOut.Printf(" %14.4f", Energy[i]);
          TestOut.Printf("\n");
          // DEBUG*/
          count++;
        }
        frame++;
        if (finalE) break;
      }
      // Check for NSTEP in minimization or post-processing. Values will be
      // on the next line. NOTE: NSTEP means something different for imin=5.
      if ((imin == 1 || imin == 5) && strncmp(ptr, "   NSTEP", 8) == 0) {
        ptr = buffer.Line(); // Get next line
        sscanf(ptr, " %6lf    %13lE  %13lE  %13lE", Energy+NSTEP, Energy+EPtot,
               Energy+RMS, Energy+GMAX);
        EnergyExists[NSTEP] = true;
        EnergyExists[EPtot] = true;
        EnergyExists[RMS] = true;
        EnergyExists[GMAX] = true;
        ptr = buffer.Line();
      }
      // Tokenize line, scan through until '=' is reached; value after is target.
      int ntokens = buffer.TokenizeLine(" ");
      if (ntokens > 0) {
        int nidx = 0;
        Name[0].clear();
        Name[1].clear();
        for (int tidx = 0; tidx < ntokens; tidx++) {
          const char* tkn = buffer.NextToken();
          if (tkn[0] == '=') {
            FieldType Eindex = getEindex(Name);
            tkn = buffer.NextToken();
            ++tidx;
            if (Eindex != N_FIELDTYPES) {
              Energy[Eindex] = atof( tkn );
              EnergyExists[Eindex] = true;
            }
            nidx = 0;
            Name[0].clear();
            Name[1].clear();
          } else {
            if (nidx > 1) break; // Two tokens, no '=' found. Not an E line.
            Name[nidx++].assign( tkn );
          }
        }
      }
      // Set time
      switch (imin) {
        case 5: time = (double)set; break;
        case 1: time = Energy[0] + lastx; break;
        case 0: time = (Energy[0] * dt) + lastx; break;
      }
      // Read in next line
      ptr = buffer.Line();
    }
/*    TestOut.CloseFile(); // DEBUG*/
    mprintf("\t%i frames\n", frame);
    lastx = time;
    buffer.CloseFile();
  } // END loop over mdout files
  
  // ----- REMOVE EMPTY DATASETS -----
  for (int i = 1; i < (int)N_FIELDTYPES; i++) { // Do not store NSTEP
    if (Esets[i]->Empty())
      datasetlist.erase( Esets[i] );
  }
  return 0;
}
