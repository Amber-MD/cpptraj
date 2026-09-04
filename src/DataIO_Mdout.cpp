#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include "DataIO_Mdout.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToDouble
#include "DataSet_double.h"
#include "AmberEterm.h"

DataIO_Mdout::DataIO_Mdout()
{

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
  using namespace Cpptraj;
  AmberEterm AEterm;
  AmberEterm::Darray Energy = AEterm.AllocEnergyArray();
  std::vector<bool> EnergyExists = AEterm.AllocExistsArray();
  DataSetList::Darray TimeVals;
  DataSetList::DataListType inputSets(AmberEterm::NenergyTerms(), 0);
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
        for (int i = 0; i < AmberEterm::NenergyTerms(); i++) {
          if (EnergyExists[i]) {
            if (inputSets[i] == 0) {
              MetaData md( dsname, AmberEterm::Ename(i) );
              md.SetLegend( dsname + "_" + AmberEterm::Ename(i) );
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
      sscanf(ptr, " %i %lE %lE %lE", &minStep, Eptr+AmberEterm::EPTOT, Eptr+AmberEterm::RMS, Eptr+AmberEterm::GMAX);
      EnergyExists[AmberEterm::EPTOT] = true;
      EnergyExists[AmberEterm::RMS] = true;
      EnergyExists[AmberEterm::GMAX] = true;
      ptr = buffer.Line();
    }
    // Tokenize line, scan through until '=' is reached; value after is target.
    if (AEterm.GetAmberEterms(buffer.CurrentLine(), Energy, EnergyExists))
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
