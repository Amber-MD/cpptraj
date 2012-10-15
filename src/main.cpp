/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * 2010 Daniel R. Roe
 */
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
#include <cstring>
#include <cstdio>
#include <cstdlib> // atoi
#ifndef CPPTRAJ_VERSION_STRING
#define CPPTRAJ_VERSION_STRING "V13.8.2b"
#define CPPTRAJ_INTERNAL_VERSION "V3.6.3b"
#endif

// Usage()
/// Print command line usage.
static void Usage(const char *programName) {
  mprinterr("\nUsage: %s [-p <Top1>, -p <Top2>, ...] [-i <Input1>, -i <Input2>, ...]\n",
            programName);
  mprinterr(  "       %s <Top> <Input>\n",programName);
  mprinterr(  "       Additional options:\n");
  mprinterr(  "         --help, -help : Print usage information and exit.\n");
  mprinterr(  "         -V, --version : Print version information and exit.\n");
  mprinterr(  "         --defines     : Print list of defines used in compilation.\n");
  mprinterr(  "         -debug <N>    : Set global debug level.\n");
}

static inline bool EndChar(char ptr) {
  if (ptr=='\n' || ptr=='\r' || ptr=='\0' || ptr==EOF) return true;
  return false;
}

// ProcessInputStream()
/// Process input from the file specified by filename. 
/** If filename is NULL process input from STDIN. Set up an input line that 
  * will be converted to an argument list and processed by the 
  * Cpptraj::Dispatch routine.
  * Leading and consectuive whitespace is skipped. \n or NULL executes command.
  * 'go' or EOF ends input read. lines ending with \ continue to the next line.
  * \return 0 on success, 1 on error, 2 on exit.
  */
static int ProcessInputStream(const char *inputFilename, Cpptraj& State) {
  FILE *infile;
  std::string inputLine;

  bool isStdin = false;
  // Open input file or STDIN
  if (inputFilename==NULL) {
    // Do not allow read from STDIN when > 1 process
    if (worldsize > 1) {
      mprintf("Error: Reading from STDIN not allowed with more than 1 thread.\n");
      mprintf("       To run cpptraj in parallel please use an input file.\n");
      return 1;
    }
    mprintf("INPUT: Reading Input from STDIN, type \"go\" to run, \"quit\" to exit:\n");
    infile=stdin;
    isStdin=true;
  } else {
    rprintf("INPUT: Reading Input from file %s\n",inputFilename);
    if ( (infile=fopen(inputFilename,"r"))==NULL ) {
      rprintf("Error: Could not open input file %s\n",inputFilename);
      return 1;
    }
  }

  // Read in each line of input. Newline or NULL terminates. \ continues line.
  int i = 0; // Index in inputLine
  char lastchar = '0';
  char ptr = 0;
  if (isStdin) fprintf(stdout,"> ");
  while ( ptr != EOF ) {
    //if (prompt) {fprintf(stdout,"> "); prompt=false;}
    ptr = (char)fgetc(infile);
    //fprintf(stdout,"DEBUG: %i %c %i\n",i,ptr,ptr);
    // If '#' is encountered, skip the rest of the line
    if (ptr=='#')
      while (!EndChar(ptr)) ptr=(char)fgetc(infile); 
    // newline, NULL, or EOF terminates the line
    if (EndChar(ptr)) {
      // If no chars in string continue
      if (inputLine.empty()) {
        if (isStdin) fprintf(stdout,"> ");
        continue;
      }
      // If "go" then done reading input
      if (inputLine.compare("go")==0) break;
      // If "quit" then abort 
      if (inputLine.compare("quit")==0) {
        if (!isStdin) fclose(infile);
        return 2;
      } else if (inputLine.compare(0,2,"ls")==0)
        system(inputLine.c_str());
      else if (inputLine.compare("pwd")==0)
        system("pwd");
      else {
        // Print the input line that will be sent to dispatch
        mprintf("  [%s]\n",inputLine.c_str());
        // Call Dispatch to convert input to arglist and process.
        State.Dispatch(inputLine.c_str());
      }
      // Reset Input line
      inputLine.clear();
      i=0;
      if (isStdin) fprintf(stdout,"> ");
      continue;
    }
    // Any consecutive whitespace is skipped
    if (i>0) lastchar=inputLine[i-1];
    if (isspace(ptr) && isspace(lastchar)) continue;
    // Skip leading whitespace
    if (i==0 && isspace(ptr)) {
      while ( (ptr = (char)fgetc(infile))!=EOF )
        if ( !isspace(ptr) ) break;
    } 
    // Backslash followed by newline continues to next line. Otherwise backslash
    // followed by next char will be inserted. 
    if (ptr=='\\') {
      ptr = (char)fgetc(infile);
      if ( ptr == EOF ) break;
      if (ptr == '\n' || ptr == '\r') continue;
      inputLine += "\\";
      /*while ( (ptr = (char)fgetc(infile))!='\n' ) 
        if ( ptr == EOF ) break;
      // NOTE: Insert a space into InputLine?
      continue;*/
    }
    // Skip any line beginning with # character
    if (i==0 && ptr=='#') {
      while ( !EndChar(ptr) ) ptr = (char)fgetc(infile);
      if (isStdin) fprintf(stdout,"> ");
      continue;
    }
    // Add character to input line
    inputLine += ptr;
    ++i;
  }

  if (!isStdin) fclose(infile);
  return 0;
}

// ProcessCmdLineArgs()
/// Process arguments on the command line. 
/** Process the input file last
 * despite its position on the command line to allow any prmtops to
 * load.
 * \return 0 on success
 * \return 1 if unrecognized input on command line
 * \return 2 if ProcessInputStream indicates we should just quit.
 * \return 3 if input from STDIN needed
 */
static int ProcessCmdLineArgs(int argc, char **argv, Cpptraj &State) {
  int debug=0;
  bool hasInput = false;

  if (argc == 1) return 3; // No command line args, interactive

  for (int i=1; i<argc; i++) {
    // --help, -help: Print usage and exit
    if (strcmp(argv[i],"--help")==0 || strcmp(argv[i],"-help")==0) {
      return 1;

    // -V, --version: Print version number and exit
    // NOTE: version number is already printed - should order be changed?
    } else if (strcmp(argv[i],"-V")==0 || strcmp(argv[i],"--version")==0) {
      return 2;

    // -p: Topology file
    } else if (strcmp(argv[i],"-p")==0 && i+1!=argc) {
      i++;
      if (debug>0) mprintf("Adding topology file from command line: %s\n", argv[i]);
      State.AddParm(argv[i]);

    // -i: Input file
    } else if (strcmp(argv[i],"-i")==0 && i+1!=argc) {
      i++;
      if (ProcessInputStream(argv[i], State)) return 2;
      hasInput = true;

    // -debug: Set overall debug level
    } else if (strcmp(argv[i],"-debug")==0 && i+1!=argc) {
      i++;
      debug = atoi(argv[i]);
      State.SetGlobalDebug( debug );

    // Print information on compiler defines used and exit
    } else if (strcmp(argv[i],"--defines")==0) {
      mprintf("\nCompiled with:");
#ifdef DEBUG
      mprintf(" -DDEBUG");
#endif
#ifdef HASBZ2
      mprintf(" -DHASBZ2");
#endif
#ifdef HASGZ
      mprintf(" -DHASGZ");
#endif
#ifdef BINTRAJ
      mprintf(" -DBINTRAJ");
#endif
#ifdef MPI
      mprintf(" -DMPI");
#endif
#ifdef _OPENMP
      mprintf(" -D_OPENMP");
#endif
#ifdef NO_MATHLIB
      mprintf(" -DNO_MATHLIB");
#endif
      mprintf("\n");
      return 2;

    // The following 2 are for backwards compatibility with PTRAJ
    // Position 1: TOP file
    } else if (i==1) {
      State.AddParm(argv[i]);
    // Position 2: INPUT file
    } else if (i==2) {
      if (ProcessInputStream(argv[i], State)) return 2;
      hasInput = true;

    // Unrecognized
    } else {
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      return 1;
    }
  }

  if (!hasInput) return 3;

  return 0;
}

// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
/** 1. Call parallel Init (does nothing if not a parallel run)
 * 2. Process input from command line/inputfiles/stdin
 * 3. Run
 */
int main(int argc, char **argv) {
  Cpptraj State;
  int err;

  // Parallel Init: NOTE Should check for err
  parallel_init(argc,argv);

  mprintf("\nCPPTRAJ: Trajectory Analysis. %s\n",CPPTRAJ_VERSION_STRING);
  mprintf("    ___  ___  ___  ___\n");
  mprintf("     | \\/ | \\/ | \\/ | \n");
  mprintf("    _|_/\\_|_/\\_|_/\\_|_\n");
#ifdef MPI
  mprintf("Running on %i processors\n",worldsize);
#endif
  err = ProcessCmdLineArgs(argc,argv,State);
  if (err == 3) // More input needed, go to interactive mode
    err = ProcessInputStream(NULL, State);
  switch ( err ) {
    case 0 : State.Run(); break;
    case 1 : Usage(argv[0]); break;
  }

  parallel_end();

  mprintf("\n");
  return 0;
}

