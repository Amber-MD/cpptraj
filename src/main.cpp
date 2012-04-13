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
#define CPPTRAJ_VERSION_STRING "V12.1b"
#define CPPTRAJ_INTERNAL_VERSION "V3.2.0b"
#endif

// Usage()
/// Print command line usage.
static void Usage(char *programName) {
  mprinterr("Usage: %s [-p <Top1>, -p <Top2>, ...] [-i <Input>] [-debug <N>]\n",programName);
  mprinterr("       %s <Top1> <Input>\n",programName);
  mprinterr("       Additional options:\n");
  mprinterr("         --help, -help: print usage information and exit.\n");
  mprinterr("         -V, --version: print version information and exit.\n");
  mprinterr("             --defines: print list of defines used in compilation.\n");
}

// ProcessInputStream()
/// Process input from the file specified by filename. 
/** If filename is NULL process input from STDIN. Set up an input line that 
  * will be converted to an argument list and processed by the 
  * Cpptraj::Dispatch routine.
  * Leading and consectuive whitespace is skipped. \n or NULL executes command.
  * 'go' or EOF ends input read. lines ending with \ continue to the next line.
  */
static int ProcessInputStream(char *inputFilename, Cpptraj &State) {
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
      while (ptr!='\n' && ptr!=EOF && ptr!='\0') ptr=(char)fgetc(infile); 
    // newline, NULL, or EOF terminates the line
    if (ptr=='\n' || ptr=='\0' || ptr==EOF) {
      // If no chars in string continue
      if (inputLine.empty()) continue;
      // If "go" then done reading input
      if (inputLine.compare("go")==0) break;
      // If "quit" then abort 
      if (inputLine.compare("quit")==0) {
        if (!isStdin) fclose(infile);
        return 1;
      }
      // Print the input line that will be sent to dispatch
      mprintf("  [%s]\n",inputLine.c_str());
      // Call Dispatch to convert input to arglist and process.
      State.Dispatch((char*)inputLine.c_str());
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
      if (ptr == '\n') continue;
      inputLine += "\\";
      /*while ( (ptr = (char)fgetc(infile))!='\n' ) 
        if ( ptr == EOF ) break;
      // NOTE: Insert a space into InputLine?
      continue;*/
    }
    // Skip any line beginning with # character
    if (i==0 && ptr=='#') {
      while ( (ptr = (char)fgetc(infile))!='\n' ) 
        if ( ptr == EOF ) break;
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
 * \return 1 if unrecognized input on command line, 2 if ProcessInputStream indicates 
 *         we should just quit.
 */
static int ProcessCmdLineArgs(int argc, char **argv, Cpptraj &State) {
  char *inputFilename = NULL;
  int debug=0; 

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
      //ProcessInputFile(argv[i]);
      inputFilename=argv[i];

    // -debug: Set overall debug level
    } else if (strcmp(argv[i],"-debug")==0 && i+1!=argc) {
      i++;
      debug = atoi(argv[i]);
      State.SetGlobalDebug( debug );

    // Print information on compiler defines used and exit
    } else if (strcmp(argv[i],"--defines")==0) {
      mprintf("Compiled with:");
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
#ifdef NO_PTRAJ_ANALYZE
      mprintf(" -DNO_PTRAJ_ANALYZE");
#endif
      return 2;

    // The following 2 are for backwards compatibility with PTRAJ
    // Position 1: TOP file
    } else if (i==1) {
      State.AddParm(argv[i]);
    // Position 2: INPUT file
    } else if (i==2) {
      inputFilename=argv[i];

    // Unrecognized
    } else {
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      return 1;
    }
  }

  if ( ProcessInputStream(inputFilename,State) ) return 2;

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
  mprintf("    _|_/\\_|_/\\_|_/\\_|_\n\n");
#ifdef MPI
  mprintf("Running on %i processors\n\n",worldsize);
#endif

  err = ProcessCmdLineArgs(argc,argv,State);
  switch ( err ) {
    case 0 : State.Run(); break;
    case 1 : Usage(argv[0]); break;
  }

  parallel_end();

  mprintf("\n");
  return 0;
}

