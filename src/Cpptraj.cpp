#include <cstdio> // stdin, fileno
#include <unistd.h> // isatty
#include "Cpptraj.h"
#include "CpptrajStdio.h"
#include "CommandList.h"
#include "ReadLine.h"
#include "Version.h"

void Cpptraj::Usage() {
  mprinterr("\n"
            "Usage: cpptraj [-p <Top0>] [-i <Input0>] [-y <trajin>] [-x <trajout>]\n"
            "               [-h | --help] [-V | --version] [--defines] [-debug <#>]\n"
            "               [--interactive] [--log <logfile>]\n"
            "       cpptraj <Top> <Input>\n"
            "\t-p <Top0>      : Load <Top0> as a topology file. May be specified more than once.\n"
            "\t-i <Input0>    : Read input from <Input0>. May be specified more than once.\n"
            "\t-y <trajin>    : Read from trajectory file <trajin>; same as input 'trajin <trajin>'.\n"
            "\t-x <trajout>   : Write trajectory file <trajout>; same as input 'trajout <trajout>'.\n"
            "\t-h | --help    : Print command line help and exit.\n"
            "\t-V | --version : Print version and exit.\n"
            "\t--defines      : Print compiler defines and exit.\n"
            "\t-debug <#>     : Set global debug level to <#>; same as input 'debug <#>'.\n"
            "\t--interactive  : Force interactive mode.\n"
            "\t--log <logfile>: Record commands to <logfile> (interactive mode only). Default is 'cpptraj.log'.\n");
}

void Cpptraj::Intro() {
  mprintf("\nCPPTRAJ: Trajectory Analysis. %s\n"
          "    ___  ___  ___  ___\n     | \\/ | \\/ | \\/ | \n    _|_/\\_|_/\\_|_/\\_|_\n",
          CPPTRAJ_VERSION_STRING);
# ifdef MPI
  mprintf("Running on %i threads\n",worldsize);
# endif
}

void Cpptraj::Finalize() {
  mprintf("--------------------------------------------------------------------------------\n"
    "To cite CPPTRAJ use:\n"
    "Daniel R. Roe and Thomas E. Cheatham, III, \"PTRAJ and CPPTRAJ: Software for\n"
    "  Processing and Analysis of Molecular Dynamics Trajectory Data\". J. Chem.\n"
    "  Theory Comput., 2013, 9 (7), pp 3084-3095.\n");
}

int Cpptraj::RunCpptraj(int argc, char** argv) {
  int err = 0;
  Cpptraj::Intro();
  Mode cmode = ProcessCmdLineArgs(argc, argv);
  if ( cmode == BATCH ) {
    // If run has not yet been called, run now.
    if (State_.Nrun() < 1)
      err = State_.Run();
  } else if ( cmode == INTERACTIVE ) {
    err = Interactive();
  } else if ( cmode == ERROR ) {
    err = 1;
  }
  if (err == 0)
    Cpptraj::Finalize();
  mprintf("\n");
  return err;
}

/** Read command line args. */
Cpptraj::Mode Cpptraj::ProcessCmdLineArgs(int argc, char** argv) {
  bool hasInput = false;
  bool interactive = false;
  std::vector<std::string> inputFiles;
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    if ( arg == "--help" || arg == "-h" ) {
      // --help, -help: Print usage and exit
      Usage();
      return QUIT;
    }
    if ( arg == "-V" || arg == "--version" ) {
      // -V, --version: Print version number and exit
      //mprintf("CPPTRAJ: Version %s\n", CPPTRAJ_VERSION_STRING);
      return QUIT;
    }
    if ( arg == "--internal-version" ) {
      // --internal-version: Print internal version number and quit.
      mprintf("Internal version #: %s\n", CPPTRAJ_INTERNAL_VERSION);
      return QUIT;
    }
    if ( arg == "--defines" ) {
      // --defines: Print information on compiler defines used and exit
      mprintf("Compiled with:");
#     ifdef DEBUG
      mprintf(" -DDEBUG");
#     endif
#     ifdef HASBZ2
      mprintf(" -DHASBZ2");
#     endif
#     ifdef HASGZ
      mprintf(" -DHASGZ");
#     endif
#     ifdef BINTRAJ
      mprintf(" -DBINTRAJ");
#     endif
#     ifdef MPI
      mprintf(" -DMPI");
#     endif
#     ifdef _OPENMP
      mprintf(" -D_OPENMP");
#     endif
#     ifdef NO_MATHLIB
      mprintf(" -DNO_MATHLIB");
#     endif
      mprintf("\n");
      return QUIT;
    }
    if ( arg == "--interactive" )
      interactive = true;
    else if ( arg == "-debug" && i+1 != argc) {
      // -debug: Set overall debug level
      ArgList dbgarg( argv[++i] );
      State_.SetListDebug( dbgarg );
    } else if ( arg == "--log" && i+1 != argc)
      // --log: Set up log file for interactive mode
      logfilename_ = argv[++i];
    else if ( arg == "-p" && i+1 != argc) {
      // -p: Topology file
      if (State_.PFL()->AddParmFile( argv[++i] )) return ERROR;
    } else if ( arg == "-y" && i+1 != argc) {
      // -y: Trajectory file in
      ArgList targ("trajin " + std::string(argv[++i]));
      if (State_.AddTrajin( targ )) return ERROR;
    } else if ( arg == "-x" && i+1 != argc) {
      // -x: Trajectory file out
      ArgList targ("trajout " + std::string(argv[++i]));
      if (State_.AddTrajout( targ )) return ERROR;
      hasInput = true; // NOTE: Why?
    } else if ( arg == "-c" && i+1 != argc) {
      // -c: Reference file
      ArgList targ("reference " + std::string(argv[++i]));
      if (State_.AddReference( targ )) return ERROR;
    } else if (arg == "-i" && i+1 != argc) {
      // -i: Input file(s)
      inputFiles.push_back( argv[++i] );
    } else if (arg == "-ms" && i+1 != argc) {
      // -ms: Mask string
      if (State_.MaskString( std::string(argv[++i]))) return ERROR;
      return QUIT;
    } else if ( i == 1 ) {
      // For backwards compatibility with PTRAJ; Position 1 = TOP file
      if (State_.PFL()->AddParmFile( argv[i])) return ERROR;
    } else if ( i == 2 ) {
      // For backwards compatibility with PTRAJ; Position 2 = INPUT file
      inputFiles.push_back( argv[i] );
    } else {
      // Unrecognized
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      Usage();
      return QUIT;
    }
  }
  // Process all input files specified on command line.
  if ( !inputFiles.empty() ) {
    hasInput = true;
    for (std::vector<std::string>::const_iterator inputFilename = inputFiles.begin();
                                                  inputFilename != inputFiles.end();
                                                  ++inputFilename)
    {
      CommandList::RetType c_err = CommandList::ProcessInput( State_, *inputFilename );
      if (c_err == CommandList::C_ERR) return ERROR;
      if (c_err == CommandList::C_QUIT) return QUIT;
    }
  }
  // Determine whether to enter interactive mode
  if (!hasInput || interactive) {
    // Test if input is really from a console
    if ( isatty(fileno(stdin)) )
      return INTERACTIVE;
    else {
      // "" means read from STDIN
      CommandList::RetType c_err = CommandList::ProcessInput( State_, "" ); 
      if (c_err == CommandList::C_ERR) return ERROR;
      if (c_err == CommandList::C_QUIT) return QUIT;
    }
  }
  return BATCH;
}

// Cpptraj::Interactive()
int Cpptraj::Interactive() {
  ReadLine inputLine;
  // By default when interactive do not exit on errors
  State_.SetNoExitOnError();
  // Open log file. If no name has been set, use default.
  CpptrajFile logfile_;
  if (!logfilename_.empty())
    logfile_.OpenWrite(logfilename_);
  else
    logfile_.OpenAppend("cpptraj.log");
  CommandList::RetType readLoop = CommandList::C_OK;
  while ( readLoop != CommandList::C_QUIT ) {
    if (inputLine.GetInput()) {
      // EOF (Ctrl-D) specified. If there are actions/analyses queued, ask 
      // user if they really want to quit.
      if (!State_.Empty()) {
        if (inputLine.YesNoPrompt("EOF (Ctrl-D) specified but there are actions/analyses"
                                  " queued. Really quit? [y/n]> "))
          break;
      } else
        break;
    }
    if (!inputLine.empty()) {
      readLoop = CommandList::Dispatch( State_, *inputLine );
      if (logfile_.IsOpen())
        logfile_.Printf("%s\n", inputLine.c_str());
    }
  }
  logfile_.CloseFile();
  if (readLoop == CommandList::C_ERR) return 1;
  return 0;
}
