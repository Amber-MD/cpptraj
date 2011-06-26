#include <cstdio>
#include <cstdlib> // for atoi
#include <cstring>
#include <cctype>
#include "PtrajState.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"

// Constructor
PtrajState::PtrajState() {
  TotalErrors=0; 
  debug=0;
  showProgress=true;
}

// Destructor
PtrajState::~PtrajState() {
  for (std::list<ArgList*>::iterator it=DF_Args.begin(); it!=DF_Args.end(); it++) 
    delete (*it); 
}

/*
 * PtrajState::SetGlobalDebug()
 * Set the debug level for all components of PtrajState.
 */
void PtrajState::SetGlobalDebug(int debugIn) {
  debug = debugIn;
  rprintf("DEBUG LEVEL SET TO %i\n",debug);
  trajinList.SetDebug(debug);
  referenceList.SetDebug(debug);
  trajoutList.SetDebug(debug);
  parmFileList.SetDebug(debug);
  ptrajActionList.SetDebug(debug);
  DFL.SetDebug(debug);
}

/* 
 * PtrajState::ProcessCmdLineArgs()
 * Process arguments on the command line. Process the input file last
 * despite its position on the command line to allow any prmtops to
 * load.
 * Return 1 if unrecognized input on command line.
 * Return 2 if ProcessInputStream indicates we should just quit.
 */
int PtrajState::ProcessCmdLineArgs(int argc, char **argv) {
  int i;
  char *inputFilename;

  inputFilename=NULL;
  for (i=1; i<argc; i++) {
    // --help, -help: Print usage and exit
    if (strcmp(argv[i],"--help")==0 || strcmp(argv[i],"-help")==0) {
      return 1;

    // -p: Topology file
    } else if (strcmp(argv[i],"-p")==0 && i+1!=argc) {
      i++;
      if (debug>0) mprintf("Adding topology file from command line: %s\n", argv[i]);
      parmFileList.Add(argv[i]);

    // -i: Input file
    } else if (strcmp(argv[i],"-i")==0 && i+1!=argc) {
      i++;
      //ProcessInputFile(argv[i]);
      inputFilename=argv[i];

    // -debug: Set overall debug level
    } else if (strcmp(argv[i],"-debug")==0 && i+1!=argc) {
      i++;
      SetGlobalDebug( atoi(argv[i]) );

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
      return 2;

    // The following 2 are for backwards compatibility with PTRAJ
    // Position 1: TOP file
    } else if (i==1) {
      parmFileList.Add(argv[i]);
    // Position 2: INPUT file
    } else if (i==2) {
      inputFilename=argv[i];

    // Unrecognized
    } else {
      mprintf("PtrajState::ProcessArgs: Unrecognized input on command line: %i: %s\n",
              i,argv[i]);
      return 1;
    }
  }

  if ( ProcessInputStream(inputFilename) ) return 2;

  return 0;
}

/*
 * PtrajState::ProcessInputStream()
 * Process input from the file specified by filename. If filename is NULL
 * process input from STDIN. Set up an argument list and execute commands
 * via the Dispatch routine. 
 * Leading and consectuive whitespace is skipped. \n or NULL executes command.
 * go or EOF ends input read. lines ending with \ continue to the next line.
 */
int PtrajState::ProcessInputStream(char *inputFilename) {
  FILE *infile;
  char ptr,lastchar;
  char inputLine[BUFFER_SIZE]; // Careful, could blow this
  int i;
  bool isStdin;

  isStdin=false;
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
  memset(inputLine,' ',BUFFER_SIZE);
  i=0; // Index in inputLine
  lastchar='0';
  ptr=0;
  if (isStdin) fprintf(stdout,"> ");
  while ( ptr != EOF ) {
    ptr = (char)fgetc(infile);
    //fprintf(stdout,"DEBUG: %i %c %i\n",i,ptr,ptr);
    // newline, NULL, or EOF terminates the line
    if (ptr=='\n' || ptr=='\0' || ptr==EOF) {
      inputLine[i]='\0';
      // If no chars in string continue
      if (strlen(inputLine)==0) continue;
      // If "go" then done reading input
      if (strncmp(inputLine,"go",2)==0) break;
      // If "quit" then abort this - only for stdin
      if (isStdin && strncmp(inputLine,"quit",4)==0) return 1;
      mprintf("  [%s]\n",inputLine);
      // Convert the input line to a list of arguments
      A=new ArgList(inputLine," "); // Space delimited only?
      //A->print();
      // Call Dispatch to perform command in arglist
      Dispatch(); 
      // Reset Input line and Argument List
      delete A;
      memset(inputLine,' ',BUFFER_SIZE);
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
    // Forward slash continues to next line. Anything after slash is ignored
    if (ptr=='\\') {
      while ( (ptr = (char)fgetc(infile))!='\n' ) 
        if ( ptr == EOF ) break;
      // NOTE: Insert a space into InputLine?
      continue;
    }
    // Skip any line beginning with # character
    if (i==0 && ptr=='#') {
      while ( (ptr = (char)fgetc(infile))!='\n' ) 
        if ( ptr == EOF ) break;
      if (isStdin) fprintf(stdout,"> ");
      continue;
    }
    inputLine[i++]=ptr;
    // Check to make sure we arent blowing buffer
    if (i==BUFFER_SIZE) {
      rprintf("Error: Input line is greater than BUFFER_SIZE (%u)\n",BUFFER_SIZE);
      if (!isStdin) fclose(infile);
      return 1;
    }
  }

  if (!isStdin) fclose(infile);
  return 0;
}

/* 
 * PtrajState::Dispatch()
 * Send commands to their appropriate classes.
 * The command is tried on each class in turn. If the class rejects command
 * move onto the next one. If command is accepted return.
 * NOTE: Should differentiate between keyword rejection and outright error.
 */
void PtrajState::Dispatch() {
  AtomMask *tempMask;  // For ParmInfo
  AmberParm *tempParm; // For ParmInfo
  //printf("    *** %s ***\n",A->ArgLine());
  // First argument is the command
  if (A->Command()==NULL) {
    if (debug>0) mprintf("NULL Command.\n");
    return;
  }

  // Check if command pertains to coordinate lists
  // If it does, get a parm based on parm/parmindex keywords in arg list
  if (A->CommandIs("trajin")) {
    tempParm = parmFileList.GetParm(A);
    trajinList.Add(NULL, A, tempParm);
    return;
  }
  if (A->CommandIs("reference")) {
    tempParm = parmFileList.GetParm(A);
    referenceList.Add(NULL, A, tempParm);
    return;
  }
  if (A->CommandIs("trajout")) {
    tempParm = parmFileList.GetParm(A);
    trajoutList.Add(NULL, A, tempParm);
    return;
  }

  if (A->CommandIs("parm")) {
    parmFileList.Add(A->getNextString());
    return;
  }

  if (A->CommandIs("noprogress")) {
    showProgress=false;
    mprintf("    noprogress: Progress bar will not be shown.\n");
    return;
  }

  if (A->CommandIs("debug")) {
    SetGlobalDebug( A->getNextInteger(0) );
    return ;
  }

  if (A->CommandIs("parminfo")) {
    if ( (tempParm=parmFileList.GetParm(A->getNextInteger(0)))!=NULL ) {
      tempMask = new AtomMask();
      tempMask->SetMaskString( A->getNextMask() );
      tempMask->SetupCharMask( tempParm, debug);
      for (int atom=0; atom < tempParm->natom; atom++) {
        if (tempMask->AtomInCharMask(atom)) tempParm->AtomInfo(atom);
      }
      delete tempMask;
    }
    return;
  }

  // Check if command pertains to datafiles
  if ( A->CommandIs("datafile") ) {
    DF_Args.push_back(A->Copy());
    return;
  }

  // Check if command pertains to an action
  if ( ptrajActionList.Add(A)==0 ) return; 

  mprintf("Warning: Unknown Command %s.\n",A->Command());
}

/*
 * PtrajState::ProcessDataFileCmd()
 * Process command relating to data files. All datafile commands have format:
 *   datafile <cmd> <datafile> ...
 */
void PtrajState::ProcessDataFileCmd() {
  char *df_cmd = NULL;
  char *name1 = NULL;
  char *name2 = NULL;
  int width,precision;
  DataFile *df;

  if (DF_Args.empty()) return;
  mprintf("DATAFILE SETUP:\n");

  // Loop through all "datafile" arguments
  for (std::list<ArgList*>::iterator it=DF_Args.begin(); it!=DF_Args.end(); it++) {
    A = (*it);
    // Next string will be the argument passed to datafile
    df_cmd = A->getNextString();
    mprintf("  [%s]\n",A->ArgLine());
    // Next string is datafile that command pertains to
    name1 = A->getNextString();
    if (name1==NULL) {
      mprintf("Error: datafile %s: No filename given.\n",df_cmd);
      continue;
    }
    df = DFL.GetDataFile(name1);

    // datafile create
    // Usage: datafile create <filename> <dataset0> <dataset1> ...
    if ( strcmp(df_cmd,"create")==0 ) {
      if (df==NULL)
        mprintf("    Creating file %s\n",name1);
      while ( (name2=A->getNextString())!=NULL ) {
        if ( DFL.Add(name1, DSL.Get(name2))==NULL ) {
          mprintf("Warning: Dataset %s does not exist in main dataset list, skipping.\n",name2);
        }
      }

    // datafile xlabel
    // Usage: datafile xlabel <filename> <label>
    } else if ( strcmp(df_cmd,"xlabel")==0 ) {
      if (df==NULL) {
        mprintf("Error: datafile xlabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetXlabel(A->getNextString());

    // datafile invert
    // Usage: datafile invert <filename>
    } else if ( strcmp(df_cmd,"invert")==0 ) {
      if (df==NULL) {
        mprintf("Error: datafile invert: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Inverting datafile %s\n",name1);
      df->SetInverted();

    // datafile noxcol
    // Usage: datafile noxcol <filename>
    } else if ( strcmp(df_cmd,"noxcol")==0 ) {
      if (df==NULL) {
        mprintf("Error: datafile noxcol: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Not printing x column for datafile %s\n",name1);
      df->SetNoXcol();
    
    // datafile precision
    // Usage: datafile precision <filename> <dataset> [<width>] [<precision>]
    //        If width/precision not specified default to 12.4
    } else if ( strcmp(df_cmd,"precision")==0 ) {
      if (df==NULL) {
        mprintf("Error: datafile precision: DataFile %s does not exist.\n",name1);
        continue;
      }
      name2 = A->getNextString();
      width = A->getNextInteger(12);
      precision = A->getNextInteger(4);
      df->SetPrecision(name2,width,precision);
    }

  } // END loop over datafile args
}  

/* 
 * PtrajState::Run()
 * Process trajectories in trajFileList. Each frame in trajFileList is sent
 * to the actions in ptrajActionList for processing.
 */
int PtrajState::Run() {
  std::list<TrajectoryFile*>::iterator traj;
  int maxFrames=0;        // Total # of frames that will be read
  int actionSet=0;        // Internal data frame
  int readSets=0;         // Number of frames actually read
  int lastPindex=-1;      // Index of the last loaded parm file
  AmberParm *CurrentParm=NULL; 
  Frame *CurrentFrame=NULL;
  FrameList refFrames;

  // ========== S E T U P   P H A S E ========== 
  // Calculate frame division among trajectories
  mprintf("\nINPUT TRAJECTORIES:\n");
  maxFrames=trajinList.SetupFrames();
  //trajinList.Info(1);
  if (maxFrames==-1)
    mprintf("  Coordinate processing will occur until EOF (unknown number of frames).\n");
  else
    mprintf("  Coordinate processing will occur on %i frames.\n",maxFrames);

  // Parameter file information
  parmFileList.Print();

  // Setup reference frames if reference files were specified 
  referenceList.SetupRefFrames(&refFrames);
  refFrames.Info();

  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.Info(0);
 
  // Set max frames in the data set list
  DSL.SetMax(maxFrames); 
  
  // Initialize actions and set up data set and data file list
  ptrajActionList.Init( &DSL, &refFrames, &DFL, &parmFileList);

  // ========== R U N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN TRAJECTORY PROCESSING:\n");
  for (traj=trajinList.begin(); traj!=trajinList.end(); traj++) {
    // Open up the trajectory file. If an error occurs, bail 
    if ((*traj)->BeginTraj(showProgress) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->TrajName());
      break;
    }
    // Set current parm from current traj.
    CurrentParm = (*traj)->TrajParm();

    // If Parm has changed set up the Action list for new topology file
    if (lastPindex != CurrentParm->pindex) {
      if (ptrajActionList.Setup( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->parmName);
        continue;
      }
      //fprintf(stdout,"DEBUG: After setup of Actions in PtrajState parm name is %s\n",
      //        CurrentParm->parmName);
      // Set up the Frame for this parm
      if (CurrentFrame!=NULL) delete CurrentFrame;
      CurrentFrame = new Frame(CurrentParm->natom, CurrentParm->mass);
      lastPindex = CurrentParm->pindex;
    }

    // Loop over every Frame in trajectory
    while ( (*traj)->GetNextFrame(CurrentFrame->X, CurrentFrame->box, &(CurrentFrame->T)) ) {
      // Perform Actions on Frame
      ptrajActionList.DoActions(&CurrentFrame, actionSet);
      // Do Output
      trajoutList.Write(actionSet, CurrentParm, CurrentFrame);
#ifdef DEBUG
      dbgprintf("\tDEBUG: %30s: %4i\n",CurrentParm->parmName,CurrentParm->outFrame);
#endif
      // Increment frame counters
      actionSet++; 
    }

    // Close the trajectory file
    (*traj)->EndTraj();
    // Update how many frames across all threads have been written for parm
    // Do this for the original parm since it may have been modified.
    (*traj)->TrajParm()->outFrame += (*traj)->Total_Read_Frames();
    readSets+=(*traj)->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  if (CurrentFrame!=NULL) delete CurrentFrame;
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output traj
  trajoutList.Close();

  // Do action output
  ptrajActionList.Print( );

  // Print Dataset information
  DSL.Info();
  // Process any datafile commands
  ProcessDataFileCmd();
  // Print Datafile information
  DFL.Info();

  // Do dataset output - first sync datasets
  DSL.Sync();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();
 
  return 0;
}
