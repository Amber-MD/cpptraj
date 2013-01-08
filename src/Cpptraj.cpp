#include <cstdio> 
#include <cstdlib> // system
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "CpptrajStdio.h"
#include "Trajin_Multi.h"
#include "FrameArray.h"
#include "ReadLine.h"
#include "ParmFile.h"
#include "DataSet_Coords.h" // CrdAction

void Cpptraj::Usage(const char* programName) {
  mprinterr("\nUsage: %s [-p <Top1>, -p <Top2>, ...] [-i <Input1>, -i <Input2>, ...]\n",
            programName);
  mprinterr(  "       %s <Top> <Input>\n",programName);
  mprinterr(  "       Additional options:\n");
  mprinterr(  "         --help, -help : Print usage information and exit.\n");
  mprinterr(  "         -V, --version : Print version information and exit.\n");
  mprinterr(  "         --defines     : Print list of defines used in compilation.\n");
  mprinterr(  "         -debug <N>    : Set global debug level.\n");
  mprinterr(  "         --interactive : Enter interactive mode.\n");
  mprinterr(  "         --log <file>  : Record commands used in interactive mode to <file>.\n");
}

void Cpptraj::Help_Help() {
  mprintf("help [<cmd>]\n");
}

static const char TypeList[] = 
  "(<type> = actions,trajin,trajout,ref,parm,analysis,datafile,dataset)";

void Cpptraj::Help_List() {
  mprintf("list <type> %s\n", TypeList);
}

void Cpptraj::Help_Debug() {
  mprintf("debug [<type>] <#> %s\n", TypeList);
}

void Cpptraj::Help_Clear() {
  mprintf("clear [ {all | <type>} ] %s\n", TypeList);
  mprintf("\tAll lists will be cleared only if 'all' is specified by itself.\n");
}

void Cpptraj::Help_ActiveRef() {
  mprintf("activeref <#>\n");
  mprintf("\tSet the reference structure to be used for coordinate-based mask parsing.\n");
}

void Cpptraj::Help_Create_DataFile() {
  mprintf("create <filename> <dataset0> <dataset1> ...\n");
}

void Cpptraj::Help_Precision() {
  mprintf("precision {<filename> | <dataset arg>} [<width>] [<precision>]\n");
  mprintf("\tSet precision for all datasets in datafile <filename> or\n");
  mprintf("dataset(s) specified by <dataset arg> to <width>.<precision>\n");
  mprintf("If width/precision not specified default to 12.4\n");
}

void Cpptraj::Help_SelectDS() {
  mprintf("selectds <dataset selection>\n");
  mprintf("\tShow results of data set selection. Data set selection format is:\n");
  mprintf("\t\t<name>[<aspect]:<idx range>\n");
  mprintf("Where '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n");
  mprintf("and <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n");
  mprintf("The aspect and index portions may be optional. An asterisk '*' may be used as\n");
  mprintf("a wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

void Cpptraj::Help_Trajin() {
  mprintf("trajin <filename> {[<start>] [<stop> | last] [offset]} | lastframe\n");
  mprintf("       [parm <parmfile> | parmindex <#>]\n");
  mprintf("       [ remdtraj [remdtrajtemp <T> | remdtrajidx <#>]\n");
  mprintf("         [trajnames <rep1>,<rep2>,...,<repN> ] ]\n");
}

void Cpptraj::Help_Ensemble() {
  mprintf("ensemble <file0> {[<start>] [<stop> | last] [offset]} | lastframe\n");
  mprintf("          [parm <parmfile> | parmindex <#>]\n");
  mprintf("         [trajnames <file1>,<file2>,...,<fileN>\n");
}

void Cpptraj::Help_Trajout() {
  mprintf("trajout <filename> [<fileformat>] [append] [nobox]\n");
  mprintf("        [parm <parmfile> | parmindex <#>] [onlyframes <range>] ]title <title>]\n");
  mprintf("        <Format Options>\n");
}

void Cpptraj::Help_Reference() {
  mprintf("reference <filename> [<frame#>] [<mask>] [TAG] [lastframe]\n");
  mprintf("          [average [<stop>] [<offset>]]\n");
}

void Cpptraj::Help_Parm() {
  mprintf("parm <filename> [<tag>] [nobondsearch | bondsearch [<offset>]]\n");
  mprintf("\tAdd <filename> to parm list\n");
}

void Cpptraj::Help_ParmInfo() {
  mprintf("parminfo [<parmindex>] [<mask>]:\n");
  mprintf("\tPrint information on parm <parmindex> (0 by default). If <mask> is given\n");
  mprintf("print info on atoms in mask. If no mask given print overall information.\n");
}

void Cpptraj::Help_ParmWrite() {
  mprintf("parmwrite out <filename> [<parmindex>]\n");
  mprintf("\tWrite parm <parmindex> to <filename>\n");
}

void Cpptraj::Help_ParmStrip() {
  mprintf("parmstrip <mask> [<parmindex>]\n");
  mprintf("\tStrip atoms in mask from parm\n");
}

void Cpptraj::Help_ParmBox() {
  mprintf("[parm]box [<parmindex>] [x <xval>] [y <yval>] [z <zval>]");
  mprintf(" [alpha <a>] [beta <b>] [gamma <g>] [nobox]\n");
  mprintf("\tSet the given parm box info to what is specified. If nobox, remove box info.\n");
}

void Cpptraj::Help_Solvent() {
  mprintf("solvent [<parmindex>] <mask>\n");
  mprintf("\tSet solvent for the given parm (default 0) based on <mask>\n");
}

void Cpptraj::Help_BondInfo() {
  mprintf("bondinfo [<parmindex>] <mask>\n");
  mprintf("\tPrint bond information for parm <parmindex> (0 by default).\n");
}

void Cpptraj::Help_ResInfo() {
  mprintf("resinfo [<parmindex>]\n");
  mprintf("\tPrint residue information for parm <parmindex> (0 by default).\n");
}

void Cpptraj::Help_MolInfo() {
  mprintf("molinfo [<parmindex>] <mask>\n");
  mprintf("\tPrint molecule information for parm <parmindex> (0 by default).\n");
}

void Cpptraj::Help_CrdAction() {
  mprintf("crdaction <crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n");
}

void Cpptraj::Help_CrdOut() {
  mprintf("crdout <crd set> <filename> [<trajout args>] [crdframes <start>,<stop>,<offset>]\n");
}

void Cpptraj::Help_RunAnalysis() {
  mprintf("runanalysis [<analysis> [<analysis args>]]\n");
  mprintf("\tIf specified alone, run all analyses in the analysis list.\n");
  mprintf("Otherwise run the specified analysis immediately.\n");
}

// -----------------------------------------------------------------------------
enum GeneralCmdTypes { LIST = 0, HELP, QUIT, RUN, DEBUG, NOPROG, NOEXITERR, 
                       SYSTEM, ACTIVEREF, READDATA, CREATE, PRECISION, DATAFILE,
                       SELECT, SELECTDS, READINPUT, RUN_ANALYSIS, WRITEDATA,
                       CLEAR, CRDACTION, CRDOUT };

const DispatchObject::Token Cpptraj::GeneralCmds[] = {
  { DispatchObject::GENERAL, "activeref",     0, Help_ActiveRef,       ACTIVEREF    },
  { DispatchObject::GENERAL, "clear",         0, Help_Clear,           CLEAR        },
  { DispatchObject::GENERAL, "crdaction",     0, Help_CrdAction,       CRDACTION    },
  { DispatchObject::GENERAL, "crdout",        0, Help_CrdOut,          CRDOUT       },
  { DispatchObject::GENERAL, "create",        0, Help_Create_DataFile, CREATE       },
  { DispatchObject::GENERAL, "datafile",      0, 0,                    DATAFILE     },
  { DispatchObject::GENERAL, "debug",         0, Help_Debug,           DEBUG        },
  { DispatchObject::GENERAL, "exit" ,         0, 0,                    QUIT         },
  { DispatchObject::GENERAL, "gnuplot" ,      0, 0,                    SYSTEM       },
  { DispatchObject::GENERAL, "go"   ,         0, 0,                    RUN          },
  { DispatchObject::GENERAL, "head" ,         0, 0,                    SYSTEM       },
  { DispatchObject::GENERAL, "help" ,         0, Help_Help,            HELP         },
  { DispatchObject::GENERAL, "list" ,         0, Help_List,            LIST         },
  { DispatchObject::GENERAL, "ls",            0, 0,                    SYSTEM       },
  { DispatchObject::GENERAL, "noexitonerror", 0, 0,                    NOEXITERR    },
  { DispatchObject::GENERAL, "noprogress",    0, 0,                    NOPROG       },
  { DispatchObject::GENERAL, "precision",     0, Help_Precision,       PRECISION    },
  { DispatchObject::GENERAL, "prnlev",        0, Help_Debug,           DEBUG        },
  { DispatchObject::GENERAL, "pwd",           0, 0,                    SYSTEM       },
  { DispatchObject::GENERAL, "quit" ,         0, 0,                    QUIT         },
  { DispatchObject::GENERAL, "readdata",      0, 0,                    READDATA     },
  { DispatchObject::GENERAL, "readinput",     0, 0,                    READINPUT    },
  { DispatchObject::GENERAL, "run"   ,        0, 0,                    RUN          },
  { DispatchObject::GENERAL, "runanalysis",   0, Help_RunAnalysis,     RUN_ANALYSIS },
  { DispatchObject::GENERAL, "select",        0, 0,                    SELECT    },
  { DispatchObject::GENERAL, "selectds",      0, Help_SelectDS,        SELECTDS     },
  { DispatchObject::GENERAL, "writedata",     0, 0,                    WRITEDATA    },
  { DispatchObject::GENERAL, "xmgrace",       0, 0,                    SYSTEM       },
  { DispatchObject::NONE,    0,               0, 0,                    0            }
};

enum CoordCmdTypes { REFERENCE, TRAJIN, TRAJOUT };

const DispatchObject::Token Cpptraj::CoordCmds[] = {
  { DispatchObject::COORD, "reference", 0, Help_Reference, REFERENCE },
  { DispatchObject::COORD, "trajin",    0, Help_Trajin,    TRAJIN    },
  { DispatchObject::COORD, "ensemble",  0, Help_Ensemble,  TRAJIN    },
  { DispatchObject::COORD, "trajout",   0, Help_Trajout,   TRAJOUT   },
  { DispatchObject::NONE,  0,           0, 0,              0         }
};

enum ParmCmdTypes { LOADPARM=0, PARMINFO, PARMWRITE, PARMSTRIP, PARMBOX,
                    SOLVENT, BONDINFO, RESINFO, MOLINFO };

const DispatchObject::Token Cpptraj::ParmCmds[] = {
  { DispatchObject::PARM, "bondinfo",     0, Help_BondInfo,  BONDINFO  },
  { DispatchObject::PARM, "molinfo",      0, Help_MolInfo,   MOLINFO   },
  { DispatchObject::PARM, "parm",         0, Help_Parm,      LOADPARM  },
  { DispatchObject::PARM, "parmbondinfo", 0, Help_BondInfo,  BONDINFO  },
  { DispatchObject::PARM, "parmbox",      0, Help_ParmBox,   PARMBOX   },
  { DispatchObject::PARM, "parminfo",     0, Help_ParmInfo,  PARMINFO  },
  { DispatchObject::PARM, "parmmolinfo",  0, Help_MolInfo,   MOLINFO   },
  { DispatchObject::PARM, "parmresinfo",  0, Help_ResInfo,   RESINFO   },
  { DispatchObject::PARM, "parmstrip",    0, Help_ParmStrip, PARMSTRIP },
  { DispatchObject::PARM, "parmwrite",    0, Help_ParmWrite, PARMWRITE },
  { DispatchObject::PARM, "resinfo",      0, Help_ResInfo,   RESINFO   },
  { DispatchObject::PARM, "solvent",      0, Help_Solvent,   SOLVENT   },
  { DispatchObject::NONE, 0,              0, 0,              0         }
};

const DispatchObject::Token Cpptraj::Deprecated[] = {
  { DispatchObject::DEPRECATED, "molsearch",    0, 0, 0 },
  { DispatchObject::DEPRECATED, "nomolsearch",  0, 0, 0 },
  { DispatchObject::DEPRECATED, "bondsearch",   0, 0, 0 },
  { DispatchObject::DEPRECATED, "nobondsearch", 0, 0, 0 },
  { DispatchObject::NONE      , 0,              0, 0, 0 }
};

// -----------------------------------------------------------------------------
// Constructor
Cpptraj::Cpptraj() : 
  debug_(0),
  showProgress_(true),
  exitOnError_(true),
  nrun_(0)
{}

/** List all commands, or call help function of specific command. */
void Cpptraj::Help(ArgList& argIn) {
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty()) {
    mprintf("General Commands:\n");
    ListAllCommands( GeneralCmds );
    mprintf("Topology Commands:\n");
    ListAllCommands( ParmCmds );
    mprintf("Coordinate Commands:\n");
    ListAllCommands( CoordCmds );
    mprintf("Action Commands:\n");
    ListAllCommands( ActionList::DispatchArray );
    mprintf("Analysis Commands:\n");
    ListAllCommands( AnalysisList::DispatchArray );
  } else {
    DispatchObject::TokenPtr dispatchToken = SearchToken( arg );
    if (dispatchToken == 0 || dispatchToken->Help == 0) 
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken->Help();
  }
}

/// Types of lists
enum ListType { L_ACTION = 0, L_TRAJIN, L_REF, L_TRAJOUT, L_PARM, L_ANALYSIS,
                L_DATAFILE, L_DATASET, N_LISTS };
/// Select lists from ArgList
static std::vector<bool> ListsFromArg( ArgList& argIn, bool allowEnableAll ) {
  std::vector<bool> enabled( (int)N_LISTS );
  enabled[L_ACTION] = argIn.hasKey("actions");
  enabled[L_TRAJIN] = argIn.hasKey("trajin");
  enabled[L_REF] = argIn.hasKey("ref");
  enabled[L_TRAJOUT] = argIn.hasKey("trajout");
  enabled[L_PARM] = argIn.hasKey("parm");
  enabled[L_ANALYSIS] = argIn.hasKey("analysis");
  enabled[L_DATAFILE] = argIn.hasKey("datafile");
  enabled[L_DATASET] = argIn.hasKey("dataset");
  if (!allowEnableAll) return enabled;
  // If nothing is enabled, set all enabled
  bool nothing_enabled = true;
  for (std::vector<bool>::iterator en = enabled.begin(); en != enabled.end(); ++en)
    if (*en) {
      nothing_enabled = false;
      break;
    }
  if (nothing_enabled) enabled.assign( (int)N_LISTS, true );
  return enabled;
}

/** List all members of specified lists. */
void Cpptraj::List(ArgList& argIn) {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  if ( enabled[L_ACTION] ) actionList.List();
  if ( enabled[L_TRAJIN] ) trajinList.List();
  if ( enabled[L_REF] ) refFrames.List();
  if ( enabled[L_TRAJOUT] ) trajoutList.List();
  if ( enabled[L_PARM] ) parmFileList.List();
  if ( enabled[L_ANALYSIS] ) analysisList.List();
  if ( enabled[L_DATAFILE] ) DFL.List();
  if ( enabled[L_DATASET] ) DSL.List();
}

/** Set debug level of specified lists */
void Cpptraj::Debug(ArgList& argIn) {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  debug_ = argIn.getNextInteger(0);
  if ( enabled[L_ACTION] ) actionList.SetDebug( debug_ );
  if ( enabled[L_TRAJIN] ) trajinList.SetDebug( debug_ );
  if ( enabled[L_REF] ) refFrames.SetDebug( debug_ );
  if ( enabled[L_TRAJOUT] ) trajoutList.SetDebug( debug_ );
  if ( enabled[L_PARM] ) parmFileList.SetDebug( debug_ );
  if ( enabled[L_ANALYSIS] ) analysisList.SetDebug( debug_ );
  if ( enabled[L_DATAFILE] ) DFL.SetDebug( debug_ );
  if ( enabled[L_DATASET] ) DSL.SetDebug( debug_ );
}

/** Clear specified lists */
void Cpptraj::Clear(ArgList& argIn) {
  std::vector<bool> enabled = ListsFromArg( argIn, argIn.hasKey("all") );
  if ( enabled[L_ACTION] ) actionList.Clear();
  if ( enabled[L_TRAJIN] ) trajinList.Clear();
  if ( enabled[L_REF] ) refFrames.Clear();
  if ( enabled[L_TRAJOUT] ) trajoutList.Clear();
  if ( enabled[L_PARM] ) parmFileList.Clear();
  if ( enabled[L_ANALYSIS] ) analysisList.Clear();
  if ( enabled[L_DATAFILE] ) DFL.Clear();
  if ( enabled[L_DATASET] ) DSL.Clear();
}

/** Add DataFile to DataFileList using specified sets. */
int Cpptraj::Create_DataFile(ArgList& dataArg) {
  // Next string is datafile that command pertains to.
  std::string name1 = dataArg.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: create: No filename given.\nError: Usage: ");
    Help_Create_DataFile();
    return 1;
  }
  DataFile* df = DFL.AddDataFile(name1, dataArg);
  if (df==0) {
    mprinterr("Error: create: Could not create file %s:",name1.c_str());
    return 1;
  }
  // Process any recognized datafile args
  df->ProcessArgs( dataArg );
  // Treat all remaining args as dataset names
  int err = 0;
  ArgList dsetArgs = dataArg.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList Sets = DSL.GetMultipleSets( *dsa );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", (*dsa).c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      mprintf(" %s", (*set)->Legend().c_str());
      if ( df->AddSet(*set) ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->Legend().c_str());
        ++err;
      }
    }
  }
  mprintf("\n");
  return err;
}

/** Set precision for specific set or all sets in specified DataFile */
int Cpptraj::Precision(ArgList& dataArg) {
  // Next string is DataSet(s)/DataFile that command pertains to.
  std::string name1 = dataArg.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: precision: No filename/setname given.\nError: Usage: ");
    Help_Precision();
    return 1;
  }
  // This will break if dataset name starts with a digit...
  int width = dataArg.getNextInteger(12);
  if (width < 1) {
    mprintf("Error: precision: Cannot set width < 1 (%i).\n", width);
    return 1;
  }
  int precision = dataArg.getNextInteger(4);
  if (precision < 0) precision = 0;
  DataFile* df = DFL.GetDataFile(name1);
  if (df != 0) {
    mprintf("\tSetting precision for all sets in %s to %i.%i\n", df->Filename(),
            width, precision);
    df->SetPrecision(width, precision);
  } else {
    DataSetList dsets = DSL.GetMultipleSets( name1 );
    mprintf("\tSetting precision for %i sets to %i.%i\n", dsets.size(),
            width, precision);
    for (DataSetList::const_iterator set = dsets.begin(); set != dsets.end(); ++set)
      (*set)->SetPrecision(width, precision);
  }
  return 0;
}

/** Read data from file into master DataSetList */
int Cpptraj::ReadData(ArgList& argIn) {
  DataFile dataIn;
  if (dataIn.ReadData( argIn, DSL )!=0) {
    mprinterr("Error: Could not read data file.\n");
    return 1;
  }
  return 0;
}

/** Show results of DataSet selection */
void Cpptraj::SelectDS(ArgList& argIn) {
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = DSL.GetMultipleSets( dsarg );
  mprintf("SelectDS: Arg [%s]:", dsarg.c_str());
  dsets.List();
}

// -----------------------------------------------------------------------------
/** Load file into TopologyList */
int Cpptraj::LoadParm(ArgList& argIn) {
  std::string parmtag = argIn.getNextTag();
  bool bondsearch = !argIn.hasKey("nobondsearch");
  double offset = argIn.getKeyDouble("bondsearch", -1.0);
  return parmFileList.AddParmFile(argIn.GetStringNext(), parmtag, bondsearch, offset);
}

/** Print information for specified parm */
int Cpptraj::ParmInfo(ArgList& argIn, int cmdidxIn) {
  ParmCmdTypes cmdidx = (ParmCmdTypes) cmdidxIn;
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parm index %i not loaded.\n",pindex);
    return 1;
  }
  std::string maskarg = argIn.GetMaskNext();
  switch (cmdidx) {
    case PARMINFO:
      if (!maskarg.empty())
        parm->PrintAtomInfo( maskarg );
      else
        parm->Summary();
      break;
    case BONDINFO: parm->PrintBondInfo( maskarg ); break;
    case RESINFO : parm->PrintResidueInfo( maskarg ); break;
    case MOLINFO : parm->PrintMoleculeInfo( maskarg ); break;
    default: return 1; // Should never get here
  }
  return 0;
}

/** Write parm to Amber Topology file. */
int Cpptraj::ParmWrite(ArgList& argIn) {
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: parmwrite: No output filename specified (use 'out <filename>').\n");
    return 1;
  }
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parmwrite: parm index %i not loaded.\n",pindex);
    return 1;
  }
  mprintf("\tWriting parm %i (%s) to Amber parm %s\n",pindex,
          parm->c_str(), outfilename.c_str());
  ParmFile pfile;
  pfile.Write( *parm, outfilename, ParmFile::AMBERPARM, debug_ );
  return 0;
}

// Cpptraj::ParmStrip()
int Cpptraj::ParmStrip(ArgList& argIn) {
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parmstrip: parm index %i not loaded.\n",pindex);
    return 1;
  }
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMask();
  parm->SetupIntegerMask( tempMask );
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           parm->Natom() - tempMask.Nselected(), parm->c_str());
  Topology* tempParm = parm->modifyStateByMask(tempMask);
  if (tempParm==0) {
    mprinterr("Error: parmstrip: Could not strip parm.\n");
    return 1;
  } else {
    // Replace parm with stripped version
    tempParm->ParmInfo();
    mprintf("\n");
    parmFileList.ReplaceParm(pindex, tempParm);
  }
  return 0;
}

/** Modify parm box information. */
int Cpptraj::ParmBox(ArgList& argIn) {
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parmbox: parm index %i not loaded.\n",pindex);
    return 1;
  }
  if ( argIn.hasKey("nobox") )
    parm->SetBox( Box() );
  else {
    Box pbox;
    pbox.SetX( argIn.getKeyDouble("x",0) );
    pbox.SetY( argIn.getKeyDouble("y",0) );
    pbox.SetZ( argIn.getKeyDouble("z",0) );
    pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
    pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
    pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( parm->ParmBox() );
    parm->SetBox( pbox );
  }
  return 0;
}

/** Modify parm solvent information */
int Cpptraj::ParmSolvent(ArgList& argIn) {
  std::string maskexpr = argIn.GetMaskNext();
  if ( maskexpr.empty() ) {
    mprinterr("Error: solvent: No mask specified.\n");
    return 1;
  }
  // Get parm index
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: solvent: parm index %i not loaded.\n",pindex);
    return 1;
  } 
  parm->SetSolvent( maskexpr );
  return 0;
}

/** Show results of mask expression */
int Cpptraj::Select(ArgList& argIn) {
  AtomMask tempMask( argIn.GetMaskNext() );
  int pindex = argIn.getNextInteger(0);
  Topology* parm = parmFileList.GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: solvent: parm index %i not loaded.\n",pindex);
    return 1;
  }
  parm->SetupIntegerMask( tempMask );
  mprintf("Selected %i atoms.\n", tempMask.Nselected());
  if (!argIn.hasKey("total"))
    tempMask.PrintMaskAtoms("Selected");
  return 0;
}

// -----------------------------------------------------------------------------
// Cpptraj::CrdAction()
/** Perform action on given COORDS dataset */
int Cpptraj::CrdAction(ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdaction: Specify COORDS dataset name.\n");
    Help_CrdAction();
    return 1;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)DSL.FindSetOfType( setname, DataSet::COORDS );
  if (CRD == 0) {
    mprinterr("Error: crdaction: No COORDS set with name %s found.\n", setname.c_str());
    return 1;
  }
  // Start, stop, offset
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  int start = crdarg.getNextInteger(1) - 1;
  int stop = crdarg.getNextInteger(CRD->Size());
  int offset = crdarg.getNextInteger(1);
  if (debug_ > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  DispatchObject::TokenPtr tkn = SearchTokenArray( ActionList::DispatchArray, actionargs);
  if ( tkn == 0 ) return 1;
  Action* act = (Action*)tkn->Alloc();
  if (act == 0) return 1;
  if ( act->Init( actionargs, &parmFileList, &refFrames, &DSL, &DFL, debug_ ) != Action::OK ) {
    delete act;
    return 1;
  }
  actionargs.CheckForMoreArgs();
  // Set up frame and parm for COORDS.
  Topology* originalParm = new Topology();
  *originalParm = CRD->Top();
  Frame* originalFrame = new Frame( CRD->Top().Atoms() );
  // Set up for this topology
  Topology* currentParm = originalParm;
  if ( act->Setup( currentParm, &currentParm ) == Action::ERR ) {
    delete act;
    return 1;
  }
  // Check if parm was modified. If so, update COORDS.
  if ( currentParm != originalParm ) {
    mprintf("Info: crdaction: Parm for %s was modified by action %s\n", 
            CRD->Legend().c_str(), actionargs.Command());
    CRD->SetTopology( *currentParm );
  }
  // Loop over all frames in COORDS.
  ProgressBar progress( stop );
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( frame );
    CRD->GetFrame( frame, *originalFrame );
    Frame* currentFrame = originalFrame;
    if (act->DoAction( frame, currentFrame, &currentFrame ) == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i\n", frame + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    // TODO: Have actions indicate whether they will modify coords
    //if ( currentFrame != originalFrame ) 
      CRD->SetCRD( frame, *currentFrame );
  }
  delete originalFrame;
  delete originalParm;
  act->Print();
  delete act;
  return 0;
}

// Cpptraj::CrdOut()
/** Write out COORDS dataset */
int Cpptraj::CrdOut(ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    Help_CrdOut();
    return 1;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)DSL.FindSetOfType( setname, DataSet::COORDS );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return 1;
  }
  setname = argIn.GetStringNext();
  // Start, stop, offset
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  int start = crdarg.getNextInteger(1) - 1;
  int stop = crdarg.getNextInteger(CRD->Size());
  int offset = crdarg.getNextInteger(1);
  if (debug_ > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  Trajout outtraj;
  Topology* currentParm = (Topology*)&(CRD->Top()); // TODO: Fix cast
  if (outtraj.SetupTrajWrite( setname, &argIn, currentParm, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: crdout: Could not set up output trajectory.\n");
    return 1;
  }
  outtraj.PrintInfo( 1 );
  Frame currentFrame( CRD->Top().Natom() );
  ProgressBar progress( stop );
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( frame );
    CRD->GetFrame( frame, currentFrame );
    if ( outtraj.WriteFrame( frame, currentParm, currentFrame ) ) {
      mprinterr("Error writing %s to output trajectory, frame %i.\n", 
                CRD->Legend().c_str(), frame + 1);
      break;
    }
  }
  return 0;
}

// Cpptraj::CrdAnalyze()
/** Run a single analysis. */
int Cpptraj::CrdAnalyze(ArgList& argIn) {
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  DispatchObject::TokenPtr tkn = SearchTokenArray( AnalysisList::DispatchArray, analyzeargs);
  if ( tkn == 0 ) return 1;
  Analysis* ana = (Analysis*)tkn->Alloc();
  if (ana == 0) return 1;
  if ( ana->Setup( analyzeargs, &DSL, &parmFileList, debug_ ) != Analysis::OK ) {
    delete ana;
    return 1;
  }
  if (ana->Analyze() == Analysis::OK)
    ana->Print(&DFL);
  delete ana;
  return 0;
}

// -----------------------------------------------------------------------------
// Cpptraj::ListAllCommands()
/** List all commands in the given token array. */
void Cpptraj::ListAllCommands(DispatchObject::TokenPtr DispatchArray) {
  int col = 0;
  mprintf("\t");
  for (DispatchObject::TokenPtr token = DispatchArray;
                                token->Type != DispatchObject::NONE; ++token)
  {
    mprintf("%s  ", token->Cmd);
    ++col;
    if (col == 8) {
      mprintf("\n\t");
      col = 0;
    }
  }
  mprintf("\n");
}

// Cpptraj::SearchTokenArray()
/** Search the given token array for command.
  * \return DispatchObject associated with command or 0 if not found.
  */
DispatchObject::TokenPtr Cpptraj::SearchTokenArray(DispatchObject::TokenPtr DispatchArray, 
                                                   ArgList const& arg)
{
  for (DispatchObject::TokenPtr token = DispatchArray;
                                token->Type != DispatchObject::NONE; ++token)
  {
    if ( arg.CommandIs( token->Cmd ) ) 
      return token;
  }
  return 0;
}

// Cpptraj::SearchToken()
/** Search each token list for the given command. If the command is found in
  * a list then dispatchToken is set by SearchTokenArray and 1 is returned.
  * \return the token if found, 0 if not.
  */
DispatchObject::TokenPtr Cpptraj::SearchToken(ArgList& argIn) {
  DispatchObject::TokenPtr tkn = 0;
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    return SearchTokenArray( AnalysisList::DispatchArray, argIn);
  } else {
    tkn = SearchTokenArray( GeneralCmds, argIn );
    if (tkn != 0) return tkn;
    tkn = SearchTokenArray( ParmCmds, argIn );
    if (tkn != 0) return tkn;
    tkn = SearchTokenArray( CoordCmds, argIn );
    if (tkn != 0) return tkn;
    tkn = SearchTokenArray( ActionList::DispatchArray, argIn);
    if (tkn != 0) return tkn;
    tkn = SearchTokenArray( AnalysisList::DispatchArray, argIn);
    if (tkn != 0) return tkn;
    tkn = SearchTokenArray( Deprecated, argIn);
    if (tkn != 0) return tkn;
  }
  mprinterr("[%s]: Command not found.\n",argIn.Command());
  return 0;
}

// -----------------------------------------------------------------------------
// Cpptraj::Interactive()
Cpptraj::Mode Cpptraj::Interactive() {
  ReadLine inputLine;
  // By default when interactive do not exit on errors
  exitOnError_ = false;
  // Open log file if name has been set
  if (!logfile_.FullFileName().empty())
    logfile_.OpenFile();
  Mode readLoop = C_OK;
  while ( readLoop == C_OK ) {
    if (inputLine.GetInput()) break; 
    if (!inputLine.empty()) {
      readLoop = Dispatch( inputLine.c_str() );
      if (logfile_.IsOpen())
        logfile_.Printf("%s\n", inputLine.c_str());
    }
  }
  logfile_.CloseFile();
  // If we broke out of loop because of EOF and Run has been called at least
  // once, indicate that we can safely quit.
  if (readLoop == C_OK && nrun_ > 0) return C_QUIT;
  return readLoop;
}

static inline bool EndChar(char ptr) {
  if (ptr=='\n' || ptr=='\r' || ptr=='\0' || ptr==EOF) return true;
  return false;
}

/** Read commands from an input file. '#' indicates the beginning of a
  * comment, backslash at the end of a line indicates continuation
  * (otherwise indicates 'literal').
  * \return 0 if successfully read, 1 on error.
  */
Cpptraj::Mode Cpptraj::ProcessInput(std::string const& inputFilename) {
  FILE *infile;
  if (inputFilename.empty()) return C_ERR;
  mprintf("INPUT: Reading Input from file %s\n",inputFilename.c_str());
  if ( (infile=fopen(inputFilename.c_str(),"r"))==0 ) {
    rprintf("Error: Could not open input file %s\n",inputFilename.c_str());
    return C_ERR;
  }
  // Read in each line of input. Newline or null terminates. \ continues line.
  std::string inputLine;
  unsigned int idx = 0;
  char lastchar = '0';
  char ptr = 0;
  Mode cmode = C_OK;
  while ( ptr != EOF ) {
    ptr = (char)fgetc(infile);
    // Skip leading whitespace
    if (idx == 0 && isspace(ptr)) {
      while ( (ptr = (char)fgetc(infile))!=EOF )
        if ( !isspace(ptr) ) break;
    }
    // If '#' is encountered, skip the rest of the line
    if (ptr=='#') 
      while (!EndChar(ptr)) ptr=(char)fgetc(infile);
    // newline, null, or EOF terminates the line
    if (EndChar(ptr)) {
      // If no chars in string continue
      if (inputLine.empty()) continue;
      // Print the input line that will be sent to dispatch
      mprintf("  [%s]\n",inputLine.c_str());
      // Call Dispatch to convert input to arglist and process.
      cmode = Dispatch(inputLine.c_str());
      if (cmode != C_OK) break;
      // Reset Input line
      inputLine.clear();
      idx = 0;
      continue;
    }
    // Any consecutive whitespace is skipped
    if (idx > 0) lastchar = inputLine[idx-1];
    if (isspace(ptr) && isspace(lastchar)) continue;
    // Backslash followed by newline continues to next line. Otherwise backslash
    // followed by next char will be inserted. 
    if (ptr=='\\') {
      ptr = (char)fgetc(infile);
      if ( ptr == EOF ) break;
      if (ptr == '\n' || ptr == '\r') continue;
      inputLine += "\\";
      inputLine += ptr;
      idx += 2;
      continue;
    }
    // Add character to input line
    inputLine += ptr;
    ++idx;
  }
  fclose(infile);
  return cmode;
} 

/** Read command line args. */
Cpptraj::Mode Cpptraj::ProcessCmdLineArgs(int argc, char** argv) {
  if (argc == 1) return C_INTERACTIVE;
  bool hasInput = false;
  bool interactive = false;
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]); 
    if ( arg == "--help" || arg == "-help" ) {
      // --help, -help: Print usage and exit
      Usage( argv[0] );
      return C_QUIT;
    }
    if ( arg == "-V" || arg == "--version" ) 
      // -V, --version: Print version number and exit
      // Since version number should be printed before this is called, quit.
      return C_QUIT;
    if ( arg == "--defines" ) {
      // --defines: Print information on compiler defines used and exit
      mprintf("\nCompiled with:");
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
      return C_QUIT;
    }
    if ( arg == "--interactive" )
      interactive = true;
    else if ( arg == "-debug" && i+1 != argc) { 
      // -debug: Set overall debug level
      ArgList dbgarg( argv[++i] );
      Debug( dbgarg );
    } else if ( arg == "--log" && i+1 != argc)
      // --log: Set up log file for interactive mode
      logfile_.SetupWrite( argv[++i], debug_ );
    else if ( arg == "-p" && i+1 != argc) {
      // -p: Topology file
      if (parmFileList.AddParmFile( argv[++i] )) return C_ERR;
    } else if (arg == "-i" && i+1 != argc) {
      // -i: Input file(s)
      Cpptraj::Mode cmode = ProcessInput( argv[++i] );
      if (cmode == C_ERR) return C_ERR;
      if (cmode == C_QUIT) return C_QUIT;
      hasInput = true;
    } else if (arg == "-ms" && i+1 != argc) {
      // -ms: Mask string
      ArgList maskArg( argv[++i] );
      ParmInfo( maskArg, PARMINFO );
      return C_QUIT; 
    } else if ( i == 1 ) {
      // For backwards compatibility with PTRAJ; Position 1 = TOP file
      if (parmFileList.AddParmFile( argv[i])) return C_ERR;
    } else if ( i == 2 ) {
      // For backwards compatibility with PTRAJ; Position 2 = INPUT file
      Cpptraj::Mode cmode = ProcessInput( argv[i] );
      if (cmode == C_ERR) return C_ERR;
      if (cmode == C_QUIT) return C_QUIT;
      hasInput = true;
    } else {
      // Unrecognized
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      Usage(argv[0]);
      return C_QUIT;
    }
  }
  if (!hasInput || interactive) return C_INTERACTIVE;
  // If Run has already been called, just quit.
  if (nrun_ > 0) return C_QUIT;
  return C_OK;
}

// Cpptraj::Dispatch()
/** The input line is converted into a whitespace-delimited array of
  * arguments, the first of which is considered the command. This command
  * is searched for and if it is recognized it is sent to the appropriate
  * class. 
  * \param inputLine null-terminated string consisting of command and arguments.
  * \return C_OK if command was accepted or no error occurred.
  * \return C_ERR if error occurred.
  * \return C_QUIT if quit requested.
  */
Cpptraj::Mode Cpptraj::Dispatch(const char* inputLine) {
  int err = 0;
  //mprintf("\t[%s]\n", inputLine);
  ArgList command( inputLine );
  if ( command.empty() ) return C_OK;
  command.MarkArg(0); // Always mark the command
  DispatchObject::TokenPtr dispatchToken = SearchToken( command );
  if ( dispatchToken != 0 ) {
    //mprintf("TOKEN FOUND. CMD=%s  TYPE=%i\n", dispatchToken->Cmd, (int)dispatchToken->Type);
    switch (dispatchToken->Type) {
      case DispatchObject::PARM :
        switch ( dispatchToken->Idx ) {
          case LOADPARM : err = LoadParm( command ); break;
          case BONDINFO :
          case RESINFO  :
          case MOLINFO  :
          case PARMINFO : err = ParmInfo( command, dispatchToken->Idx ); break;
          case PARMWRITE: err = ParmWrite( command ); break;
          case PARMSTRIP: err = ParmStrip( command ); break;
          case PARMBOX  : err = ParmBox( command ); break;
          case SOLVENT  : err = ParmSolvent(command); break;
        } 
        break;
      case DispatchObject::COORD :
        switch ( dispatchToken->Idx ) {
          case TRAJIN :
            // Update # of sets to be read in for master DSL
            err = trajinList.AddTrajin(command, parmFileList); 
            if (err == 0) DSL.SetMax( trajinList.MaxFrames() );
            break;
          case TRAJOUT :
            // For setting up ensemble, save trajout arg
            trajoutArgs_.push_back(command);
            err = trajoutList.AddTrajout(command, parmFileList);
            break;
          case REFERENCE : err = refFrames.AddReference(command, parmFileList); break;
        }
        break;
      case DispatchObject::ACTION : 
        // For setting up ensemble, save action arg
        actionArgs_.push_back(command);
        err = actionList.AddAction( dispatchToken->Alloc, command, &parmFileList, 
                                    &refFrames, &DSL, &DFL );
        break;
      case DispatchObject::ANALYSIS :
        err = analysisList.AddAnalysis( dispatchToken->Alloc, command, &parmFileList, &DSL );
        break;
      case DispatchObject::GENERAL :
        switch ( dispatchToken->Idx ) {
          case HELP      : Help(command); break;
          case LIST      : List(command); break;
          case DEBUG     : Debug(command); break;
          case CLEAR     : Clear(command); break;
          case CRDACTION : err = CrdAction(command); break;
          case CRDOUT    : err = CrdOut(command); break;
          case SELECT    : err = Select(command); break; 
          case SELECTDS  : SelectDS(command); break;
          case NOPROG    : 
            showProgress_ = false; 
            mprintf("\tProgress bar will not be shown.\n");
            break;
          case NOEXITERR:
            exitOnError_ = false;
            mprintf("\tcpptraj will attempt to ignore errors if possible.\n");
            break;
          case ACTIVEREF : refFrames.SetActiveRef( command.getNextInteger(0) ); break;
          case READDATA  : err = ReadData( command ); break;
          case READINPUT :
            switch (ProcessInput( command.GetStringNext() )) {
              case C_ERR  : if ( exitOnError_ ) return C_ERR; break;
              case C_QUIT : return C_QUIT; break;
              default     : break;
            } 
            break;
          case CREATE      : err = Create_DataFile( command ); break;
          case PRECISION   : err = Precision( command ); break;
          case DATAFILE    : err = DFL.ProcessDataFileArgs( command ); break;
          case SYSTEM      : system( command.ArgLine() ); break;
          case RUN         : Run(); break;
          case RUN_ANALYSIS:
            // If only 1 arg (the command) run all analyses in list
            if (command.Nargs() == 1)  
              analysisList.DoAnalyses(&DFL);
            else
              err = CrdAnalyze(command);
            mprintf("Analysis complete. Use 'writedata' to write datafiles to disk.\n");
            break;
          case WRITEDATA   : if (worldrank == 0) DFL.Write(); break;
          case QUIT        : return C_QUIT; break;
        }
        break;
      case DispatchObject::DEPRECATED: 
        mprintf("Warning: %s is deprecated.\n", command.Command()); 
        break;
      default: mprintf("Dispatch type is currently not handled.\n");
    }
    if (err != 0 && exitOnError_) return C_ERR;
  } else { // Command not recognized
    if (exitOnError_) return C_ERR;
  }
  return C_OK;
}

// -----------------------------------------------------------------------------
// Cpptraj::Run()
int Cpptraj::Run() {
  int err = 0;
  ++nrun_;
  switch ( trajinList.Mode() ) {
    case TrajinList::NORMAL   : err = RunNormal(); break;
    case TrajinList::ENSEMBLE : err = RunEnsemble(); break;
    default: 
      mprinterr("No trajectories loaded. Exiting.\n");
      err = 1; 
  }
  return err;
}

// Cpptraj::RunEnsemble()
int Cpptraj::RunEnsemble() {
  FrameArray FrameEnsemble;

  mprintf("\nINPUT ENSEMBLE:\n");
  // Ensure all ensembles are of the same size
  int ensembleSize = -1;
  for (TrajinList::const_iterator traj = trajinList.begin(); traj != trajinList.end(); ++traj) 
  {
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    if (ensembleSize == -1)
      ensembleSize = mtraj->EnsembleSize();
    else if (ensembleSize != mtraj->EnsembleSize()) {
      mprinterr("Error: Ensemble size (%i) does not match first ensemble size (%i).\n",
                mtraj->EnsembleSize(), ensembleSize);
      return 1;
    }
    // Perform ensemble setup - this also resizes FrameEnsemble
    if ( mtraj->EnsembleSetup( FrameEnsemble ) ) return 1;
  }
  mprintf("  Ensemble size is %i\n", ensembleSize);

  // Calculate frame division among trajectories
  trajinList.List();
  int maxFrames = trajinList.MaxFrames();
  // Parameter file information
  parmFileList.List();
  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames.List();

  // Allocate an ActionList, TrajoutList, and DataSetList for each
  // member of the ensemble.
  std::vector<ActionList> ActionEnsemble( ensembleSize );
  std::vector<TrajoutList> TrajoutEnsemble( ensembleSize );
  std::vector<DataSetList> DataSetEnsemble( ensembleSize );

  // Set up output trajectories for each member of the ensemble
  for (ArgsArray::iterator targ = trajoutArgs_.begin(); targ != trajoutArgs_.end(); ++targ)
  {
    for (int member = 0; member < ensembleSize; ++member) 
      TrajoutEnsemble[member].AddEnsembleTrajout( *targ, parmFileList, member );
  }
  mprintf("\n");
  for (int member = 0; member < ensembleSize; ++member) {
    mprintf("OUTPUT TRAJECTORIES Member %i:\n", member);
    TrajoutEnsemble[member].List();
  }

  // TODO: One loop over member?
  for (int member = 0; member < ensembleSize; ++member) {
    mprintf("***** ENSEMBLE MEMBER %i: ", member);
    // Set max frames in the data set list and allocate
    DataSetEnsemble[member].SetMax( maxFrames );
    DataSetEnsemble[member].AllocateSets();
    // Initialize actions 
    for (ArgsArray::iterator aarg = actionArgs_.begin(); aarg != actionArgs_.end(); ++aarg)
    {
      DispatchObject::TokenPtr dispatchToken = SearchToken( *aarg );
      if ( dispatchToken != 0 ) {
        // Create copy of arg list so that args remain unmarked for next member
        ArgList command = *aarg;
        if (ActionEnsemble[member].AddAction( dispatchToken->Alloc, command, &parmFileList,
                                          &refFrames, &(DataSetEnsemble[member]), &DFL ))
          return 1;
      }
    }
  }
      
  // ========== A C T I O N  P H A S E ==========
  int lastPindex=-1;          // Index of the last loaded parm file
  int readSets = 0;
  int actionSet = 0;
  bool hasVelocity = false;
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN ENSEMBLE PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList.begin();
                                   traj != trajinList.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->FullTrajStr());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());

    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (hasVelocity != (*traj)->HasVelocity()))
      FrameEnsemble.SetupFrames(CurrentParm->Atoms(), (*traj)->HasVelocity());
    hasVelocity = (*traj)->HasVelocity();

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( refFrames.ActiveReference() );
      // Set up actions for this parm
      bool setupOK = true;
      for (int member = 0; member < ensembleSize; ++member) {
        if (ActionEnsemble[member].SetupActions( &CurrentParm )) {
          mprintf("Warning: Ensemble member %i: Could not set up actions for %s: skipping.\n",
                  member, CurrentParm->c_str());
          setupOK = false;
        }
      }
      if (!setupOK) continue;
      lastPindex = CurrentParm->Pindex();
    }

    // Loop over every collection of frames in the ensemble
    (*traj)->PrintInfoLine();
    Trajin_Multi* mtraj = (Trajin_Multi*)*traj;
    while ( mtraj->GetNextEnsemble(FrameEnsemble) ) {
      if (!mtraj->BadEnsemble()) {
        // Loop over all members of the ensemble
        for (int member = 0; member < ensembleSize; ++member) {
          // Get this members current position
          int pos = mtraj->EnsemblePosition( member );
          // Since Frame can be modified by actions, save original and use CurrentFrame
          Frame* CurrentFrame = &(FrameEnsemble[member]);
          // Perform Actions on Frame
          bool suppress_output = ActionEnsemble[pos].DoActions(&CurrentFrame, actionSet);
          // Do Output
          if (!suppress_output) 
            TrajoutEnsemble[pos].Write(actionSet, CurrentParm, CurrentFrame);
        } // END loop over ensemble
      } else {
        mprinterr("Error: Could not read frame %i for ensemble.\n", actionSet + 1);
      }
      // Increment frame counter
      ++actionSet;
    }

    // Close the trajectory file
    (*traj)->EndTraj();
    // Update how many frames have been processed.
    readSets += (*traj)->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output trajectories
  for (int member = 0; member < ensembleSize; ++member)
    TrajoutEnsemble[member].Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  for (int member = 0; member < ensembleSize; ++member)
    actionList.Print( );

  // Sync DataSets and print DataSet information
  // TODO - Also have datafilelist call a sync??
  for (int member = 0; member < ensembleSize; ++member) {
    DataSetEnsemble[member].Sync();
    DataSetEnsemble[member].sort();
    mprintf("\nENSEMBLE MEMBER %i DATASETS:\n",member);
    DataSetEnsemble[member].List();
  }

  // Print Datafile information
  DFL.List();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();

  return 0;
}

// Cpptraj::RunNormal()
/** Process trajectories in trajinList. Each frame in trajinList is sent
 *  to the actions in actionList for processing.
 */
int Cpptraj::RunNormal() {
  int actionSet=0;            // Internal data frame
  int readSets=0;             // Number of frames actually read
  int lastPindex=-1;          // Index of the last loaded parm file
  Frame TrajFrame;            // Original Frame read in from traj

  // ========== S E T U P   P H A S E ========== 
  // Parameter file information
  parmFileList.List();
  // Input coordinate file information
  trajinList.List();
  // Print reference information 
  mprintf("\nREFERENCE COORDS:\n");
  refFrames.List();
  // Output traj
  mprintf("\nOUTPUT TRAJECTORIES:\n");
  trajoutList.List();
  // Allocate DataSets in the master DataSetList based on # frames to be read
  DSL.AllocateSets(); 
  
  // ========== A C T I O N  P H A S E ==========
  // Loop over every trajectory in trajFileList
  rprintf("BEGIN TRAJECTORY PROCESSING:\n");
  for ( TrajinList::const_iterator traj = trajinList.begin();
                                   traj != trajinList.end(); ++traj)
  {
    // Open up the trajectory file. If an error occurs, bail 
    if ( (*traj)->BeginTraj(showProgress_) ) {
      mprinterr("Error: Could not open trajectory %s.\n",(*traj)->FullTrajStr());
      break;
    }
    // Set current parm from current traj.
    Topology* CurrentParm = (*traj)->TrajParm();
    // Check if parm has changed
    bool parmHasChanged = (lastPindex != CurrentParm->Pindex());

    // If Parm has changed or trajectory velocity status has changed,
    // reset the frame.
    if (parmHasChanged || (TrajFrame.HasVelocity() != (*traj)->HasVelocity()))
      TrajFrame.SetupFrameV(CurrentParm->Atoms(), (*traj)->HasVelocity());

    // If Parm has changed, reset actions for new topology.
    if (parmHasChanged) {
      // Set active reference for this parm
      CurrentParm->SetReferenceCoords( refFrames.ActiveReference() );
      // Set up actions for this parm
      if (actionList.SetupActions( &CurrentParm )) {
        mprintf("WARNING: Could not set up actions for %s: skipping.\n",
                CurrentParm->c_str());
        continue;
      }
      lastPindex = CurrentParm->Pindex();
    }

    // Loop over every Frame in trajectory
    (*traj)->PrintInfoLine();
    while ( (*traj)->GetNextFrame(TrajFrame) ) {
      // Since Frame can be modified by actions, save original and use CurrentFrame
      Frame* CurrentFrame = &TrajFrame;
      // Perform Actions on Frame
      bool suppress_output = actionList.DoActions(&CurrentFrame, actionSet);
      // Do Output
      if (!suppress_output)
        trajoutList.Write(actionSet, CurrentParm, CurrentFrame);
      // Increment frame counter
      ++actionSet; 
    }

    // Close the trajectory file
    (*traj)->EndTraj();
    // Update how many frames have been processed.
    readSets += (*traj)->NumFramesProcessed();
    mprintf("\n");
  } // End loop over trajin
  rprintf("Read %i frames and processed %i frames.\n",readSets,actionSet);

  // Close output traj
  trajoutList.Close();

  // ========== A C T I O N  O U T P U T  P H A S E ==========
  actionList.Print( );

  // Sync DataSets and print DataSet information
  DSL.Sync();

  // ========== A N A L Y S I S  P H A S E ==========
  mprintf("\nDATASETS:\n");
  if (!analysisList.Empty()) {
    DSL.List();
    analysisList.DoAnalyses(&DFL);
    // DEBUG: DataSets, post-Analysis
    mprintf("\nDATASETS AFTER ANALYSIS:\n");
  }
  DSL.List();

  // ========== D A T A  W R I T E  P H A S E ==========
  // Print Datafile information
  DFL.List();
  // Only Master does DataFile output
  if (worldrank==0)
    DFL.Write();
 
  return 0;
}
