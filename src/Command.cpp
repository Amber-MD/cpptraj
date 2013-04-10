#include "Command.h"
#include "CpptrajStdio.h"
// INC_ACTION==================== ALL ACTION CLASSES GO HERE ===================
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_Angle.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_DSSP.h"
#include "Action_Center.h"
#include "Action_Hbond.h"
#include "Action_Image.h"
#include "Action_Surf.h"
#include "Action_Radgyr.h"
#include "Action_Mask.h"
#include "Action_Closest.h"
#include "Action_NAstruct.h"
#include "Action_Pucker.h"
#include "Action_Outtraj.h"
#include "Action_Average.h"
#include "Action_Radial.h"
#include "Action_DistRmsd.h"
#include "Action_Jcoupling.h"
#include "Action_Pairwise.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
#include "Action_DihedralScan.h"
#include "Action_Rotdif.h"
#include "Action_RunningAvg.h"
#include "Action_AtomicFluct.h"
#include "Action_Watershell.h"
#include "Action_Contacts.h"
#include "Action_Vector.h"
#include "Action_Principal.h"
#include "Action_Matrix.h"
#include "Action_LIE.h"
#include "Action_Grid.h"
#include "Action_GridFreeEnergy.h"
#include "Action_Dipole.h"
#include "Action_Projection.h"
#include "Action_ClusterDihedral.h"
#include "Action_Unwrap.h"
#include "Action_Diffusion.h"
#include "Action_DNAionTracker.h"
#include "Action_Scale.h"
#include "Action_RandomizeIons.h"
#include "Action_AutoImage.h"
#include "Action_STFC_Diffusion.h"
#include "Action_AtomicCorr.h"
#include "Action_Bounds.h"
#include "Action_Rotate.h"
#include "Action_Translate.h"
#include "Action_Box.h"
#include "Action_CreateCrd.h"
#include "Action_MultiDihedral.h"
#include "Action_MakeStructure.h"
#include "Action_SymmetricRmsd.h"
#include "Action_Volmap.h"
#include "Action_Spam.h"
#include "Action_Temperature.h"

// INC_ANALYSIS================= ALL ANALYSIS CLASSES GO HERE ==================
#include "Analysis_Hist.h"
#include "Analysis_Corr.h"
#include "Analysis_Matrix.h"
#include "Analysis_Timecorr.h"
#include "Analysis_IRED.h"
#include "Analysis_Modes.h"
#include "Analysis_CrankShaft.h"
#include "Analysis_Statistics.h"
#include "Analysis_CrossCorr.h"
#include "Analysis_AutoCorr.h"
#include "Analysis_Lifetime.h"
#include "Analysis_FFT.h"
#include "Analysis_CrdFluct.h"
#include "Analysis_RmsAvgCorr.h"
#include "Analysis_Rms2d.h"
#include "Analysis_Clustering.h"
#include "Analysis_RunningAvg.h"

// ====================== CPPTRAJ COMMANDS HELP ================================
static void Help_Help() {
  mprintf("\t{[<cmd>] | General | Action | Analysis | Topology | Trajectory}\n");
  mprintf("\tWith no arguments list all known commands, otherwise display help for\n");
  mprintf("\tcommand <cmd>. If General/Action/Analysis/Topology/Trajectory specified\n");
  mprintf("\tlist commands only in that category.\n");
}

static void Help_System() {
  mprintf("\tCall command from system.\n");
}

static void Help_NoProgress() {
  mprintf("\tDo not print progress while reading in trajectories.\n");
}

static void Help_NoExitOnError() {
  mprintf("\tDo not exit when errors are encountered. This is the default\n");
  mprintf("\tin interactive mode.\n");
}

static void Help_Run() {
  mprintf("\tProcess all trajectories currently in input trajectory list.\n");
  mprintf("\tAll actions in action list will be run on each frame.\n");
  mprintf("\tIf not processing ensemble input, all analyses in analysis\n");
  mprintf("\tlist will be run after trajectory processing.\n");
}

static void Help_Quit() {
  mprintf("\tExit CPPTRAJ\n");
}

static const char TypeList[] =
  "(<type> = actions,trajin,trajout,ref,parm,analysis,datafile,dataset)";

static void Help_List() {
  mprintf("\t[<type>] %s\n", TypeList);
  mprintf("\tList currently loaded objects of the specified type. If no type is given\n");
  mprintf("\tlist all loaded objects.\n");
}

static void Help_Debug() {
  mprintf("\t[<type>] <#> %s\n", TypeList);
  mprintf("\tSet debug level for new objects of the specified type. If no type is given\n");
  mprintf("\tset debug level for all new objects. Does not affect current objects.\n");
}

static void Help_Clear() {
  mprintf("\t[ {all | <type>} ] %s\n", TypeList);
  mprintf("\tClear currently loaded objects of the specified type. If 'all' is specified\n");
  mprintf("\tclear all loaded objects.\n");
}

static void Help_ActiveRef() {
  mprintf("\t<#>\n");
  mprintf("\tSet the reference structure to be used for coordinate-based mask parsing.\n");
  mprintf("\t<#> starts from 0 (first loaded reference).\n");
}

static void Help_Create_DataFile() {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n");
  mprintf("\tAdd a file with specified data sets to the data file list. Does not\n");
  mprintf("\timmediately write the data.\n");
}

static void Help_DataFile() {
  mprintf("\t<data filename> <datafile cmd>\n");
  mprintf("\tPass <datafile cmd> to specified data file currently in data file list.\n");
}

static void Help_ReadData() {
  mprintf("\t<filename>\n");
  mprintf("\tRead data from <filename> into data sets.\n");
}

static void Help_ReadInput() {
  mprintf("\t<filename>\n");
  mprintf("\tRead commands from <filename>\n");
}

static void Help_Write_DataFile() {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n");
  mprintf("\tWrite specified data sets to <filename> immediately.\n");
}

static void Help_WriteData() {
  mprintf("\tWrite all files currently in the data file list.\n");
}

static void Help_Precision() {
  mprintf("\t{<filename> | <dataset arg>} [<width>] [<precision>]\n");
  mprintf("\tSet precision for all datasets in datafile <filename> or dataset(s)\n");
  mprintf("\tspecified by <dataset arg> to <width>.<precision>. If width/precision\n");
  mprintf("\tnot specified default to 12.4\n");
}

static void Help_Select() {
  mprintf("\t[<parmindex>] <mask>\n");
  mprintf("\tShow atom numbers selected by <mask> for parm <parmindex>\n");
  mprintf("\t(default first parm)\n");
}

static void Help_SelectDS() {
  mprintf("\t<dataset selection>\n");
  mprintf("\tShow results of data set selection. Data set selection format is:\n");
  mprintf("\t\t<name>[<aspect]:<idx range>\n");
  mprintf("\tWhere '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n");
  mprintf("\tand <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n");
  mprintf("\tThe aspect and index portions may be optional. An asterisk '*' may be used as\n");
  mprintf("\ta wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

static void Help_Trajin() {
  mprintf("\t<filename> {[<start>] [<stop> | last] [offset]} | lastframe\n");
  mprintf("\t           [parm <parmfile> | parmindex <#>]\n");
  mprintf("\t           [ remdtraj [remdtrajtemp <T> | remdtrajidx <#>]\n");
  mprintf("\t           [trajnames <rep1>,<rep2>,...,<repN> ] ]\n");
  mprintf("\tLoad trajectory specified by <filename> to the input trajectory list.\n");
}

static void Help_Ensemble() {
  mprintf("\t<file0> {[<start>] [<stop> | last] [offset]} | lastframe\n");
  mprintf("\t        [parm <parmfile> | parmindex <#>]\n");
  mprintf("\t        [trajnames <file1>,<file2>,...,<fileN>\n");
  mprintf("\tLoad an ensemble of trajectories starting with <file0> that will be processed together.\n");
}

static void Help_Trajout() {
  mprintf("\t<filename> [<fileformat>] [append] [nobox]\n");
  mprintf("\t           [parm <parmfile> | parmindex <#>] [onlyframes <range>] [title <title>]\n");
  mprintf("\t           %s\n", ActionFrameCounter::HelpText);
  mprintf("\t           [ <Format Options> ]\n");
  mprintf("\tSpecify output trajectory.\n");
}

static void Help_Reference() {
  mprintf("\t<filename> [<frame#>] [<mask>] [TAG] [lastframe]\n");
  mprintf("\t           [average [<stop>] [<offset>]]\n");
  mprintf("\tLoad trajectory <filename> as a reference frame.\n");
}

static void Help_Parm() {
  mprintf("\t<filename> [<tag>] [nobondsearch | bondsearch [<offset>]]\n");
  mprintf("\tAdd <filename> to the topology list.\n");
}
static void Help_ParmInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint information on topology <parmindex> (0 by default). If <mask> is given\n");
  mprintf("\tprint info on atoms in mask. If no mask given print overall information.\n");
}

static void Help_ParmWrite() {
  mprintf("\tout <filename> [<parmindex>]\n");
  mprintf("\tWrite topology <parmindex> to <filename> as an Amber Topology file.\n");
}

static void Help_ParmStrip() {
  mprintf("\t<mask> [<parmindex>]\n");
  mprintf("\tStrip atoms in mask from topology <parmindex>.\n");
}

static void Help_ParmBox() {
  mprintf("\t[<parmindex>] [x <xval>] [y <yval>] [z <zval>]");
  mprintf("\t              [alpha <a>] [beta <b>] [gamma <g>] [nobox]\n");
  mprintf("\tSet the specified topology box info to what is specified. If nobox, remove box info.\n");
}

static void Help_Solvent() {
  mprintf("\t[<parmindex>] <mask>\n");
  mprintf("\tSet solvent for the specified topology (default 0) based on <mask>\n");
}

static void Help_BondInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint bond information of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_ChargeInfo() {
  mprintf("\t[<parmindex>] <mask>\n");
  mprintf("\tPrint the total charge of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_ResInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint information for residues in <mask> for topology <parmindex> (0 by default).\n");
}
static void Help_MolInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n");
  mprintf("\tPrint information for molecules in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_LoadCrd() {
  mprintf("\t<filename> [parm <parm> | parmindex<#>] [<trajin args>] [<name>]\n");
  mprintf("\tLoad trajectory <filename> as a COORDS data set named <name> (default <filename>).\n");
}

static void Help_CrdAction() {
  mprintf("\t<crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n");
  mprintf("\tPerform action <actioncmd> on COORDS data set <crd set>.\n");
}

static void Help_CrdOut() {
  mprintf("\t<crd set> <filename> [<trajout args>] [crdframes <start>,<stop>,<offset>]\n");
  mprintf("\tWrite COORDS data set <crd set> to trajectory file <filename>\n");
}

static void Help_RunAnalysis() {
  mprintf("\t[<analysis> [<analysis args>]]\n");
  mprintf("\tIf specified alone, run all analyses in the analysis list.\n");
  mprintf("\tOtherwise run the specified analysis immediately.\n");
}

// ================ LIST OF ALL COMMANDS =======================================
/** Ideally keep this array first sorted by type (1st field), then 
  * alphabetically by command string (2nd field).
  */
const DispatchObject::Token Command::Commands[] = {
  // GENERAL COMMANDS
  { DispatchObject::GENERAL, "activeref",     0, Help_ActiveRef,       ACTIVEREF    },
  { DispatchObject::GENERAL, "clear",         0, Help_Clear,           CLEAR        },
  { DispatchObject::GENERAL, "crdaction",     0, Help_CrdAction,       CRDACTION    },
  { DispatchObject::GENERAL, "crdout",        0, Help_CrdOut,          CRDOUT       },
  { DispatchObject::GENERAL, "create",        0, Help_Create_DataFile, CREATE       },
  { DispatchObject::GENERAL, "datafile",      0, Help_DataFile,        DATAFILE     },
  { DispatchObject::GENERAL, "debug",         0, Help_Debug,           DEBUG        },
  { DispatchObject::GENERAL, "exit" ,         0, Help_Quit,            QUIT         },
  { DispatchObject::GENERAL, "gnuplot",       0, Help_System,          SYSTEM       },
  { DispatchObject::GENERAL, "go",            0, Help_Run,             RUN          },
  { DispatchObject::GENERAL, "head",          0, Help_System,          SYSTEM       },
  { DispatchObject::GENERAL, "help",          0, Help_Help,            HELP         },
  { DispatchObject::GENERAL, "list",          0, Help_List,            LIST         },
  { DispatchObject::GENERAL, "loadcrd",       0, Help_LoadCrd,         LOADCRD      },
  { DispatchObject::GENERAL, "ls",            0, Help_System,          SYSTEM       },
  { DispatchObject::GENERAL, "noexitonerror", 0, Help_NoExitOnError,   NOEXITERR    },
  { DispatchObject::GENERAL, "noprogress",    0, Help_NoProgress,      NOPROG       },
  { DispatchObject::GENERAL, "precision",     0, Help_Precision,       PRECISION    },
  { DispatchObject::GENERAL, "prnlev",        0, Help_Debug,           DEBUG        },
  { DispatchObject::GENERAL, "pwd",           0, Help_System,          SYSTEM       },
  { DispatchObject::GENERAL, "quit" ,         0, Help_Quit,            QUIT         },
  { DispatchObject::GENERAL, "readdata",      0, Help_ReadData,        READDATA     },
  { DispatchObject::GENERAL, "readinput",     0, Help_ReadInput,       READINPUT    },
  { DispatchObject::GENERAL, "run"   ,        0, Help_Run,             RUN          },
  { DispatchObject::GENERAL, "runanalysis",   0, Help_RunAnalysis,     RUN_ANALYSIS },
  { DispatchObject::GENERAL, "select",        0, Help_Select,          SELECT       },
  { DispatchObject::GENERAL, "selectds",      0, Help_SelectDS,        SELECTDS     },
  { DispatchObject::GENERAL, "write",         0, Help_Write_DataFile,  WRITE        },
  { DispatchObject::GENERAL, "writedata",     0, Help_WriteData,       WRITEDATA    },
  { DispatchObject::GENERAL, "xmgrace",       0, Help_System,          SYSTEM       },
  // TRAJECTORY COMMANDS
  { DispatchObject::TRAJ,   "ensemble",      0, Help_Ensemble,        TRAJIN     },
  { DispatchObject::TRAJ,   "reference",     0, Help_Reference,       REFERENCE  },
  { DispatchObject::TRAJ,   "trajin",        0, Help_Trajin,          TRAJIN     },
  { DispatchObject::TRAJ,   "trajout",       0, Help_Trajout,         TRAJOUT    },
  // TOPOLOGY COMMANDS
  { DispatchObject::PARM,    "bondinfo",      0, Help_BondInfo,        BONDINFO   },
  { DispatchObject::PARM,    "charge",        0, Help_ChargeInfo,      CHARGEINFO },
  { DispatchObject::PARM,    "molinfo",       0, Help_MolInfo,         MOLINFO    },
  { DispatchObject::PARM,    "parm",          0, Help_Parm,            LOADPARM   },
  { DispatchObject::PARM,    "parmbondinfo",  0, Help_BondInfo,        BONDINFO   },
  { DispatchObject::PARM,    "parmbox",       0, Help_ParmBox,         PARMBOX    },
  { DispatchObject::PARM,    "parminfo",      0, Help_ParmInfo,        PARMINFO   },
  { DispatchObject::PARM,    "parmmolinfo",   0, Help_MolInfo,         MOLINFO    },
  { DispatchObject::PARM,    "parmresinfo",   0, Help_ResInfo,         RESINFO    },
  { DispatchObject::PARM,    "parmstrip",     0, Help_ParmStrip,       PARMSTRIP  },
  { DispatchObject::PARM,    "parmwrite",     0, Help_ParmWrite,       PARMWRITE  },
  { DispatchObject::PARM,    "resinfo",       0, Help_ResInfo,         RESINFO    },
  { DispatchObject::PARM,    "solvent",       0, Help_Solvent,         SOLVENT    },
  // INC_ACTION: ACTION COMMANDS
  { DispatchObject::ACTION, "angle", Action_Angle::Alloc, Action_Angle::Help, 0 },
  { DispatchObject::ACTION, "atomiccorr", Action_AtomicCorr::Alloc, Action_AtomicCorr::Help, 0 },
  { DispatchObject::ACTION, "atomicfluct", Action_AtomicFluct::Alloc, Action_AtomicFluct::Help, 0 },
  { DispatchObject::ACTION, "atommap", Action_AtomMap::Alloc, Action_AtomMap::Help, 0 },
  { DispatchObject::ACTION, "autoimage", Action_AutoImage::Alloc, Action_AutoImage::Help, 0 },
  { DispatchObject::ACTION, "average", Action_Average::Alloc, Action_Average::Help, 0 },
  { DispatchObject::ACTION, "bounds", Action_Bounds::Alloc, Action_Bounds::Help, 0 },
  { DispatchObject::ACTION, "box", Action_Box::Alloc, Action_Box::Help, 0 },
  { DispatchObject::ACTION, "center", Action_Center::Alloc, Action_Center::Help, 0 },
  { DispatchObject::ACTION, "check", Action_CheckStructure::Alloc, Action_CheckStructure::Help, 0 },
  { DispatchObject::ACTION, "checkstructure", Action_CheckStructure::Alloc, Action_CheckStructure::Help, 0 },
  { DispatchObject::ACTION, "closest", Action_Closest::Alloc, Action_Closest::Help, 0 },
  { DispatchObject::ACTION, "clusterdihedral", Action_ClusterDihedral::Alloc, Action_ClusterDihedral::Help, 0 },
  { DispatchObject::ACTION, "contacts", Action_Contacts::Alloc, Action_Contacts::Help, 0 },
  { DispatchObject::ACTION, "createcrd", Action_CreateCrd::Alloc, Action_CreateCrd::Help, 0 },
  { DispatchObject::ACTION, "diffusion", Action_Diffusion::Alloc, Action_Diffusion::Help, 0 },
  { DispatchObject::ACTION, "dihedral", Action_Dihedral::Alloc, Action_Dihedral::Help, 0 },
  { DispatchObject::ACTION, "dihedralscan", Action_DihedralScan::Alloc, Action_DihedralScan::Help, 0 },
  { DispatchObject::ACTION, "dipole", Action_Dipole::Alloc, Action_Dipole::Help, 0 },
  { DispatchObject::ACTION, "distance", Action_Distance::Alloc, Action_Distance::Help, 0 },
//  { DispatchObject::ACTION, "dnaiontracker", Action_DNAionTracker::Alloc, Action_DNAionTracker::Help, 0 },
  { DispatchObject::ACTION, "drms", Action_DistRmsd::Alloc, Action_DistRmsd::Help, 0 },
  { DispatchObject::ACTION, "drmsd", Action_DistRmsd::Alloc, Action_DistRmsd::Help, 0 },
//  { DispatchObject::ACTION, "gfe", Action_GridFreeEnergy::Alloc, Action_GridFreeEnergy::Help, 0 },
  { DispatchObject::ACTION, "grid", Action_Grid::Alloc, Action_Grid::Help, 0 },
  { DispatchObject::ACTION, "hbond", Action_Hbond::Alloc, Action_Hbond::Help, 0 },
  { DispatchObject::ACTION, "image", Action_Image::Alloc, Action_Image::Help, 0 },
  { DispatchObject::ACTION, "jcoupling", Action_Jcoupling::Alloc, Action_Jcoupling::Help, 0 },
  { DispatchObject::ACTION, "lie", Action_LIE::Alloc, Action_LIE::Help, 0 },
  { DispatchObject::ACTION, "makestructure", Action_MakeStructure::Alloc, Action_MakeStructure::Help, 0 },
  { DispatchObject::ACTION, "mask", Action_Mask::Alloc, Action_Mask::Help, 0 },
  { DispatchObject::ACTION, "matrix", Action_Matrix::Alloc, Action_Matrix::Help, 0 },
  { DispatchObject::ACTION, "molsurf", Action_Molsurf::Alloc, Action_Molsurf::Help, 0 },
  { DispatchObject::ACTION, "multidihedral", Action_MultiDihedral::Alloc, Action_MultiDihedral::Help, 0 },
  { DispatchObject::ACTION, "nastruct", Action_NAstruct::Alloc, Action_NAstruct::Help, 0 },
  { DispatchObject::ACTION, "outtraj", Action_Outtraj::Alloc, Action_Outtraj::Help, 0 },
  { DispatchObject::ACTION, "pairwise", Action_Pairwise::Alloc, Action_Pairwise::Help, 0 },
  { DispatchObject::ACTION, "principal", Action_Principal::Alloc, Action_Principal::Help, 0 },
  { DispatchObject::ACTION, "projection", Action_Projection::Alloc, Action_Projection::Help, 0 },
  { DispatchObject::ACTION, "pucker", Action_Pucker::Alloc, Action_Pucker::Help, 0 },
  { DispatchObject::ACTION, "radgyr", Action_Radgyr::Alloc, Action_Radgyr::Help, 0 },
  { DispatchObject::ACTION, "radial", Action_Radial::Alloc, Action_Radial::Help, 0 },
  { DispatchObject::ACTION, "randomizeions", Action_RandomizeIons::Alloc, Action_RandomizeIons::Help, 0 },
  { DispatchObject::ACTION, "rms", Action_Rmsd::Alloc, Action_Rmsd::Help, 0 },
  { DispatchObject::ACTION, "rmsd", Action_Rmsd::Alloc, Action_Rmsd::Help, 0 },
  { DispatchObject::ACTION, "rog", Action_Radgyr::Alloc, Action_Radgyr::Help, 0 },
  { DispatchObject::ACTION, "rotate", Action_Rotate::Alloc, Action_Rotate::Help, 0 },
  { DispatchObject::ACTION, "rotdif", Action_Rotdif::Alloc, Action_Rotdif::Help, 0 },
  { DispatchObject::ACTION, "runavg", Action_RunningAvg::Alloc, Action_RunningAvg::Help, 0 },
  { DispatchObject::ACTION, "runningaverage", Action_RunningAvg::Alloc, Action_RunningAvg::Help, 0 },
  { DispatchObject::ACTION, "scale", Action_Scale::Alloc, Action_Scale::Help, 0 },
  { DispatchObject::ACTION, "secstruct", Action_DSSP::Alloc, Action_DSSP::Help, 0 },
  { DispatchObject::ACTION, "spam", Action_Spam::Alloc, Action_Spam::Help, 0 },
  { DispatchObject::ACTION, "stfcdiffusion", Action_STFC_Diffusion::Alloc, Action_STFC_Diffusion::Help, 0 },
  { DispatchObject::ACTION, "strip", Action_Strip::Alloc, Action_Strip::Help, 0 },
  { DispatchObject::ACTION, "surf", Action_Surf::Alloc, Action_Surf::Help, 0 },
  { DispatchObject::ACTION, "symmrmsd", Action_SymmetricRmsd::Alloc, Action_SymmetricRmsd::Help, 0 },
  { DispatchObject::ACTION, "temperature", Action_Temperature::Alloc, Action_Temperature::Help, 0 },
  { DispatchObject::ACTION, "trans", Action_Translate::Alloc, Action_Translate::Help, 0 },
  { DispatchObject::ACTION, "translate", Action_Translate::Alloc, Action_Translate::Help, 0 },
  { DispatchObject::ACTION, "unstrip", Action_Unstrip::Alloc, Action_Unstrip::Help, 0 },
  { DispatchObject::ACTION, "unwrap", Action_Unwrap::Alloc, Action_Unwrap::Help, 0 },
  { DispatchObject::ACTION, "vector", Action_Vector::Alloc, Action_Vector::Help, 0 },
  { DispatchObject::ACTION, "watershell", Action_Watershell::Alloc, Action_Watershell::Help, 0 },
  { DispatchObject::ACTION, "volmap", Action_Volmap::Alloc, Action_Volmap::Help, 0},
  // INC_ANALYSIS: ANALYSIS COMMANDS
  { DispatchObject::ANALYSIS, "2drms", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, 0 },
  { DispatchObject::ANALYSIS, "autocorr", Analysis_AutoCorr::Alloc, Analysis_AutoCorr::Help, 0 },
  { DispatchObject::ANALYSIS, "cluster", Analysis_Clustering::Alloc, Analysis_Clustering::Help, 0 },
  { DispatchObject::ANALYSIS, "corr", Analysis_Corr::Alloc, Analysis_Corr::Help, 0 },
  { DispatchObject::ANALYSIS, "correlationcoe", Analysis_Corr::Alloc, Analysis_Corr::Help, 0 },
  { DispatchObject::ANALYSIS, "crank", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, 0 },
  { DispatchObject::ANALYSIS, "crankshaft", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, 0 },
  { DispatchObject::ANALYSIS, "crdfluct", Analysis_CrdFluct::Alloc, Analysis_CrdFluct::Help, 0 },
  { DispatchObject::ANALYSIS, "crosscorr", Analysis_CrossCorr::Alloc, Analysis_CrossCorr::Help, 0 },
  { DispatchObject::ANALYSIS, "diagmatrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, 0 },
  { DispatchObject::ANALYSIS, "fft", Analysis_FFT::Alloc, Analysis_FFT::Help, 0 },
  { DispatchObject::ANALYSIS, "hist", Analysis_Hist::Alloc, Analysis_Hist::Help, 0 },
  { DispatchObject::ANALYSIS, "histogram", Analysis_Hist::Alloc, Analysis_Hist::Help, 0 },
  { DispatchObject::ANALYSIS, "ired", Analysis_IRED::Alloc, Analysis_IRED::Help, 0 },
  { DispatchObject::ANALYSIS, "lifetime", Analysis_Lifetime::Alloc, Analysis_Lifetime::Help, 0 },
  { DispatchObject::ANALYSIS, "matrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, 0 },
  { DispatchObject::ANALYSIS, "modes", Analysis_Modes::Alloc, Analysis_Modes::Help, 0 },
  { DispatchObject::ANALYSIS, "rms2d", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, 0 },
  { DispatchObject::ANALYSIS, "rmsavgcorr", Analysis_RmsAvgCorr::Alloc, Analysis_RmsAvgCorr::Help, 0 },
  { DispatchObject::ANALYSIS, "stat", Analysis_Statistics::Alloc, Analysis_Statistics::Help, 0 },
  { DispatchObject::ANALYSIS, "statistics", Analysis_Statistics::Alloc, Analysis_Statistics::Help, 0 },
  { DispatchObject::ANALYSIS, "timecorr", Analysis_Timecorr::Alloc, Analysis_Timecorr::Help, 0 },
  { DispatchObject::ANALYSIS, "runningavg", Analysis_RunningAvg::Alloc, Analysis_RunningAvg::Help, 0 },
  // DEPRECATED COMMANDS
  { DispatchObject::DEPRECATED, "molsearch",    0, 0, 0 },
  { DispatchObject::DEPRECATED, "nomolsearch",  0, 0, 0 },
  { DispatchObject::DEPRECATED, "bondsearch",   0, 0, 0 },
  { DispatchObject::DEPRECATED, "nobondsearch", 0, 0, 0 },
  { DispatchObject::NONE      , 0,              0, 0, 0 }
};

/// Strings that correspond to enumerated type DispatchObject::DispatchType
static const char* CommandTitle[] = { 0, "Topology", "Trajectory", "Action",
  "Analysis", "General", "Deprecated" };

/** List all commands of the given type, or all commands if type
  * is NONE.
  */
void Command::List(DispatchObject::DispatchType dtype) {
  DispatchObject::DispatchType lastType = DispatchObject::NONE;
  int col = 0;
  for (DispatchObject::TokenPtr token = Commands; 
                                token->Type != DispatchObject::DEPRECATED; ++token)
  {
    DispatchObject::DispatchType currentType = token->Type;
    if (dtype != DispatchObject::NONE && dtype != currentType) continue;
    // Command type title
    if (currentType != lastType) {
      if (col != 0) mprintf("\n");
      mprintf("%s Commands:\n", CommandTitle[currentType]);
      lastType = currentType;
      col = 0;
    }
    if (col == 0) mprintf("\t");
    mprintf("%s  ", token->Cmd);
    ++col;
    if (col == 8) {
      mprintf("\n");
      col = 0;
    }
  }
  mprintf("\n");
}

/** Search Commands list for a specific type of command. */
DispatchObject::TokenPtr Command::SearchTokenType(DispatchObject::DispatchType dtype,
                                                  ArgList& argIn)
{
  for (DispatchObject::TokenPtr token = Commands;
                                token->Type != DispatchObject::NONE; ++token)
  {
    if (dtype != token->Type) continue;
    if (argIn.CommandIs( token->Cmd )) return token;
  }
  mprintf("[%s]: Command not found.\n", argIn.Command());
  return 0;
}

/** Search the Commands list for given command.
  * \return the token if found, 0 if not.
  */
DispatchObject::TokenPtr Command::SearchToken(ArgList& argIn) {
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    return (SearchTokenType(DispatchObject::ANALYSIS, argIn));
  }
  // Search for command.
  for (DispatchObject::TokenPtr token = Commands; 
                                token->Type != DispatchObject::NONE; ++token)
    if (argIn.CommandIs( token->Cmd )) return token;
  mprintf("[%s]: Command not found.\n", argIn.Command());
  return 0;
}
