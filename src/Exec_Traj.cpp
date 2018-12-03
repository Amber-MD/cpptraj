#include "Exec_Traj.h"
#include "CpptrajStdio.h"

void Exec_Trajin::Help() const {
  mprintf("\t<filename> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t           [%s]\n", DataSetList::TopArgs);
  mprintf("\t           [mdvel <velocities>] [mdfrc <forces>]\n"
          "\t           [as <format keyword>] [ <Format Options> ]\n"
          "\t           [ remdtraj {remdtrajtemp <T> |\n"
          "\t                       remdtrajidx <indices list> |\n"
          "\t                       remdtrajvalues <values list>}\n"
          "\t             [trajnames <rep1>,<rep2>,...,<repN>] ]\n"
          "  Load trajectory specified by <filename> to the input trajectory list.\n"
          "  If desired, additional velocity or force information can be read from\n"
          "  files specified by 'mdvel' and/or 'mdfrc'.\n"
          "  The 'remdtraj' keyword can be used to extract frames for a specific replica\n"
          "  from an ensemble of replica trajectories. In this case, if only <filename> is\n"
          "  specified it is assumed <filename> has format <name>.<ext> where <ext> is\n"
          "  a numerical suffix; other members of the ensemble will be automatically\n"
          "  searched for. Otherwise, additional members can be specified with 'trajnames'.\n"
          "  'remdtrajtemp' can be used to extract frames from 1D T-REMD simulations.\n"
          "  'remdtrajidx' and 'remdtrajvalues' can be used to extract frames from multi-\n"
          "  dimensional REMD simulations; <indices list> and <values list> are comma-\n"
          "  separated lists.\n"
          "  Use 'help Formats trajin' for help with specific formats.\n");
}

// -----------------------------------------------------------------------------
void Exec_Ensemble::Help() const {
  mprintf("\t<file0> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t        [%s]\n", DataSetList::TopArgs);
  mprintf("\t        [trajnames <file1>,<file2>,...,<fileN>]\n"
          "\t        [nosort | [remlog <remlogfile> [nstlim <nstlim> ntwx <ntwx>]]]\n"
          "  Load an ensemble of trajectories starting with <file0> that will be\n"
          "  processed together as an ensemble.\n"
          "  The default behavior is to sort the ensemble by replica. If 'nosort'\n"
          "  is specified the incoming ensemble will not be sorted; this is useful\n"
          "  if the ensemble is already sorted or the trajectories are independent.\n"
          "  If 'remlog' is specified the ensemble will be sorted by coordinate index;\n"
          "  'nstlim' specifies the number of steps between exchanges and 'ntwx'\n"
          "  specifies the number of steps between trajectory writes.\n"
          "  When running in parallel, the 'ensemblesize' command can be used to specify\n"
          "  the number of members in the ensemble, which may improve set-up performance.\n"
          "  Use 'help Formats trajin' for help with specific formats.\n");
}

// -----------------------------------------------------------------------------
void Exec_Reference::Help() const {
  mprintf("\t<name> [<frame#>] [<mask>] [{[TAG] | name <setname>}] [lastframe] [crdset]\n"
          "\t       [%s]\n", DataSetList::TopArgs);
  mprintf("  Load trajectory file <name> as a reference frame.\n"
          "  If 'crdset' is specified use COORDS data set specified by <name> as reference.\n"
          "  Use 'help Formats trajin' for help with specific formats.\n");
}

// -----------------------------------------------------------------------------
void Exec_Trajout::Help() const {
  mprintf("\t<filename> [<fileformat>] [append] [nobox] [novelocity]\n"
          "\t           [notemperature] [notime] [noforce] [noreplicadim]\n"
          "\t           [%s] [onlyframes <range>] [title <title>]\n"
          "\t           [onlymembers <memberlist>]\n", DataSetList::TopArgs);
  mprintf("\t           %s\n", ActionFrameCounter::HelpText);
  mprintf("\t           [ <Format Options> ]\n"
          "  Write frames after all actions have been processed to output trajectory\n"
          "  specified by <filename>.\n"
          "  Use 'help Formats trajout' for help with specific formats.\n");
}

// -----------------------------------------------------------------------------
void Exec_EnsembleSize::Help() const {
  mprintf("\t<#>\n  Set expected ensemble size to <#> in order to improve ensemble setup efficiency in parallel.\n");
}

Exec::RetType Exec_EnsembleSize::Execute(CpptrajState& State, ArgList& argIn) {
  int eSize = argIn.getNextInteger(-1);
  if (eSize > 0)
# ifdef MPI
    mprintf("\tSetting expected ensemble size to %i\n", eSize);
  if (Parallel::SetupComms( eSize )) return CpptrajState::ERR;
# else
    mprintf("Warning: This command has no effect when not running in parallel.\n");
# endif
  return CpptrajState::OK;
}
