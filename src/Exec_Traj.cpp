#include "Exec_Traj.h"
#include "CpptrajStdio.h"

void Exec_Trajin::Help() const {
  mprintf("\t<filename> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t           [%s]\n", DataSetList::TopArgs);
  mprintf("\t           [ <Format Options> ]\n"
          "\t           [ remdtraj [remdtrajtemp <T> | remdtrajidx <#>]\n"
          "\t             [trajnames <rep1>,<rep2>,...,<repN> ] ]\n"
          "  Load trajectory specified by <filename> to the input trajectory list.\n");
  TrajectoryFile::ReadOptions();
}
// -----------------------------------------------------------------------------
void Exec_Ensemble::Help() const {
  mprintf("\t<file0> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t        [%s]\n", DataSetList::TopArgs);
  mprintf("\t        [trajnames <file1>,<file2>,...,<fileN>\n"
          "\t        [remlog <remlogfile> [nstlim <nstlim> ntwx <ntwx>]]\n"
          "  Load an ensemble of trajectories starting with <file0> that will be\n"
          "  processed together as an ensemble.\n");
}
// -----------------------------------------------------------------------------
void Exec_Reference::Help() const {
  mprintf("\t<name> [<frame#>] [<mask>] [TAG] [lastframe] [crdset]\n"
          "\t       [%s]\n", DataSetList::TopArgs);
  mprintf("  Load trajectory file <name> as a reference frame.\n"
          "  If 'crdset' is specified use COORDS data set specified by <name> as reference.\n");
}
// -----------------------------------------------------------------------------
void Exec_Trajout::Help() const {
  mprintf("\t<filename> [<fileformat>] [append] [nobox]\n"
          "\t           [%s] [onlyframes <range>] [title <title>]\n", DataSetList::TopArgs);
  mprintf("\t           %s\n", ActionFrameCounter::HelpText);
  mprintf("\t           [ <Format Options> ]\n"
          "  Write frames after all actions have been processed to output trajectory\n"
          "  specified by <filename>.\n");
  TrajectoryFile::WriteOptions();
}
