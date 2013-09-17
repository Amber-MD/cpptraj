#ifndef INC_CPPTRAJ_H
#define INC_CPPTRAJ_H
#include "CpptrajState.h"
// Class: Cpptraj:
/// Hold state information.
/** This is the main class for cpptraj. It holds all data and controls the 
 *  overall flow of the program. It exists in main.cpp.
 */
class Cpptraj {
  public:
    Cpptraj() {}
    int RunCpptraj(int, char**);
  private:
    enum Mode { BATCH = 0, ERROR, QUIT, INTERACTIVE, SILENT_EXIT };
    static void Usage();
    static void Intro();
    static void Finalize();
    Mode ProcessCmdLineArgs(int, char**);
    int Interactive();

    CpptrajState State_;
    std::string logfilename_; // TODO: Put in CpptrajState?
};
#endif
