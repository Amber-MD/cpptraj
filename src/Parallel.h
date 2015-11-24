#ifndef INC_PARALLEL_H
#define INC_PARALLEL_H
/** In C++, MPI is a reserved namespace, but since Amber uses it to define
  * parallel build we need to work with it. Change MPI to CPPTRAJ_MPI for
  * this header file only.
  */
#ifdef MPI
# undef MPI
# define CPPTRAJ_MPI
# include <mpi.h>
#endif
/// Static class, Cpptraj C++ interface to C MPI routines.
class Parallel {
  public:
    /// C++ class wrapper around MPI comm routines.
    class Comm;
    /// \return MPI_WORLD_COMM Comm
    static Comm const& World() { return world_; }
    /// Initialize parallel environment
    static int Init(int, char**);
    /// Stop parallel environment
    static int End();
  private:
#   ifdef CPPTRAJ_MPI
    static void printMPIerr(int, const char*, int);
    static int checkMPIerr(int, const char*, int);
#   endif
    static Comm world_;
};

class Parallel::Comm {
  public:
    Comm() : rank_(0), size_(1) {}
    int Rank() const { return rank_; }
    int Size() const { return size_; }
#   ifdef CPPTRAJ_MPI
    Comm(MPI_Comm);
#   endif
  private:
#   ifdef CPPTRAJ_MPI
    MPI_Comm comm_;
#   endif
    int rank_;
    int size_;
};
#ifdef CPPTRAJ_MPI
# undef CPPTRAJ_MPI
# define MPI
#endif
#endif
