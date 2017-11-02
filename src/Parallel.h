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
# include <sys/types.h> // off_t FIXME necessary?
# ifdef PARALLEL_DEBUG_VERBOSE
#   include <cstdio> // for FILE
# endif
#endif
/// Static class, Cpptraj C++ interface to C MPI routines.
/** NOTE: I have decided to use a class instead of a namespace to have the
  *       option of making things private. Not sure if that is really
  *       necessary.
  * Send/Receive Tags: All Comm Send()/Recv() tags should be noted here in order.
  *     1212  : Frame::X_
  *     1213  : Frame::box_
  *     1214  : Frame::T_
  *     1215  : Frame::V_
  *     1216  : Frame::remd_indices_
  *     1217  : Frame::time_
  *     1218  : Frame::F_
  *     1234  : Used by Parallel::SendMaster()
  *     1300  : Action_Hbond: Array containing hbond double info on rank.
  *     1301  :   Array containing hbond integer info on rank.
  *     1302  :   Number of bridges to expect from rank.
  *     1303  :   Array containing bridge integer info on rank.
  *     1304+X:   Array of hbond X series info from rank.
  *     1400+X: Action_AtomicCorr: Atomic movement vectors
  *     1500  : Action_NAstruct: Array containing BP step info on rank.
  *     1501+X:   Array of step series data from rank.
  */
class Parallel {
  public:
    /// C++ class wrapper around MPI comm routines.
    class Comm;
    /// C++ class wrapper around MPI file routines.
    class File;
    /// \return MPI_WORLD_COMM Comm
    static Comm const& World() { return world_; }
    /// Initialize parallel environment
    static int Init(int, char**);
    /// Stop parallel environment
    static int End();
#   ifdef CPPTRAJ_MPI
    static int Abort(int);
    /// Set up ensemble and trajectory communicators for given ensemble size
    static int SetupComms(int, bool);
    /// Set up communicators - do not allow fewer threads than groups. TODO remove
    static int SetupComms(int n) { return SetupComms(n, false); }
    /// For DEBUG: infinite loop, gives time to attach a debugger.
    static void Lock();
    /// \return Across-ensemble communicator; includes trajectory comm. masters.
    static Comm const& EnsembleComm()   { return ensembleComm_;   }
    /// \return total ensemble size.
    static int Ensemble_Size()          { return ensemble_size_;  }
    /// \return First ensemble member this thread is responsible for.
    static int Ensemble_Beg()           { return ensemble_beg_;   }
    /// \return Last+1 ensemble member this thread is responsible for.
    static int Ensemble_End()           { return ensemble_end_;   }
    /// \return Total number of ensemble members this thread is responsible for.
    static int N_Ens_Members()          { return n_ens_members_;  }
    /// \return Rank in ensemble comm. for given member.
    static int MemberEnsCommRank(int i) { return memberEnsRank_[i];}
    /// \return Across-trajectory communicator.
    static Comm const& TrajComm()       { return trajComm_;       }
    /// \return trajectory comm. if active, world comm. otherwise.
    static Comm const& ActiveComm();
#   ifdef PARALLEL_DEBUG_VERBOSE
    static FILE* mpidebugfile_;
#   endif /* PARALLEL_DEBUG_VERBOSE */
#   endif /* CPPTRAJ_MPI */
  private:
#   ifdef CPPTRAJ_MPI
    static void printMPIerr(int, const char*, int);
    static int checkMPIerr(int, const char*, int);
    static int ensemble_size_;  ///< Total number of ensemble members.
    static int ensemble_beg_;   ///< Starting member for this ensemble thread.
    static int ensemble_end_;   ///< Ending member for this ensemble thread.
    static int n_ens_members_;  ///< Number of ensemble members thread is responsible for.
    static int* memberEnsRank_; ///< Rank in ensemble comm for each member.
#   ifdef PARALLEL_DEBUG_VERBOSE
    static void dbgprintf(const char*, ...);
    static int debug_init();
    static int debug_end();
#   endif /* PARALLEL_DEBUG_VERBOSE */
    static Comm ensembleComm_;   ///< Communicator across ensemble.
    static Comm trajComm_;       ///< Communicator across single trajectory.
#   endif /* CPPTRAJ_MPI */
    static Comm world_;          ///< World communicator.
};
/// Wrapper around MPI communicator.
class Parallel::Comm {
  public:
    int Rank()    const { return rank_;      }
    int Size()    const { return size_;      }
    bool Master() const { return rank_ == 0; }
#   ifdef CPPTRAJ_MPI
    Comm() : comm_(MPI_COMM_NULL), rank_(0), size_(1) {}
    Comm(MPI_Comm);
    //~Comm();
    Comm(Comm const&);
    Comm& operator=(Comm const&);
    /// \return Internal MPI_Comm
    MPI_Comm MPIcomm() const { return comm_; }
    bool IsNull() const { return comm_ == MPI_COMM_NULL; }
    void Barrier() const;
    /// Split this Comm into a new Comm, give current rank the given ID
    Comm Split(int) const;
    void Reset();
    /// my_start, my_stop, maxElts
    int DivideAmongThreads(int&, int&, int) const;
    /// RecvBuffer, SendBuffer, Count, DataType, Op
    int ReduceMaster(void*, void*, int, MPI_Datatype, MPI_Op) const;
    /// Rank, RecvBuffer, SendBuffer, Count, DataType, Op
    int Reduce(int, void*, void*, int, MPI_Datatype, MPI_Op) const;
    /// Buffer, Count, Rank, DataType 
    int SendMaster(void*, int, int, MPI_Datatype) const;
    /// Return, Input, Count, DataType, Op
    int AllReduce(void*, void*, int, MPI_Datatype, MPI_Op) const;
    /// SendBuffer, Count, DataType, RecvBuffer
    int GatherMaster(void*, int, MPI_Datatype, void*) const;
    /// SendBuffer, Count, DataType, RecvBuffer
    int AllGather(void*, int, MPI_Datatype, void*) const;
    /// Buffer, Count, DataType, Destination Rank, Tag
    int Send(void*, int, MPI_Datatype, int, int) const;
    /// Buffer, Count, DataType, Source Rank, Tag
    int Recv(void*, int, MPI_Datatype, int, int) const;
    /// Buffer, Count, DataType
    int MasterBcast(void*, int, MPI_Datatype) const;
    int CheckError(int) const;
#   else
    Comm() : rank_(0), size_(1) {}
#   endif
  private:
#   ifdef CPPTRAJ_MPI
    MPI_Comm comm_;
#   endif
    int rank_;
    int size_;
};
/// Wrapper around MPI file routines.
class Parallel::File {
  public:
    File() {}
#   ifdef CPPTRAJ_MPI
    int OpenFile(const char*, const char*, Comm const&);
    int Flush();
    off_t Position();
    int CloseFile();
    int Fread(void*, int, MPI_Datatype);
    int Fwrite(const void*, int, MPI_Datatype);
    int Fseek(off_t, int);
    char* Fgets(char*, int);
    int SetSize(long int);
#   endif
  private:
#   ifdef CPPTRAJ_MPI
    MPI_File file_;
    Comm comm_;
#   endif
};
// Restore MPI definition
#ifdef CPPTRAJ_MPI
# undef CPPTRAJ_MPI
# define MPI
#endif
#endif
