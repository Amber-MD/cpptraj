#ifndef INC_TRAJ_XYZ_H
#define INC_TRAJ_XYZ_H
#include "TrajectoryIO.h"
#include "BufferedLine.h"
/// Read simple XYZ trajectories. 
class Traj_XYZ : public TrajectoryIO {
  public:
    Traj_XYZ();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_XYZ(); }
    static void WriteHelp();
    static void ReadHelp();
  private:
    // ----- Inherited functions -----------------
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int processWriteArgs(ArgList&, DataSetList const&);
    int processReadArgs(ArgList&);
    // -------------------------------------------
#   ifdef MPI
    // ----- Parallel functions ------------------
    int parallelOpenTrajin(Parallel::Comm const&);
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelReadFrame(int, Frame&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj();
    // -------------------------------------------
#   endif
    /// Data line types
    enum Type { UNKNOWN=0,
                XYZ,       ///< <X> <Y> <Z>
                ATOM_XYZ,  ///< <#> <X> <Y> <Z>
                NAME_XYZ,  ///< <Name> <X> <Y> <Z>
              };
    /// Specify how/when headers should appear
    enum TitleType { NO_TITLE = 0,  ///< Never
                     SINGLE,        ///< First frame only
                     MULTIPLE,      ///< Every frame
                     NATOM_COMMENT, ///< Every frame has # atoms followed by a comment
                     UNKNOWN_TITLE
                   };
    /// Line format types
    enum LineFmtType { SINGLE_INTEGER,
                       THREE_DOUBLES,
                       INTEGER_AND_THREE_DOUBLES,
                       STRING_AND_THREE_DOUBLES,
                       UNKNOWN_LINE_FORMAT };
    /// \return Format of given line
    static LineFmtType DetermineLineFormat(std::string const&);
    /// \return File format based on first three lines
    static Type DetermineFormat(std::string&, std::string const&, std::string const&,
                                std::string const&);
    /// Set given box from box line
    static int parseBoxLine(Box&, ArgList const&);
    /// Read past header, process box info
    inline void ReadTitle(Box&);
    /// Read past header
    inline void ReadTitle();
    /// Read xyz coords
    int readXYZ(int, int, double*, Box&);

    static const char* FMT_XYZ_;
    static const char* FMT_ATOM_XYZ_;
    static const char* FMT_NAME_XYZ_;

    BufferedLine file_;
    std::string ofmt_;
    TitleType titleType_;
    Type ftype_;
    int set_;
    int width_;
    int prec_;
    const char* fmt_; ///< Format for reading
    bool hasBox_;     ///< True if comment line contains box information
};
#endif
