#ifndef INC_TRAJ_CHARMMCOR_H
#define INC_TRAJ_CHARMMCOR_H
#include "TrajectoryIO.h"
#include "CpptrajFile.h"
/// Read CHARMM Cor file
class Traj_CharmmCor : public TrajectoryIO {
  public:
    Traj_CharmmCor() : corAtom_(0), extendedFmt_(false) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_CharmmCor(); }
    //static void ReadHelp();
    static void WriteHelp();
    /** CORWRITEMODE: Indicate how the COR should be written.
      *  SINGLE: Writing only a single frame.
      *  MULTI: Each frame written to a different file with name filename.frame
      */
    enum CORWRITEMODE {NONE = 0, SINGLE, MULTI};
  private:
    // TrajectoryIO functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&, DataSetList const&);
    int readVelocity(int, Frame&)  { return 1; }
    int readForce(int, Frame&)     { return 1; }
    int processReadArgs(ArgList&)  { return 0; }

    typedef std::vector<std::string> Sarray;

    static const char* EXTENDED_FORMAT_;
    static const char* REGULAR_FORMAT_;

    CpptrajFile file_;
    int corAtom_;               ///< # of atoms in Cor file.
    bool extendedFmt_;          ///< True for wide columns
    CORWRITEMODE corWriteMode_; ///< Determine if writing single or multiple files
    Topology* corTop_;          ///< Corresponding topology
    Sarray MaskSegPairs_;       ///< Hold masks and user-specified segment IDs.
    Sarray SegmentIds_;         ///< Hold segment ID of each residue
    const char* outputFmt_;     ///< Hold output format string
    bool prependExt_;           ///< True if prepending file extension with frame #
};
#endif
