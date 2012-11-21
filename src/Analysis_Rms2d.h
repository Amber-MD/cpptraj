#ifndef INC_ANALYSIS_RMS2D_H
#define INC_ANALYSIS_RMS2D_H
#include "Analysis.h"
#include "Trajin_Single.h"
#include "DataSet_Coords.h"
#include "TriangleMatrix.h"
// Class: Analysis_Rms2d
/// Calculate the RMSD between two sets of frames.
/** Perform RMS calculation between each input frame and each other input 
  * frame, or each frame read in from a separate reference traj and each 
  * input frame. 
  */
class Analysis_Rms2d: public Analysis {
  public:
    Analysis_Rms2d();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Rms2d(); }
    static void Help();
    static int Calc2drms(DataSet_Coords&, TriangleMatrix&,bool,bool,std::string const&);

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,int);
    Analysis::RetType Analyze();
    void Print(DataFileList*);
  private:
    enum ModeType { NORMAL = 0, REFTRAJ };
    ModeType mode_;
    DataSet_Coords* coords_;   ///< Hold coords from input frames.
    bool nofit_;               ///< Do not perform rms fitting
    bool useMass_;
    std::string maskexpr_;
    AtomMask RefMask_;         ///< Reference Traj atom mask
    Trajin_Single RefTraj_;    ///< Reference trajectory, each frame used in turn
    Topology* RefParm_;        ///< Reference trajectory Parm
    DataSet* rmsdataset_;
    DataSet* Ct_;
    std::string rmsdFile_;     ///< Output file
    std::string corrfilename_; ///< Correlation output file.

    int CalcRmsToTraj();
    void CalcAutoCorr(TriangleMatrix&);
};
#endif
