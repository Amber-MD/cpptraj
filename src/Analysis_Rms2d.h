#ifndef INC_ANALYSIS_RMS2D_H
#define INC_ANALYSIS_RMS2D_H
#include "Analysis.h"
#include "Trajin_Single.h"
#include "DataSet_Coords.h"
#include "DataSet_MatrixFlt.h"
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
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    int CalcRmsToTraj();
    int Calculate_2D();
    void CalcAutoCorr();

    static const char* ModeStrings_[];
    enum ModeType { RMS_FIT = 0, RMS_NOFIT, DME };
    ModeType mode_;
    DataSet_Coords* coords_;   ///< Hold coords from input frames.
    bool useReferenceTraj_;    ///< If true read from reference trajectory.
    bool useMass_;             ///< If true, mass-weight.
    AtomMask TgtMask_;         ///< Target atom mask.
    AtomMask RefMask_;         ///< Reference atom mask.
    Trajin_Single RefTraj_;    ///< Reference trajectory, each frame used in turn.
    Topology* RefParm_;        ///< Reference trajectory Parm.
    DataSet_MatrixFlt* rmsdataset_;
    DataSet* Ct_;
};
#endif
