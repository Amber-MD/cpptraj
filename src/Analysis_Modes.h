#ifndef INC_ANALYSIS_MODES_H
#define INC_ANALYSIS_MODES_H
#include "Analysis.h"
#include "DataSet_Modes.h"
#include "Trajout_Single.h"
class Analysis_Modes : public Analysis {
  public:
    Analysis_Modes();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Modes(); }
    void Help() const;

    ~Analysis_Modes();

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    static const double CONSQ;
    static const double TKBC2;
    static const double AVO;
    static const double CNST;
    static const double CMTOA;
    static const double CONT;

    enum modeAnalysisType { FLUCT=0, DISPLACE, CORR, TRAJ, EIGENVAL, RMSIP };
    static const char* analysisTypeString[];
    typedef std::vector< std::pair<int,int> > modestackType;
    typedef modestackType::const_iterator modestack_it;

    void CheckDeprecated(ArgList&,std::string&, const char*);
    /// Calc rms atomic fluctuations for modes
    void CalcFluct(DataSet_Modes const&);
    /// Calc displacement of coordinates along normal mode directions
    void CalcDisplacement(DataSet_Modes const&);
    /// Calc dipole-dipole correlation functions
    void CalcDipoleCorr(DataSet_Modes const&);
    /// Create pseudo-traj along mode
    int ProjectCoords(DataSet_Modes const&);
    /// DEBUG: calculate projection of coordintaes along eigenvector
    void CalculateProjection(int,Frame const&,int) const;
    /// Calc Eigenvalue fraction
    void CalcEvalFrac(DataSet_Modes const&);
    /// Calc Root mean square inner product
    int CalcRMSIP(DataSet_Modes const&, DataSet_Modes const&);

    int debug_;
    modeAnalysisType type_;
    int beg_;
    int end_;
    bool bose_;
    bool calcAll_;
    double factor_;
    DataSet_Modes* modinfo_;
    DataSet_Modes* modinfo2_;
    CpptrajFile* outfile_;
    modestackType atompairStack_;
    Topology* tOutParm_;
    Trajout_Single trajout_;
    int tMode_;
    double pcmin_;
    double pcmax_;
};
#endif
