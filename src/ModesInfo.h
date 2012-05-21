#ifndef INC_MODESINFO_H
#define INC_MODESINFO_H
#include "DataSet.h" 
/** Information relating to modes. This is a DataSet so that it can be
  * added to the master DataSetList by matrix analysis for referencing 
  * by other analyses.
  *
  *  eigenvectors are stored in evec as:
  *  [evec(1,1,x),evec(1,1,y),evec(1,1,z), ..., evec(1,n,x),
  *   ...,
  *   evec(n,1,x), ...,                         evec(n,n,x)]
  */
class ModesInfo : public DataSet {
  public:
    // NOTE: THIS MUST BE IN SAME ORDER AS MatrixType::matrixMode. 
    //       Is it even necessary?
    enum modesType {
      MT_UNKNOWN=0, MT_DIST,      MT_COVAR, MT_MWCOVAR,
      MT_CORREL,    MT_DISTCOVAR, MT_IDEA,  MT_IRED
    };
    enum modesSource {
      MS_UNKNOWN=0, MS_STACK, MS_FILE
    };

    ModesInfo();
    ModesInfo(modesType,modesSource,std::string&);
    ~ModesInfo();

    int ReadEvecFile(std::string&, int, int);
    double calc_spectral_density(double *, int, double);
    double* CalcRMSfluct(int, int, bool);

    // NOTE: Replace all these with a constructor eventually?
    int SetNavgElem(int);
    void SetAvg( double* avgIn ) { avg_ = avgIn; }
    void SetNvect( int nvIn )    { nvect_ = nvIn; }
    void SetNvectElem( int nveIn ) { nvectelem_ = nveIn; }
    void SetFreq( double* fIn )  { freq_ = fIn; }
    void SetEvec( double* eIn )  { evec_ = eIn; }

    int Nvect()    { return nvect_; }
    int NvectElem() { return nvectelem_; }
    int Navgelem() { return navgelem_; }
    double* Freq() { return freq_; }
    double* Evec() { return evec_; }
    modesSource Source() { return source_; }
    modesType Mtype() { return type_; }

    double Evec(int veci, int npair) { 
      return evec_[veci * nvectelem_ + npair]; 
    }

  private:
    static const size_t BUFSIZE_;
    /// hc/2kT in cm, with T=300K; use for quantum Bose statistics)
    static const double CONSQ;
    static const double TKBC2;
    static const double AVO;
    static const double CNST;
    static const double CMTOA;
    static const double TWOPI;
    static const double CONT; 

    modesType type_;
    modesSource source_;
    int navgelem_;
    double *avg_;
    int nvect_;
    int nvectelem_;
    double *freq_;
    double *evec_;
};
#endif
