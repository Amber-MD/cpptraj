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
    enum modesType {
      MT_UNKNOWN=0,   MT_DIST,   MT_COVAR, MT_MWCOVAR,
      MT_DISTCOVAR, MT_CORREL, MT_IDEA,  MT_IRED
    };
    enum modesSource {
      MS_UNKNOWN=0, MS_STACK, MS_FILE
    };

    ModesInfo();
    ~ModesInfo();
    int ReadEvecFile(std::string&, int, int);

    int Nvect() { return nvect_; }

    double Evec(int veci, int npair) { 
      return evec_[veci * nvectelem_ + npair]; 
    }

  private:
    static const size_t BUFSIZE_;
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
