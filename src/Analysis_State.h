#ifndef INC_ANALYSIS_STATE_H
#define INC_ANALYSIS_STATE_H
#include <map>
#include "Analysis.h"
#include "DataSet_double.h"
/// Analyze transitions between states
class Analysis_State : public Analysis {
  public:
    Analysis_State() : state_data_(0), curveOut_(0), stateOut_(0), transOut_(0), debug_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_State(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    class StateType;
    class Transition;
    typedef std::vector<StateType> StateArray;
    typedef std::pair<int,int> StatePair;
    typedef std::map<StatePair, Transition> TransMapType;
    typedef std::pair<StatePair, Transition> TransPair;
    typedef std::vector<int> Iarray;

    std::string const& StateName(int) const;
    const char* stateName(int i) const { return StateName(i).c_str(); }

    static std::string UNDEFINED_;
    StateArray States_;
    TransMapType TransitionMap_;
    DataSet* state_data_;
    DataSetList* masterDSL_;
    DataFile* curveOut_;
    CpptrajFile* stateOut_;
    CpptrajFile* transOut_;
    int debug_;
    bool normalize_;
};
// ----- PRIVATE CLASS DEFINITIONS ---------------------------------------------
/// Hold the definition of a state and associated data
class Analysis_State::StateType {
  public:
    StateType() : set_(0), min_(0.0), max_(0.0) {}
    StateType(std::string const& i, DataSet_1D* d, double m, double x) :
              id_(i), set_(d), min_(m), max_(x) {}
    std::string const& ID() const { return id_; }
    const char* id()        const { return id_.c_str(); }
    DataSet_1D const& DS()  const { return *set_; }
    double Min()            const { return min_; }
    double Max()            const { return max_; }
  private:
    std::string id_;
    DataSet_1D* set_;
    double min_;
    double max_;
};
/// Hold information about a transition from state0 to state1
class Analysis_State::Transition {
  public:
    Transition() : maxLifetime_(0), sumLifetimes_(0), Nlifetimes_(0), curve_(0) {}
    Transition(DataSet_double* ds) :
       maxLifetime_(0), sumLifetimes_(0), Nlifetimes_(0), curve_(ds) {}
    Transition(int length, DataSet_double* ds) :
       maxLifetime_(length), sumLifetimes_(length), Nlifetimes_(1), curve_(ds)
    {
      DataSet_double& CRV = static_cast<DataSet_double&>( *curve_ );
      CRV.Resize( length );
      for (int lc = 0; lc != length; lc++) CRV[lc]++;
    }
    void Update(int length) {
      if (length > maxLifetime_) maxLifetime_ = length;
      sumLifetimes_ += length;
      ++Nlifetimes_;
      DataSet_double& CRV = static_cast<DataSet_double&>( *curve_ );
      if (length > (int)CRV.Size())
        CRV.Resize( length );
      for (int lc = 0; lc != length; lc++)
        CRV[lc]++;
    }
    int Max() const { return maxLifetime_; }
    int Sum() const { return sumLifetimes_; }
    int Nlifetimes() const { return Nlifetimes_; }
    double Avg() const {
      if (Nlifetimes_ > 0) return (double)sumLifetimes_ / (double)Nlifetimes_;
      return 0.0;
    }
    DataSet_double const& DS() const { return *curve_; }
    void NormalizeCurve() const {
      if (curve_->Size() > 0) {
        DataSet_double& CRV = static_cast<DataSet_double&>( *curve_ );
        double norm = 1.0 / CRV[0];
        for (unsigned int i = 0; i != CRV.Size(); i++)
          CRV[i] *= norm;
      }
    }
  private:
    int maxLifetime_; ///< Longest time state0 was active before going to state1
    int sumLifetimes_; ///< Sum of all state0 lifetimes before going to state1
    int Nlifetimes_; ///< Number of state0 lifetimes before going to state1
    DataSet_double* curve_; ///< Lifetime curve for state0->state1
};
#endif
