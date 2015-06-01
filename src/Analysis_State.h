#ifndef INC_ANALYSIS_STATE_H
#define INC_ANALYSIS_STATE_H
#include <map>
#include "Analysis.h"
/// Analyze transitions between states
class Analysis_State : public Analysis {
  public:
    Analysis_State() : stateOut_(0) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_State(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    class StateType;
    class Transition;
    typedef std::vector<StateType> StateArray;
    typedef std::pair<int,int> StatePair;
    typedef std::map<StatePair, Transition> TransMapType;
    typedef std::pair<StatePair, Transition> TransPair;
    typedef std::vector<int> Iarray;
    StateArray States_;
    TransMapType TransitionMap_;
    DataSet* stateOut_;
};
// ----- PRIVATE CLASS DEFINITIONS ---------------------------------------------
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
    Transition() : maxLifetime_(0), sumLifetimes_(0), Nlifetimes_(0) {}
    Transition(int length) : maxLifetime_(length), sumLifetimes_(length), Nlifetimes_(1) {}
    void Update(int length) {
      if (length > maxLifetime_) maxLifetime_ = length;
      sumLifetimes_ += length;
      ++Nlifetimes_;
      if (length > (int)curve_.size())
        curve_.resize( length, 0 );
      for (int lc = 0; lc != length; lc++)
        curve_[lc]++;
    }
    int Max() const { return maxLifetime_; }
    int Sum() const { return sumLifetimes_; }
    int Nlifetimes() const { return Nlifetimes_; }
  private:
    int maxLifetime_; ///< Longest time state0 was active before going to state1
    int sumLifetimes_; ///< Sum of all state0 lifetimes before going to state1
    int Nlifetimes_; ///< Number of state0 lifetimes before going to state1
    Iarray curve_; ///< Lifetime curve for state0->state1
};
#endif
