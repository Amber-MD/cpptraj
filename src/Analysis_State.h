#ifndef INC_ANALYSIS_STATE_H
#define INC_ANALYSIS_STATE_H
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
    typedef std::vector<StateType> StateArray;
    StateArray States_;
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
#endif
