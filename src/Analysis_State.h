#ifndef INC_ANALYSIS_STATE_H
#define INC_ANALYSIS_STATE_H
#include <map>
#include "Analysis.h"
#include "Array1D.h"
#include "DataSet_double.h"
/// Analyze transitions between states
class Analysis_State : public Analysis {
  public:
    Analysis_State();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_State(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    // ----- PRIVATE CLASS DEFINITIONS ---------------------------------------------
    /// Hold the definition of a state and associated data
    class StateType {
      public:
        StateType() {}
        StateType(std::string const& i) : id_(i) {}
        /// Add criterion for state
        void AddCriterion(DataSet_1D* d, double m, double x) {
          Sets_.push_back( d );
          Min_.push_back( m );
          Max_.push_back( x );
        }
        /// \return smallest number of frames among all criterion DataSets
        size_t Nframes() const;
        /// \return State ID
        std::string const& ID() const { return id_; }
        /// \return State ID as const char*
        const char* id()        const { return id_.c_str(); }
        /// \return true if given frame satifies all criteria
        bool InState(int n) const {
          for (unsigned int idx = 0; idx != Sets_.size(); idx++) {
            double dval = Sets_[idx]->Dval(n);
            if (dval < Min_[idx] || dval >= Max_[idx]) // TODO periodic
              return false;
          }
          return true;
        }
        /// Print state info to stdout
        void PrintState() const;
      private:
        typedef std::vector<double> Darray;
        std::string id_;
        Array1D Sets_;   ///< DataSets used to determine if we are in State
        Darray Min_;    ///< Above this value we are in state.
        Darray Max_;    ///< Below this value we are in state.
    };
    /// Hold information about a transition from state0 to state1
    class Transition {
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
    CpptrajFile* countOut_;
    int debug_;
    bool normalize_;
};

#endif
