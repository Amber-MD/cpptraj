#ifndef INC_DATASET_PH
#define INC_DATASET_PH
#include "DataSet_1D.h"
#include "CphResidue.h"
/// Hold data from constant pH simulations; protonation states of a single residue.
class DataSet_pH : public DataSet_1D {
    typedef std::vector<int> Iarray;
  public:
    DataSet_pH();
    static DataSet* Alloc() { return (DataSet*)new DataSet_pH(); }

    // ----- DataSet functions -------------------
    size_t Size() const { return states_.size(); }
    void Info()   const;
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    /// Reserve space for states of each residue
    int Allocate(SizeArray const&);
    void Add( size_t, const void* ) { return; }
    int Append(DataSet*)            { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
    void Consolidate(Parallel::Comm const&, int);
#   endif
    // ----- DataSet_1D functions ----------------
    double Dval(size_t i)           const { return (double)states_[i];         }
    double Xcrd(size_t idx)         const { return Dim(0).Coord(idx);          }
    const void* VoidPtr(size_t idx) const { return (void*)(&(states_[0])+idx); }
    // -------------------------------------------
    void Resize(size_t);
    void SetResidueInfo(CphResidue const& r) { res_ = r; }
    void Set_Solvent_pH( float p )              { solvent_pH_ = p; }
    void SetState(unsigned int n, int s, int r) { states_[n] = s; recType_[n] = r; }
    void AddState(int s, int r)                 { states_.push_back( s ); recType_.push_back( r ); }
    int State(unsigned int idx)       const { return states_[idx];    }
    int RecordType(unsigned int idx)  const { return recType_[idx];      }
    float Solvent_pH()                const { return solvent_pH_;     }
    CphResidue const& Res()           const { return res_;            }
    /// Set Monte Carlo step size, initial time, and time step
    void SetTimeValues(int m, float t0, float dt) { mc_stepsize_ = m; t0_ = t0; dt_ = dt; }
    int MonteCarloStepSize()   const { return mc_stepsize_; }
    float InitialTime()        const { return t0_; }
    float TimeStep()           const { return dt_; }
  private:
    typedef std::vector<float> Farray;

    float solvent_pH_;     ///< Solvent pH
    CphResidue res_;       ///< Hold titratable residue characteristics.
    Iarray states_;        ///< Hold protonation state for each frame.
    Iarray recType_;       ///< Hold record type each frame.
    float t0_;             ///< Initial time
    float dt_;             ///< Time step
    int mc_stepsize_;      ///< Monte Carlo step size
};
#endif
