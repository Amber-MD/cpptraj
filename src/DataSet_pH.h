#ifndef INC_DATASET_PH
#define INC_DATASET_PH
#include "DataSet_1D.h"
#include "CphResidue.h"
/// Hold data from constant pH simulations; protonation states of each residue.
class DataSet_pH : public DataSet_1D {
    typedef std::vector<int> Iarray;
  public:
    DataSet_pH();
    static DataSet* Alloc() { return (DataSet*)new DataSet_pH(); }

    // ----- DataSet functions -------------------
    size_t Size() const { return states_.size(); }
    void Info()   const;
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
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
    void SetResidueInfo(NameType const&, int, Iarray const&, int);
    void SetResidueInfo(DataSet_pH const&);
    void Set_Solvent_pH( float p )       { solvent_pH_ = p; }
    void SetState(unsigned int n, int s) { states_[n] = s; }
    void AddState(int s)          { states_.push_back( s ); }
    void AddState(int s, float p) { states_.push_back( s ); pH_Values_.push_back( p ); }
    bool Has_pH()               const { return !pH_Values_.empty(); }
    int State(unsigned int idx) const { return states_[idx];    }
    float pH(unsigned int idx)  const { return pH_Values_[idx]; }
    float Solvent_pH()          const { return solvent_pH_;     }
    CphResidue const& Res()     const { return res_;            }
  private:
    typedef std::vector<bool> Barray;
    typedef std::vector<float> Farray;

    float solvent_pH_;     ///< Solvent pH
    
    NameType resname_;     ///< Residue name.
    int resid_;            ///< Residue number.
    Iarray protcnts_;      ///< Hold protonation count for each state.
    Barray protonated_;    ///< True if state is protonated.
    Iarray states_;        ///< Hold protonation state for each frame.
};
#endif
