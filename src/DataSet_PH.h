#ifndef INC_DATASET_PH
#define INC_DATASET_PH
#include "DataSet.h"
#include "NameType.h"
/// Hold data from constant pH simulations; protonation states of each residue.
class DataSet_PH : public DataSet {
    typedef std::vector<int> Iarray;
    typedef std::vector<bool> Barray;
    typedef std::vector<float> Farray;
  public:
    DataSet_PH();
    static DataSet* Alloc() { return (DataSet*)new DataSet_PH(); }
    // -------------------------------------------
    class Residue {
      public:
        Residue() : resid_(-1) {}
        /// Res name, num, protcnts, max prot
        Residue(NameType const&, int, Iarray const&, int);
        Residue(int nframes) : states_(nframes, 0) {}
        /// COPY
        Residue(Residue const&);
        /// ASSIGN
        Residue& operator=(Residue const&);
        int& operator[](int n)          { return states_[n];        }
        int operator[](int n) const     { return states_[n];        }
        void push_back(int state)       { states_.push_back(state); }
        unsigned int Nframes()    const { return states_.size();    }
        typedef Iarray::const_iterator const_iterator;
        const_iterator begin()    const { return states_.begin(); }
        const_iterator end()      const { return states_.end();   }
        void Print() const;
        int State(unsigned int n) const { return states_[n];       }
        void Allocate(unsigned int n)   { states_.reserve( n );    }
        void Resize(unsigned int n)     { states_.assign( n, 0 );  }
#       ifdef MPI
        Iarray const& States()    const { return states_; }
        int* StatesPtr()                { return &states_[0]; }
#       endif
        bool IsProtonated(int s)  const { return protonated_[s]; }
        int Nprotons(int s)       const { return protcnts_[s];   }
        NameType const& Name()    const { return resname_;       }
        int Num()                 const { return resid_;         }
      private:
        NameType resname_;  ///< Residue name.
        int resid_;         ///< Residue number.
        Iarray protcnts_;   ///< Hold protonation count for each state.
        Barray protonated_; ///< True if state is protonated.
        Iarray states_;     ///< Hold protonation state for each frame.
    };
    typedef std::vector<Residue> Rarray;
    typedef Rarray::const_iterator const_iterator;
    // -------------------------------------------
    // Residue functions
    void SetResidueInfo(Rarray const& r) { residues_ = r; }
    Rarray const& Residues()       const { return residues_; }
    const_iterator begin()         const { return residues_.begin(); }
    const_iterator end()           const { return residues_.end();   }
    Residue const& Res(int idx)    const { return residues_[idx]; }

    typedef Farray::const_iterator ph_iterator;
    Farray const& pH_Values() const { return solvent_pH_; }

    // NOTE: Bounds check should be done outside of here.
    void AddState(Iarray const& resStates, float pH) {
      for (unsigned int res = 0; res < resStates.size(); res++)
        residues_[res].push_back( resStates[res] );
      solvent_pH_.push_back( pH );
    }
    /// Set pH and state for specified frame/residue
    void SetState(unsigned int res, int frame, int state, float pH) {
      residues_[res][frame] = state;
      solvent_pH_[frame] = pH; // TODO should pass in an array for res?
    }
    /// Resize pH array and state array for each residue.
    void Resize(size_t);
    /// \return number of frames
    unsigned int Nframes() const {
      if (residues_.empty())
        return 0;
      else
        return residues_.front().Nframes();
    }
    // ----- DataSet functions -------------------
    size_t Size() const { return residues_.size(); }
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
  private:
    Farray solvent_pH_;    ///< Solvent pH values each frame
    Rarray residues_;      ///< Array of residues.
};
#endif
