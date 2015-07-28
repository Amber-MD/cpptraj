#ifndef INC_METADATA_H
#define INC_METADATA_H
#include <string>
class MetaData {
  public:
    MetaData() : ensembleNum_(-1), idx_(-1), scalarmode_(UNKNOWN_MODE), scalartype_(UNDEFINED) {}
    /// Source of data stored in DataSet, used by Analysis_Statistics. Must match Smodes.
    enum scalarMode {
      M_DISTANCE=0, M_ANGLE, M_TORSION, M_PUCKER, M_RMS, M_MATRIX, UNKNOWN_MODE
    };
    /// Type of DataSet, used by Analysis_Statistics. Must match Stypes, TypeModes
    enum scalarType {
      ALPHA=0, BETA, GAMMA, DELTA, EPSILON, ZETA,  PUCKER,
      CHI,     H1P,  C2P,   PHI,   PSI,     PCHI,  OMEGA,
      NOE,
      DIST,   COVAR,     MWCOVAR, //FIXME: May need more descriptive names 
      CORREL, DISTCOVAR, IDEA,
      IRED,   DIHCOVAR,
      UNDEFINED
    };
    /// Comparison for sorting name/aspect/idx
    inline bool operator<(const MetaData&) const;
    /// \return string containing scalar mode and type if defined.
    std::string ScalarDescription() const;
    /// \return scalarMode that matches input keyword.
    static scalarMode ModeFromKeyword(std::string const&);
    /// \return scalarType that matches keyword; check that mode is valid if specified.
    static scalarType TypeFromKeyword(std::string const&, scalarMode&);
    static scalarType TypeFromKeyword(std::string const&, scalarMode const&);
    const char* ModeString()    const { return Smodes[scalarmode_];}
    const char* TypeString()    const { return Stypes[scalartype_];}
    static const char* ModeString(scalarMode m) { return Smodes[m]; }
    static const char* TypeString(scalarType t) { return Stypes[t]; }
    /// Set default legend from name/aspect/index/ensemble num
    void SetDefaultLegend();
    /// \return string containing name based on metadata
    std::string PrintName() const;
    /// \return true if given MetaData matches this exactly.
    bool Match_Exact(MetaData const&) const;

    std::string const& Name()   const { return name_;        }
    std::string const& Aspect() const { return aspect_;      }
    std::string const& Legend() const { return legend_;      }
    int Idx()                   const { return idx_;         }
    int EnsembleNum()           const { return ensembleNum_; }
  private:
    static const char* Smodes[]; ///< String for each scalar mode
    static const char* Stypes[]; ///< String for each scalar type
    static const scalarMode TypeModes[]; ///< The mode each type is associated with.
    std::string name_;        ///< Name of the DataSet
    std::string aspect_;      ///< DataSet aspect.
    std::string legend_;      ///< DataSet legend.
    int ensembleNum_;         ///< DataSet ensemble number.
    int idx_;                 ///< DataSet index
    scalarMode scalarmode_;   ///< Source of data in DataSet.
    scalarType scalartype_;   ///< Specific type of data in DataSet (if any).
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
bool MetaData::operator<(const MetaData& rhs) const {
  if ( name_ == rhs.name_ ) {
    if ( aspect_ == rhs.aspect_ ) {
      return ( idx_ < rhs.idx_ );
    } else {
      return ( aspect_ < rhs.aspect_ );
    }
  } else {
    return ( name_ < rhs.name_ );
  }
}

#endif
