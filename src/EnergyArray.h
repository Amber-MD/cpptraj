#ifndef INC_ENERGYARRAY_H
#define INC_ENERGYARRAY_H
#include <vector>
// Forward declares
class CpptrajFile;
/// Hold energy terms
class EnergyArray {
  public:
    EnergyArray();
    /// Energy term types. Keep in sync with TypeStr_.
    enum Type { E_BOND = 0, E_ANGLE, E_DIHEDRAL, E_VDW14, E_Q14, E_VDW, E_COULOMB, E_OPENMM, N_E_TERMS };
    /// \return Pointer to specified part of the energy array.
    double* AddType(Type);
    /// Clear all terms.
    void clear() { activeTerms_.clear(); }
    /// Zero all active terms.
    void zero() {
      for (Tarray::const_iterator it = activeTerms_.begin(); it != activeTerms_.end(); ++it)
        ene_[*it] = 0.0;
    }
    /// \return Total of all active terms
    double Total() const {
      double etotal = 0.0;
      for (Tarray::const_iterator it = activeTerms_.begin(); it != activeTerms_.end(); ++it)
        etotal += ene_[*it];
      return etotal;
    }
    /// \return Specified term
    double Ene(Type t) const { return ene_[(int)t]; }
    /// \return Label for specified term
    const char* label(Type) const;
    /// Print active terms to specified file (optionally including single terms).
    void PrintActiveTerms(CpptrajFile&, bool) const;
    /// Print active labels to specified file (optionally including single terms).
    void PrintActiveLabels(CpptrajFile&, bool) const;
    /// Print active term labels to stdout
    void PrintActiveTerms() const;
  private:
    static const char* TypeStr_[];

    typedef std::vector<double> Darray;
    typedef std::vector<Type> Tarray;

    Darray ene_;         ///< Energy array for all terms in Type
    Tarray activeTerms_; ///< Energy terms in ene_ that are currently active
    
};
#endif
