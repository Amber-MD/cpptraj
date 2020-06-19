#ifndef INC_ENERGYARRAY_H
#define INC_ENERGYARRAY_H
#include <vector>
class EnergyArray {
  public:
    EnergyArray();
    /// Energy term types
    enum Type { E_BOND = 0, N_E_TERMS };
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
  private:
    static const char* TypeStr_[];

    typedef std::vector<double> Darray;
    typedef std::vector<Type> Tarray;

    Darray ene_;         ///< Energy array
    Tarray activeTerms_; ///< Active energy terms
    
};
#endif
