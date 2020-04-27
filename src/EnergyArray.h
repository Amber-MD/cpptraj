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
  private:
    static const char* TypeStr_[];

    typedef std::vector<double> Darray;
    typedef std::vector<Type> Tarray;

    Darray ene_;         ///< Energy array
    Tarray activeTerms_; ///< Active energy terms
    
};
#endif
