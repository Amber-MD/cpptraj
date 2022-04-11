#ifndef INC_EXEC_ADDMISSINGRES_PRES_H
#define INC_EXEC_ADDMISSINGRES_PRES_H
#include <string>
/// Placeholder for residues
class Exec_AddMissingRes::Pres {
  public:
    Pres() {}
    /// CONSTRUCTOR - name, original num, icode, chain
    Pres(std::string const& n, int o, char i, char c) :
      name_(n), oNum_(o), tNum_(-1), icode_(i), chain_(c) {}
    /// CONSTRUCTOR - name, original num, top num, icode, chain
    Pres(std::string const& n, int o, int t, char i, char c) :
      name_(n), oNum_(o), tNum_(t), icode_(i), chain_(c) {}

    std::string const& Name() const { return name_; }
    int Onum()                const { return oNum_; }
    int Tnum()                const { return tNum_; }
    char Icode()              const { return icode_; }
    char Chain()              const { return chain_; }
  private:
    std::string name_; ///< Residue name
    int oNum_;         ///< Original (e.g. PDB) residue number
    int tNum_;         ///< Topology (start from 0) residue index
    char icode_;       ///< Residue insertion code
    char chain_;       ///< Residue chain
};
#endif
