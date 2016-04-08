#include <algorithm> // sort
#include "Exec_CompareTop.h"
#include "CpptrajStdio.h"

void Exec_CompareTop::Help() const {
  mprintf("\t{%s} {%s} [out <file>] [atype] [lj] [bnd] [ang] [dih]\n",
          DataSetList::TopArgs, DataSetList::TopArgs);
}

typedef std::vector<int> Iarray;
typedef std::vector<NameType> Narray;

// -----------------------------------------------------------------------------
/// Function for printing NameType
static void PrintAtype(CpptrajFile& output, Topology const& parm, NameType const& name, char dir)
{
  output.Printf("%c %s\n", dir, *name);
}

/// \return An array of unique atom types
static Narray AtypeArray( Topology const& parm ) {
  std::set<NameType> atypes;
  for (Topology::atom_iterator atom = parm.begin(); atom != parm.end(); ++atom)
    atypes.insert( atom->Type() );
  Narray out;
  for (std::set<NameType>::const_iterator it = atypes.begin(); it != atypes.end(); ++it)
    out.push_back( *it );
  return out;
}

// -----------------------------------------------------------------------------
/// Hold LJ params for an atom
class LJatom {
  public:
    LJatom() : rmin_(0.0), eps_(0.0) {}
    LJatom(NameType const& n, double r, double e) : name_(n), rmin_(r), eps_(e) {}
    NameType name_;
    double rmin_;
    double eps_;
    bool operator<(const LJatom& rhs) const {
      if (name_ == rhs.name_) {
        if (rmin_ == rhs.rmin_) {
          return (eps_ < rhs.eps_);
        } else return (rmin_ < rhs.rmin_);
      } else return (name_ < rhs.name_);
    }
};

/// Array of LJatom
typedef std::vector<LJatom> LJarrayType;

/// \return Array containing LJ params for each atom.
static LJarrayType LJarray( Topology const& parm ) {
  std::set<LJatom> temp;
  for (int idx = 0; idx != parm.Natom(); idx++) {
    NonbondType const& NB = parm.GetLJparam(idx, idx);
    double eps;
    if (NB.A() > 0.0)
      eps = (NB.B() * NB.B()) / (4.0 * NB.A());
    else
      eps = 0.0;
    temp.insert( LJatom( parm[idx].Type(), parm.GetVDWradius(idx), eps ) );
  }
  LJarrayType out;
  for (std::set<LJatom>::const_iterator it = temp.begin(); it != temp.end(); ++it)
    out.push_back( *it );
  return out;
}

/// Function for printing LJatom
static void PrintLJatom(CpptrajFile& output, Topology const& parm,
                        LJatom const& first, char dir)
{
  output.Printf("%c %s Rmin= %g  Eps= %g\n", dir, *(first.name_), first.rmin_, first.eps_);
}

// -----------------------------------------------------------------------------
/// Print direction, atom mask names, and atom types
static inline void PrintHeader(char dir, Iarray const& atoms, Topology const& topIn,
                               CpptrajFile& output)
{
  output.Printf("%c", dir);
  for (unsigned int idx = 0; idx != atoms.size(); idx++) {
    if (idx > 0) output.Printf(" -");
    output.Printf(" %s", topIn.AtomMaskName(atoms[idx]).c_str());
  }
  output.Printf(" {");
  for (unsigned int idx = 0; idx != atoms.size(); idx++) {
    if (idx > 0) output.Printf("-");
    output.Printf("%s", *(topIn[atoms[idx]].Type()));
  }
  output.Printf("}");
}

// -----------------------------------------------------------------------------
/// Used to hold and sort parameters of given parameter type.
template <class T> class ParmT {
  public:
    ParmT() : natom_(0) {}
    ParmT(Iarray const& a, Iarray const& r, Narray const& n, T const& p) :
      atoms_(a), rnums_(r), names_(n), natom_(a.size()), parms_(p) {}
    bool operator<(const ParmT& rhs) const {
    // First sort on residue number
      for (int i = 0; i != natom_; i++)
        if (rnums_[i] < rhs.rnums_[i])
          return true;
        else if (rnums_[i] > rhs.rnums_[i])
          return false;
      // All residue numbers are equal. Sort by atom names.
      for (int i = 0; i != natom_; i++)
        if (names_[i] < rhs.names_[i])
          return true;
        else if (names_[i] > rhs.names_[i])
          return false;
      // All res numbers and atom names are equal. Sort by parameters.
      if (parms_ < rhs.parms_)
        return true;
      else
        return false;
    }
    Iarray const& Atoms() const { return atoms_; }
    T      const& Parms() const { return parms_; }
  private:
    Iarray atoms_; ///< Atom indices involved in parameter
    Iarray rnums_; ///< Residues corresponding to atoms involved in parameter
    Narray names_; ///< Atom names involved in parameter
    int natom_;    ///< Number of atoms involved in parameter
    T parms_;      ///< Parameters
};

// -----------------------------------------------------------------------------
static inline void SetDihParms(Topology const& topIn, DihedralType const& T,
                               Iarray& atoms, Iarray& rnums, Narray& names)
{
  atoms[0] = T.A1();
  atoms[1] = T.A2();
  atoms[2] = T.A3();
  atoms[3] = T.A4();
  rnums[0] = topIn[T.A1()].ResNum();
  rnums[1] = topIn[T.A2()].ResNum();
  rnums[2] = topIn[T.A3()].ResNum();
  rnums[3] = topIn[T.A4()].ResNum();
  names[0] = topIn[T.A1()].Name();
  names[1] = topIn[T.A2()].Name();
  names[2] = topIn[T.A3()].Name();
  names[3] = topIn[T.A4()].Name();
}

/// Hold dihedral parameters
typedef ParmT<DihedralParmType> DihT;

/// Array of dihedral atoms and parameters
typedef std::vector<DihT> DihArrayT;

/// \return Array of dihedral atoms and parameters from given Topology
static DihArrayT DihArray( Topology const& topIn ) {
  DihArrayT array;
  Iarray atoms(4);
  Iarray rnums(4);
  Narray names(4);

  for (DihedralArray::const_iterator it = topIn.Dihedrals().begin();
                                     it != topIn.Dihedrals().end(); ++it)
  {
    SetDihParms( topIn, *it, atoms, rnums, names );
    array.push_back( DihT(atoms, rnums, names, topIn.DihedralParm()[it->Idx()]) );
  }
  for (DihedralArray::const_iterator it = topIn.DihedralsH().begin();
                                     it != topIn.DihedralsH().end(); ++it)
  {
    SetDihParms( topIn, *it, atoms, rnums, names );
    array.push_back( DihT(atoms, rnums, names, topIn.DihedralParm()[it->Idx()]) );
  }
  return array;
}

/// Function for printing DihT
static void PrintDihT(CpptrajFile& output, Topology const& parm,
                      DihT const& DIH, char dir)
{
  PrintHeader(dir, DIH.Atoms(), parm, output);
  DihedralParmType const& DP = DIH.Parms(); 
    output.Printf(" Pk=%g Pn=%g Phase=%g SCEE=%g SCNB=%g\n",
                  DP.Pk(), DP.Pn(), DP.Phase(), DP.SCEE(), DP.SCNB());
}

/// Function for printing DihedralParmType
static void PrintDihP(CpptrajFile& output, Topology const& parm,
                      DihedralParmType const& DP, char dir)
{
  output.Printf("%c Pk= %g  Pn= %g  Phase= %g  SCEE= %g  SCNB= %g\n", dir,
                DP.Pk(), DP.Pn(), DP.Phase(), DP.SCEE(), DP.SCNB());
}

// -----------------------------------------------------------------------------
static inline void SetAngParms(Topology const& topIn, AngleType const& T,
                               Iarray& atoms, Iarray& rnums, Narray& names)
{
  atoms[0] = T.A1();
  atoms[1] = T.A2();
  atoms[2] = T.A3();
  rnums[0] = topIn[T.A1()].ResNum();
  rnums[1] = topIn[T.A2()].ResNum();
  rnums[2] = topIn[T.A3()].ResNum();
  names[0] = topIn[T.A1()].Name();
  names[1] = topIn[T.A2()].Name();
  names[2] = topIn[T.A3()].Name();
}

/// Hold angle parameters
typedef ParmT<AngleParmType> AngT;

/// Array of angle atoms and parameters
typedef std::vector<AngT> AngArrayT;

/// \return Array of angle atoms and parameters from given Topology
static AngArrayT AngArray( Topology const& topIn ) {
  AngArrayT array;
  Iarray atoms(3);
  Iarray rnums(3);
  Narray names(3);

  for (AngleArray::const_iterator it = topIn.Angles().begin();
                                  it != topIn.Angles().end(); ++it)
  {
    SetAngParms( topIn, *it, atoms, rnums, names );
    array.push_back( AngT(atoms, rnums, names, topIn.AngleParm()[it->Idx()]) );
  }
  for (AngleArray::const_iterator it = topIn.AnglesH().begin();
                                  it != topIn.AnglesH().end(); ++it)
  {
    SetAngParms( topIn, *it, atoms, rnums, names );
    array.push_back( AngT(atoms, rnums, names, topIn.AngleParm()[it->Idx()]) );
  }
  return array;
}

/// Function for printing AngT
static void PrintAngT(CpptrajFile& output, Topology const& parm,
                      AngT const& ANG, char dir)
{
  PrintHeader(dir, ANG.Atoms(), parm, output);
  AngleParmType const& AP = ANG.Parms(); 
  output.Printf(" Tk=%g Teq=%g\n", AP.Tk(), AP.Teq());
}

/// Function for printing AngleParmType
static void PrintAngP(CpptrajFile& output, Topology const& parm,
                      AngleParmType const& first, char dir)
{
  output.Printf("%c Tk= %g  Teq= %g\n", dir, first.Tk(), first.Teq());
}

// -----------------------------------------------------------------------------
static inline void SetBndParms(Topology const& topIn, BondType const& T,
                               Iarray& atoms, Iarray& rnums, Narray& names)
{
  atoms[0] = T.A1();
  atoms[1] = T.A2();
  rnums[0] = topIn[T.A1()].ResNum();
  rnums[1] = topIn[T.A2()].ResNum();
  names[0] = topIn[T.A1()].Name();
  names[1] = topIn[T.A2()].Name();
}

/// Hold bond parameters
typedef ParmT<BondParmType> BndT;

/// Array of bond atoms and parameters
typedef std::vector<BndT> BndArrayT;

/// \return Array of angle atoms and parameters from given Topology
static BndArrayT BndArray( Topology const& topIn ) {
  BndArrayT array;
  Iarray atoms(2);
  Iarray rnums(2);
  Narray names(2);

  for (BondArray::const_iterator it = topIn.Bonds().begin();
                                 it != topIn.Bonds().end(); ++it)
  {
    SetBndParms( topIn, *it, atoms, rnums, names );
    array.push_back( BndT(atoms, rnums, names, topIn.BondParm()[it->Idx()]) );
  }
  for (BondArray::const_iterator it = topIn.BondsH().begin();
                                 it != topIn.BondsH().end(); ++it)
  {
    SetBndParms( topIn, *it, atoms, rnums, names );
    array.push_back( BndT(atoms, rnums, names, topIn.BondParm()[it->Idx()]) );
  }
  return array;
}

/// Function for printing BndT
static void PrintBndT(CpptrajFile& output, Topology const& parm,
                      BndT const& BND, char dir)
{
  PrintHeader(dir, BND.Atoms(), parm, output);
  BondParmType const& BP = BND.Parms(); 
  output.Printf(" Rk=%g Req=%g\n", BP.Rk(), BP.Req());
}

/// Function for printing BondParmType
static void PrintBndP(CpptrajFile& output, Topology const& parm,
                      BondParmType const& first, char dir)
{
  output.Printf("%c Rk= %g  Req= %g\n", dir, first.Rk(), first.Req());
}

// -----------------------------------------------------------------------------
/// Class template for comparing two arrays of a given parameter type
template <class T> class Diff {
  public:
    /// Parameter type array
    typedef std::vector< T > ArrayType;
    /// Function pointer for printing parameter type
    typedef void (*PrintFxnType)(CpptrajFile&, Topology const&, T const&, char);
    Diff() {}
    /** Compare two parameter arrays and report differences.
      * \param a1_in First parameter array ('<').
      * \param a2_in Second parameter array ('>').
      */
    void Compare(const ArrayType& a1_in, const ArrayType& a2_in, PrintFxnType fxnIn,
                 CpptrajFile& output, Topology const& parm1, Topology const& parm2)
    {
      ArrayType a1 = a1_in;
      ArrayType a2 = a2_in;
      std::sort( a1.begin(), a1.end() );
/*      CpptrajFile out1;
      out1.OpenWrite("out1");
      out1.Printf("'%s': %zu parameters.\n", parm1.c_str(), a1.size());
      for (typename ArrayType::const_iterator it = a1.begin(); it != a1.end(); ++it)
        fxnIn(out1, parm1, *it, ' ');
      out1.CloseFile();*/

      std::sort( a2.begin(), a2.end() );
/*      CpptrajFile out2;
      out2.OpenWrite("out2");
      out2.Printf("'%s': %zu parameters.\n", parm2.c_str(), a2.size());
      for (typename ArrayType::const_iterator it = a2.begin(); it != a2.end(); ++it)
        fxnIn(out2, parm2, *it, ' ');
      out2.CloseFile();*/

      typename ArrayType::const_iterator first1 = a1.begin();
      typename ArrayType::const_iterator first2 = a2.begin();
      while (first1 != a1.end() && first2 != a2.end()) {
        if (*first1 < *first2) {
          fxnIn(output, parm1, *first1, '<');
          ++first1;
        } else if (*first2 < *first1) {
          fxnIn(output, parm2, *first2, '>');
          ++first2;
        } else { ++first1; ++first2; }
      }
      while (first1 != a1.end())
        fxnIn(output, parm1, *(first1++), '<');
      while (first2 != a2.end())
        fxnIn(output, parm2, *(first2++), '>');
    }
};

// -----------------------------------------------------------------------------
/// Compare two topologies, find differences
Exec::RetType Exec_CompareTop::Execute(CpptrajState& State, ArgList& argIn)
{
  mprintf("Warning: THIS COMMAND IS NOT FULLY IMPLEMENTED.\n");
  Topology* parm1 = State.DSL().GetTopology( argIn );
  Topology* parm2 = State.DSL().GetTopology( argIn );
  if (parm1 == 0 || parm2 == 0) {
    mprinterr("Error: Specify two topologies.\n");
    return CpptrajState::ERR;
  }
  Topology const& p1 = static_cast<Topology const&>( *parm1 );
  Topology const& p2 = static_cast<Topology const&>( *parm2 );
  CpptrajFile output;
  output.OpenWrite( argIn.GetStringKey("out") );
  mprintf("\tOutput to '%s'\n", output.Filename().full());
  output.Printf("#< %s\n#> %s\n", p1.c_str(), p2.c_str());
  bool cmp_atype = argIn.hasKey("atype");
  bool cmp_lj = argIn.hasKey("lj");
  bool cmp_bnd = argIn.hasKey("bnd");
  bool cmp_ang = argIn.hasKey("ang");
  bool cmp_dih = argIn.hasKey("dih");
  if (!cmp_atype && !cmp_lj && !cmp_bnd && !cmp_ang && !cmp_dih) {
    cmp_atype = cmp_lj = cmp_bnd = cmp_ang = cmp_dih = true;
  }
  if (cmp_atype) {
    // Atom Types
    output.Printf("# Atom types\n");
    Diff<NameType> diff_atype;
    diff_atype.Compare( AtypeArray(p1), AtypeArray(p2), PrintAtype, output, p1, p2 );
  }
  if (cmp_lj) {
    // LJ params
    output.Printf("# LJ params\n");
    Diff<LJatom> diff_lj;
    diff_lj.Compare( LJarray(p1), LJarray(p2), PrintLJatom, output, p1, p2 );
  }
  if (cmp_bnd) {
    // Bonds
    output.Printf("# Bonds\n");
    Diff<BndT> diff_bnd;
    diff_bnd.Compare( BndArray(p1), BndArray(p2), PrintBndT, output, p1, p2 );
    // Bond parameters
    output.Printf("# Bond Parameters\n");
    Diff<BondParmType> diff_bndP;
    diff_bndP.Compare( p1.BondParm(), p2.BondParm(), PrintBndP, output, p1, p2 );
  }
  if (cmp_ang) {
    // Angles
    output.Printf("# Angles\n");
    Diff<AngT> diff_ang;
    diff_ang.Compare( AngArray(p1), AngArray(p2), PrintAngT, output, p1, p2 );
    // Angle parameters
    output.Printf("# Angle Parameters\n");
    Diff<AngleParmType> diff_angP;
    diff_angP.Compare( p1.AngleParm(), p2.AngleParm(), PrintAngP, output, p1, p2 );
  }
  if (cmp_dih) {
    // Dihedrals
    output.Printf("# Dihedrals\n");
    Diff<DihT> diff_dih;
    diff_dih.Compare( DihArray(p1), DihArray(p2), PrintDihT, output, p1, p2 );
    // Dihedral parameters
    output.Printf("# Dihedral Parameters\n");
    Diff<DihedralParmType> diff_dihP;
    diff_dihP.Compare( p1.DihedralParm(), p2.DihedralParm(), PrintDihP, output, p1, p2 );
  }
  output.CloseFile();
  return CpptrajState::OK;
}
