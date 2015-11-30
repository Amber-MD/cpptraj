#include <algorithm> // sort
#include "Cmd_CompareTop.h"
#include "CpptrajStdio.h"

void Help_CompareTop() {
  mprintf("\t{%s} {%s} [out <file>]\n", DataSetList::TopArgs, DataSetList::TopArgs);
}

/// Function for printing DihedralType
static inline void PrintDihT(CpptrajFile& output, Topology* parm,
                             DihedralType const& first, char dir)
{
  output.Printf("%c %s - %s - %s - %s\n", dir,
                parm->AtomMaskName(first.A1()).c_str(),
                parm->AtomMaskName(first.A2()).c_str(),
                parm->AtomMaskName(first.A3()).c_str(),
                parm->AtomMaskName(first.A4()).c_str());
}

/// Function for printing AngleType
static inline void PrintAngT(CpptrajFile& output, Topology* parm,
                             AngleType const& first, char dir)
{
  output.Printf("%c %s - %s - %s\n", dir,
                parm->AtomMaskName(first.A1()).c_str(),
                parm->AtomMaskName(first.A2()).c_str(),
                parm->AtomMaskName(first.A3()).c_str());
}

/// Class template for comparing two arrays of a given parameter type
template <class T> class Diff {
  public:
    /// Parameter type array
    typedef std::vector< T > ArrayType;
    /// Function pointer for printing parameter type
    typedef void (*PrintFxnType)(CpptrajFile&, Topology*, T const&, char);
    Diff() {}
    /** Compare two parameter arrays and report differences.
      * \param a1_in First parameter array ('<').
      * \param a2_in Second parameter array ('>').
      */
    void Compare(const ArrayType& a1_in, const ArrayType& a2_in, PrintFxnType fxnIn,
                 CpptrajFile& output, Topology* parm1, Topology* parm2)
    {
      ArrayType a1 = a1_in;
      ArrayType a2 = a2_in;
      std::sort( a1.begin(), a1.end() );
      std::sort( a2.begin(), a2.end() );
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
        fxnIn(output, parm2, *(first2++), '<');
    }
};

int CompareTop(CpptrajState& State, ArgList& argIn)
{
  Topology* parm1 = State.DSL()->GetTopology( argIn );
  Topology* parm2 = State.DSL()->GetTopology( argIn );
  if (parm1 == 0 || parm2 == 0) {
    mprinterr("Error: Specify two topologies.\n");
    return 1;
  }
  CpptrajFile output;
  output.OpenWrite( argIn.GetStringKey("out") );
  mprintf("\tOutput to '%s'\n", output.Filename().full());
  output.Printf("#< %s\n#> %s\n", parm1->c_str(), parm2->c_str());
  // Angles
  output.Printf("# Angles\n");
  Diff<AngleType> diff_ang;
  diff_ang.Compare( parm1->Angles(), parm2->Angles(), PrintAngT, output, parm1, parm2 );
  // Dihedrals
  output.Printf("# Dihedrals\n");
  Diff<DihedralType> diff_dih;
  diff_dih.Compare( parm1->Dihedrals(), parm2->Dihedrals(), PrintDihT, output, parm1, parm2 );
  diff_dih.Compare( parm1->DihedralsH(), parm2->DihedralsH(), PrintDihT, output, parm1, parm2 );

  output.CloseFile();
  return 0;
}
