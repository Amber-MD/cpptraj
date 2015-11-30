#include <algorithm> // sort
#include "Cmd_CompareTop.h"
#include "CpptrajStdio.h"

void Help_CompareTop() {
  mprintf("\t{%s} {%s} [out <file>]\n", DataSetList::TopArgs, DataSetList::TopArgs);
}

static inline void PrintDih(CpptrajFile& output, Topology* parm,
                            DihedralArray::const_iterator const& first, char dir)
{
  output.Printf("%c %s - %s - %s - %s\n", dir,
                parm->AtomMaskName(first->A1()).c_str(),
                parm->AtomMaskName(first->A2()).c_str(),
                parm->AtomMaskName(first->A3()).c_str(),
                parm->AtomMaskName(first->A4()).c_str());
}

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
  // Dihedrals
  DihedralArray d1 = parm1->Dihedrals();
  DihedralArray d2 = parm2->Dihedrals();
  std::sort( d1.begin(), d1.end() );
  std::sort( d2.begin(), d2.end() );
  DihedralArray::const_iterator first1 = d1.begin();
  DihedralArray::const_iterator first2 = d2.begin();
  while (first1 != d1.end() && first2 != d2.end()) {
    if (*first1 < *first2) {
      PrintDih(output, parm1, first1, '<');
      ++first1;
    } else if (*first2 < *first1) {
      PrintDih(output, parm2, first2, '>');
      ++first2;
    } else { ++first1; ++first2; }
  }
  while (first1 != d1.end())
    PrintDih(output, parm1, first1++, '<');
  while (first2 != d2.end())
    PrintDih(output, parm2, first2++, '>');
  output.CloseFile();
  return 0;
}
