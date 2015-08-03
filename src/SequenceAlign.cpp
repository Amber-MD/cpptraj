#include <cstring>
#include "SequenceAlign.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
// EXPERIMENTAL ALPHA CODE
/** Idea is to take a reference PDB and a BLAST-type alignment
  * with format:
  * Query  1    MIVNVIQKDRL-KEQKLQFIRNHQQAFDVEPIRPLPLFEDFVTSIEGDCSLEASCKIESD  59
  *             MI+     D   KE+++Q +R+H ++FDVE   PLPLFE  V S++    LE S K++
  * Sbjct  1    MIMTTTWPDSYAKERRIQRLRHHFESFDVERAFPLPLFEQAVLSLDSCPLLEPSFKVQEG  60
  * and split the reference PDB into segments with the correct residue IDs so
  * it can then be put together.
  *
  * Usage: sequencealign <ref keyword> blastfile <file>
  */
int SequenceAlign(CpptrajState& State, ArgList& argIn) {
  std::string blastfile = argIn.GetStringKey("blastfile");
  if (blastfile.empty()) {
    mprinterr("Error: 'blastfile' must be specified.\n");
    return 1;
  }
  ReferenceFrame qref = State.DSL()->GetReferenceFrame(argIn);
  if (qref.error() || qref.empty()) {
    mprinterr("Error: Must specify reference structure for query.\n");
    return 1;
  }

  // Load blast file
  mprintf("\tReading BLAST alignment from '%s'\n", blastfile.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead( blastfile )) return 1;
  // Seek down to first Query line.
  const char* ptr = infile.Line();
  bool atFirstQuery = false;
  while (ptr != 0) {
    if (*ptr == 'Q') {
      if ( strncmp(ptr, "Query", 5) == 0 ) {
        atFirstQuery = true;
        break;
      }
    }
    ptr = infile.Line();
  }
  if (!atFirstQuery) {
    mprinterr("Error: 'Query' not found.\n");
    return 1;
  }

  // Read alignment. Replacing query with subject.
  typedef std::vector<char> Carray;
  typedef std::vector<int> Iarray;
  Carray Query; // Query residues
  Carray Sbjct; // Sbjct residues
  Iarray Smap;  // Smap[Sbjct index] = Query index
  while (ptr != 0) {
    const char* qline = ptr;           // query line
    const char* aline = infile.Line(); // alignment line
    const char* sline = infile.Line(); // subject line
    if (aline == 0 || sline == 0) {
      mprinterr("Error: Missing alignment line or subject line after Query:\n");
      mprinterr("Error:  %s", qline);
      return 1;
    }
    for (int idx = 12; qline[idx] != ' '; idx++) {
      if (qline[idx] == '-') {
        // Sbjct does not have corresponding res in Query
        Smap.push_back(-1);
        Sbjct.push_back( sline[idx] );
      } else if (sline[idx] == '-') {
        // Query does not have a corresponding res in Sbjct
        Query.push_back( qline[idx] );
      } else {
        // Direct Query to Sbjct map
        Smap.push_back( Query.size() );
        Sbjct.push_back( sline[idx] );
        Query.push_back( qline[idx] );
      }
    }
    // Scan to next Query 
    ptr = infile.Line();
    while (ptr != 0) {
      if (*ptr == 'Q') {
        if ( strncmp(ptr, "Query", 5) == 0 ) break;
      }
      ptr = infile.Line();
    }
  }
  // DEBUG
  mprintf("  Map of Sbjct to Query:\n");
  for (unsigned int i = 0; i != Sbjct.size(); i++) {
    mprintf("%-u %3s %i", i+1, Residue::ConvertResName(Sbjct[i]), Smap[i]+1);
    const char* qres = "";
    if (Smap[i] != -1)
      qres = Residue::ConvertResName(Query[Smap[i]]);
    mprintf(" %3s\n", qres);
  }
  // Check that query residues match reference.
  for (unsigned int i = 0; i != Sbjct.size(); i++) {
    int qres = Smap[i];
    if (qres != -1) {
      if (Query[qres] != qref.Parm().Res(qres).SingleCharName()) {
        mprintf("Warning: Potential residue mismatch: Query %s reference %s\n",
                Residue::ConvertResName(Query[qres]), qref.Parm().Res(qres).c_str());
      }
    }
  }
 
  return 0;
}
