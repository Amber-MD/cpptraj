#include <cstring>
#include "SequenceAlign.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "Trajout_Single.h"
// EXPERIMENTAL ALPHA CODE
/** Idea is to take a reference PDB and a BLAST-type alignment
  * with format:
  * Query  1    MIVNVIQKDRL-KEQKLQFIRNHQQAFDVEPIRPLPLFEDFVTSIEGDCSLEASCKIESD  59
  *             MI+     D   KE+++Q +R+H ++FDVE   PLPLFE  V S++    LE S K++
  * Sbjct  1    MIMTTTWPDSYAKERRIQRLRHHFESFDVERAFPLPLFEQAVLSLDSCPLLEPSFKVQEG  60
  * and split the reference PDB into segments with the correct residue IDs so
  * it can then be put together.
  *
  * Usage: sequencealign <ref keyword> blastfile <file> out <file> [{pdb | mol2}]
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
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: Must specify output file.\n");
    return 1;
  }
  TrajectoryFile::TrajFormatType fmt = TrajectoryFile::GetFormatFromArg(argIn);
  if (fmt != TrajectoryFile::PDBFILE && fmt != TrajectoryFile::MOL2FILE)
    fmt = TrajectoryFile::PDBFILE; // Default to PDB
  Trajout_Single trajout;

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
  for (unsigned int sres = 0; sres != Sbjct.size(); sres++) {
    mprintf("%-u %3s %i", sres+1, Residue::ConvertResName(Sbjct[sres]), Smap[sres]+1);
    const char* qres = "";
    if (Smap[sres] != -1)
      qres = Residue::ConvertResName(Query[Smap[sres]]);
    mprintf(" %3s\n", qres);
  }
  // Check that query residues match reference.
  for (unsigned int sres = 0; sres != Sbjct.size(); sres++) {
    int qres = Smap[sres];
    if (qres != -1) {
      if (Query[qres] != qref.Parm().Res(qres).SingleCharName()) {
        mprintf("Warning: Potential residue mismatch: Query %s reference %s\n",
                Residue::ConvertResName(Query[qres]), qref.Parm().Res(qres).c_str());
      }
    }
  }
  // Build subject using coordinate from reference.
  //AtomMask sMask; // Contain atoms that should be in sTop
  Topology sTop;
  Frame sFrame;
  Iarray placeHolder; // Atom indices of placeholder residues.
  for (unsigned int sres = 0; sres != Sbjct.size(); sres++) {
    int qres = Smap[sres];
    NameType SresName( Residue::ConvertResName(Sbjct[sres]) );
    if (qres != -1) {
      Residue const& QR = qref.Parm().Res(qres);
      Residue SR(SresName, sres+1, ' ', QR.ChainID());
      if (Query[qres] == Sbjct[sres]) { // Exact match. All non-H atoms.
        for (int qat = QR.FirstAtom(); qat != QR.LastAtom(); qat++)
        {
          if (qref.Parm()[qat].Element() != Atom::HYDROGEN)
            sTop.AddTopAtom( qref.Parm()[qat], SR, qref.Coord().XYZ(qat) );
            sFrame.AddXYZ( qref.Coord().XYZ(qat) );
            //sMask.AddAtom(qat);
        }
      } else { // Partial match. Copy only backbone and CB.
        for (int qat = QR.FirstAtom(); qat != QR.LastAtom(); qat++)
        {
          if ( qref.Parm()[qat].Name().Match("N" ) ||
               qref.Parm()[qat].Name().Match("CA") ||
               qref.Parm()[qat].Name().Match("CB") ||
               qref.Parm()[qat].Name().Match("C" ) ||
               qref.Parm()[qat].Name().Match("O" ) )
          {
            sTop.AddTopAtom( qref.Parm()[qat], SR, qref.Coord().XYZ(qat) );
            sFrame.AddXYZ( qref.Coord().XYZ(qat) );
          }
        }
      }
    } else {
      // Residue in query does not exist for subject. Just put placeholder CA for now.
      Vec3 Zero(0.0);
      placeHolder.push_back( sTop.Natom() );
      sTop.AddTopAtom( Atom("CA", "C "), Residue(SresName, sres+1, ' ', ' '), Zero.Dptr() );
      sFrame.AddXYZ( Zero.Dptr() );
    }
  }
  sTop.SetRes(sTop.Nres()-1).SetLastAtom( sTop.Natom() ); // FIXME to not call CommonSetup
  //sTop.PrintAtomInfo("*");
  mprintf("\tPlaceholder residue indices:");
  for (Iarray::const_iterator p = placeHolder.begin(); p != placeHolder.end(); ++p)
    mprintf(" %i", *p + 1);
  mprintf("\n");
  // Try to give placeholders more reasonable coordinates.
  if (!placeHolder.empty()) {
    Iarray current_indices;
    unsigned int pidx = 0;
    while (pidx < placeHolder.size()) {
      if (current_indices.empty()) {
        current_indices.push_back( placeHolder[pidx++] );
        // Search for the end of this segment
        for (; pidx != placeHolder.size(); pidx++) {
          if (placeHolder[pidx] - current_indices.back() > 1) break;
          current_indices.push_back( placeHolder[pidx] );
        }
        // DEBUG
        mprintf("\tSegment:");
        for (Iarray::const_iterator it = current_indices.begin();
                                    it != current_indices.end(); ++it)
          mprintf(" %i", *it + 1);
        // Get coordinates of residues bordering segment.
        int prev_res = sTop[current_indices.front()].ResNum() - 1;
        int next_res = sTop[current_indices.back() ].ResNum() + 1;
        mprintf(" (prev_res=%i, next_res=%i)\n", prev_res+1, next_res+1);
        Vec3 prev_crd(sFrame.XYZ(current_indices.front() - 1));
        Vec3 next_crd(sFrame.XYZ(current_indices.back()  + 1));
        prev_crd.Print("prev_crd");
        next_crd.Print("next_crd");
        Vec3 crd_step = (next_crd - prev_crd) / (double)(current_indices.size()+1);
        crd_step.Print("crd_step");
        double* xyz = sFrame.xAddress() + (current_indices.front() * 3);
        for (unsigned int i = 0; i != current_indices.size(); i++, xyz += 3) {
          prev_crd += crd_step;
          xyz[0] = prev_crd[0];
          xyz[1] = prev_crd[1];
          xyz[2] = prev_crd[2];
        }
        current_indices.clear();
      }
    }
  }
  //Topology* sTop = qref.Parm().partialModifyStateByMask( sMask );
  //if (sTop == 0) return 1;
  //Frame sFrame(qref.Coord(), sMask);
  // Write output traj
  if (trajout.InitTrajWrite(outfilename, argIn, &sTop, fmt)) return 1;
  if (trajout.SetupTrajWrite(&sTop)) return 1;
  if (trajout.WriteSingle(0, sFrame)) return 1;
  trajout.EndTraj();
  return 0;
}
