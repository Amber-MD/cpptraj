#include <cstring>
#include "Exec_SequenceAlign.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "Trajout_Single.h"
#include "StringRoutines.h"
// EXPERIMENTAL ALPHA CODE
void Exec_SequenceAlign::Help() const {
  mprintf("\t%s blastfile <file> out <file>\n"
          "\t[{pdb | mol2}] [<trajout args>] [smaskoffset <#>] [qmaskoffset <#>]\n"
          "  blastfile: File containing sequence alignment.\n"
          "  out: File to write resulting structure to.\n"
          "  [{pdb | mol2}] Format of structure (default pdb).\n"
          "  [smaskoffset]|[qmaskoffset]: Resdiue offset for output.\n"
          "Given a reference structure and BLAST-like sequence alignment of format:\n"
          "    Query  1    MIV...\n"
          "                MI+...\n"
          "    Sbjct  1    MIM...\n"
          "where the Query sequence corresponds to the reference, output a structure\n"
          "with the Sbjct sequence. All hydrogens will be removed. Non-matching residues\n"
          "will include only backbone and CB (if applicable) atoms.\n",
          DataSetList::RefArgs);
}
    
Exec::RetType Exec_SequenceAlign::Execute(CpptrajState& State, ArgList& argIn) {
  mprintf("Warning: THIS COMMAND IS NOT FULLY IMPLEMENTED.\n");
  std::string blastfile = argIn.GetStringKey("blastfile");
  if (blastfile.empty()) {
    mprinterr("Error: 'blastfile' must be specified.\n");
    return CpptrajState::ERR;
  }
  ReferenceFrame qref = State.DSL().GetReferenceFrame(argIn);
  if (qref.error() || qref.empty()) {
    mprinterr("Error: Must specify reference structure for query.\n");
    return CpptrajState::ERR;
  }
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: Must specify output file.\n");
    return CpptrajState::ERR;
  }
  TrajectoryFile::TrajFormatType fmt = TrajectoryFile::GetFormatFromArg(argIn);
  if (fmt != TrajectoryFile::PDBFILE && fmt != TrajectoryFile::MOL2FILE)
    fmt = TrajectoryFile::PDBFILE; // Default to PDB
  int smaskoffset = argIn.getKeyInt("smaskoffset", 0) + 1;
  int qmaskoffset = argIn.getKeyInt("qmaskoffset", 0) + 1;

  // Load blast file
  mprintf("\tReading BLAST alignment from '%s'\n", blastfile.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead( blastfile )) return CpptrajState::ERR;
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
    return CpptrajState::ERR;
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
      return CpptrajState::ERR;
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
  std::string SmaskExp, QmaskExp;
  if (State.Debug() > 0) mprintf("  Map of Sbjct to Query:\n");
  for (int sres = 0; sres != (int)Sbjct.size(); sres++) {
    if (State.Debug() > 0)
      mprintf("%-i %3s %i", sres+smaskoffset, Residue::ConvertResName(Sbjct[sres]),
              Smap[sres]+qmaskoffset);
    const char* qres = "";
    if (Smap[sres] != -1) {
      qres = Residue::ConvertResName(Query[Smap[sres]]);
      if (SmaskExp.empty())
        SmaskExp.assign( integerToString(sres+smaskoffset) );
      else
        SmaskExp.append( "," + integerToString(sres+smaskoffset) );
      if (QmaskExp.empty())
        QmaskExp.assign( integerToString(Smap[sres]+qmaskoffset) );
      else
        QmaskExp.append( "," + integerToString(Smap[sres]+qmaskoffset) );

    }
    if (State.Debug() > 0) mprintf(" %3s\n", qres);
  }
  mprintf("Smask: %s\n", SmaskExp.c_str());
  mprintf("Qmask: %s\n", QmaskExp.c_str());
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
            sTop.AddTopAtom( qref.Parm()[qat], SR );
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
            sTop.AddTopAtom( qref.Parm()[qat], SR );
            sFrame.AddXYZ( qref.Coord().XYZ(qat) );
          }
        }
      }
    } else {
      // Residue in query does not exist for subject. Just put placeholder CA for now.
      Vec3 Zero(0.0);
      placeHolder.push_back( sTop.Natom() );
      sTop.AddTopAtom( Atom("CA", "C "), Residue(SresName, sres+1, ' ', ' ') );
      sFrame.AddXYZ( Zero.Dptr() );
    }
  }
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
  //if (sTop == 0) return CpptrajState::ERR;
  //Frame sFrame(qref.Coord(), sMask);
  // Write output traj
  Trajout_Single trajout;
  if (trajout.PrepareTrajWrite(outfilename, argIn, &sTop, CoordinateInfo(), 1, fmt))
    return CpptrajState::ERR;
  if (trajout.WriteSingle(0, sFrame)) return CpptrajState::ERR;
  trajout.EndTraj();
  return CpptrajState::OK;
}
