#include <cstdlib> // atoi
#include <cctype> // isdigit
#include "Exec_ViewRst.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

void Exec_ViewRst::Help() const {
  mprintf("\tcrdset <coords set> rstfile <filename> mol2out <mol2 name> [noe]\n");
}

Exec::RetType Exec_ViewRst::Execute(CpptrajState& State, ArgList& argIn)
{
  mprintf("Warning: This command is experimental.\n");
  int Ntops = 1;
  enum TopType {TOP_STRONG=0, TOP_MEDIUM, TOP_WEAK, TOP_VERYWEAK};
  enum OutputType { ALL=0, BY_STRENGTH };
  OutputType outputMode = ALL;
  // Get output type
  if (argIn.hasKey("noe")) {
    mprintf("\tOutput to mol2 by restraint distances (NOE strength).\n");
    outputMode = BY_STRENGTH;
    Ntops = 4; // strong, med, weak, v. weak
  }
  // Get output mol2 name
  std::string mol2out = argIn.GetStringKey("mol2out");
  if (mol2out.empty()) {
    mprinterr("Error: Output mol2 name ('mol2out') not specfied.\n");
    return CpptrajState::ERR;
  }
  std::vector< std::string > OutNames;
  switch (outputMode) {
    case BY_STRENGTH:
      OutNames.push_back("strong."   + mol2out);
      OutNames.push_back("medium."   + mol2out);
      OutNames.push_back("weak."     + mol2out);
      OutNames.push_back("veryweak." + mol2out);
      break;
    case ALL:
      OutNames.push_back( mol2out );
      break;
  }
  // Attempt to get coords dataset from datasetlist
  std::string name = argIn.GetStringKey("crdset");
  DataSetList sets = State.DSL().SelectGroupSets( name, DataSet::COORDINATES );
  if (sets.empty()) {
    mprinterr("Error: Could not locate COORDS set corresponding to '%s'\n", name.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords* coords = (DataSet_Coords*)sets[0];
  mprintf("\tCoords set: '%s'\n", coords->legend());
  // Get file with amber restraints
  name = argIn.GetStringKey("rstfile");
  if (name.empty()) {
    mprinterr("Error: No Amber restraint file given.\n");
    return CpptrajState::ERR;
  }
  BufferedLine infile;
  if (infile.OpenFileRead( name )) return CpptrajState::ERR;
  // Add all atoms to each pseudo topology.
  std::vector<Topology> Pseudo( Ntops );
  for (int nt = 0; nt != Ntops; nt++) {
    for (Topology::atom_iterator atm = coords->Top().begin(); atm != coords->Top().end(); ++atm)
      Pseudo[nt].AddTopAtom( *atm, coords->Top().Res( atm->ResNum() ) );
  }
  // Process Amber restraint file. Create pseudo topology bonds.
  const char* ptr = infile.Line();
  const char* SEP = " ,=";
  int Ncols = infile.TokenizeLine(SEP);
  int currentCol = 0;
  bool inRst = false;
  std::vector<int> bondIndices;
  double Rvals[4]; // r1, r2, r3, r4
  std::fill(Rvals, Rvals+4, -1.0);
  enum StateType {PROCESS_TOKEN=0, NEXT_TOKEN, READ_ATOMS, READ_R1, READ_R2, READ_R3, READ_R4};
  //const char* STATE[] = {"PROCESS_TOKEN", "NEXT_TOKEN", "READ_ATOMS", "READ_R1", "READ_R2", "READ_R3", "READ_R4"};
  StateType currentState = PROCESS_TOKEN;
  while (ptr != 0) {
    // Load next line if necessary
    if (currentCol == Ncols) {
      ptr = infile.Line();
      if (ptr == 0)
        break;
      else {
        Ncols = infile.TokenizeLine(SEP);
        currentCol = 0;
      }
    }
    const char* currentToken = infile.NextToken();
    currentCol++;
    // If state is NEXT_TOKEN, now change to process token.
    if (currentState == NEXT_TOKEN)
      currentState = PROCESS_TOKEN;
    //mprintf("DEBUG: State='%s'  Token='%s'  Col=%i\n",STATE[currentState],currentToken,currentCol);
    // Skip comment lines
    if (currentCol == 1 && currentToken[0] == '#') {
      currentCol = Ncols;
      continue;
    }

    // Check state
    if (currentState == READ_ATOMS) {
      // Reading list of atoms.
      if (currentToken[0] != '0' && isdigit(currentToken[0]))
        bondIndices.push_back( atoi(currentToken) - 1 );
      else {
        // If atoms ended on '0' we need next token. Otherwise at next token.
        if (currentToken[0] == '0')
          currentState = NEXT_TOKEN;
        else
          currentState = PROCESS_TOKEN;
      }
    } else if (currentState == READ_R1) {
      Rvals[0] = atof(currentToken);
      currentState = NEXT_TOKEN;
    } else if (currentState == READ_R2) {
      Rvals[1] = atof(currentToken);
      currentState = NEXT_TOKEN;
    } else if (currentState == READ_R3) {
      Rvals[2] = atof(currentToken);
      currentState = NEXT_TOKEN;
    } else if (currentState == READ_R4) {
      Rvals[3] = atof(currentToken);
      currentState = NEXT_TOKEN;
    }

    if (currentState == PROCESS_TOKEN) {
      // Process token
      std::string tkn(currentToken);
      if (tkn == "&rst") {
        if (inRst) {
          mprinterr("Error: Line %i: New restraint started before end of old restraint.\n",
                    infile.LineNumber());
          return CpptrajState::ERR;
        }
        inRst = true;
      } else if (tkn =="&end") {
        if (!inRst) {
          mprinterr("Error: Line %i: '&end' encountered before '&rst'\n", infile.LineNumber());
          return CpptrajState::ERR;
        }
        if (bondIndices.size() != 2) {
          mprinterr("Error: Line %i: Restraint had %zu atom indices, expected only 2.\n",
                    infile.LineNumber(), bondIndices.size());
          return CpptrajState::ERR;
        }
        int ntop = 0;
        if (outputMode == BY_STRENGTH) {
          // Rvals[2] is r3
          if (Rvals[2] < 3.0)
            ntop = (int)TOP_STRONG;
          else if (Rvals[2] < 4.1)
            ntop = (int)TOP_MEDIUM;
          else if (Rvals[2] < 5.5)
            ntop = (int)TOP_WEAK;
          else
            ntop = (int)TOP_VERYWEAK;
        }
        Pseudo[ntop].AddBond( bondIndices[0], bondIndices[1] );
        if (State.Debug() > 0)
          mprintf("Bonding indices: %i %i  Rvals: %g %g %g %g {%i}\n",
                  bondIndices[0], bondIndices[1], Rvals[0], Rvals[1], Rvals[2], Rvals[3], ntop);
        bondIndices.clear();
        inRst = false;
      } else if (tkn == "iat") { // TODO check inRst?
        currentState = READ_ATOMS;
      } else if (tkn == "r1") {
        currentState = READ_R1;
      } else if (tkn == "r2") {
        currentState = READ_R2;
      } else if (tkn == "r3") {
        currentState = READ_R3;
      } else if (tkn == "r4") {
        currentState = READ_R4;
      }
    }
  } // END loop over file
  // Write output mol2. Use first frame of input coords.
  Frame frameOut = coords->AllocateFrame();
  coords->GetFrame(0, frameOut);
  for (int nt = 0; nt != Ntops; nt++) {
    Trajout_Single trajout;
    if (trajout.PrepareTrajWrite(OutNames[nt], ArgList(), &(Pseudo[nt]), CoordinateInfo(), 1,
                                 TrajectoryFile::MOL2FILE))
      return CpptrajState::ERR;
    if (trajout.WriteSingle(0, frameOut)) return CpptrajState::ERR;
    trajout.EndTraj();
  }

  return CpptrajState::OK;
}
