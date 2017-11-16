#include <cstdio>
#include <cstdlib>
#include <cctype>
#include "DataIO_Cpout.h"
#include "CpptrajStdio.h"
#include "DataSet_pH.h"
#include "DataSet_pH_REMD.h"

/// CONSTRUCTOR
DataIO_Cpout::DataIO_Cpout() :
  type_(NONE),
  original_pH_(0.0)
{
  SetValid( DataSet::PH );
  SetValid( DataSet::PH_REMD );
}

const char* DataIO_Cpout::FMT_REDOX_ = "Redox potential: %f V";

const char* DataIO_Cpout::FMT_PH_ = "Solvent pH: %f";

// DataIO_Cpout::ID_DataFormat()
bool DataIO_Cpout::ID_DataFormat(CpptrajFile& infile)
{
  bool iscpout = false;
  type_ = NONE;
  if (!infile.OpenFile()) {
    const char* ptr = infile.NextLine();
    if (ptr != 0) {
      if (sscanf(ptr, FMT_REDOX_, &original_pH_) == 1) {
        type_ = REDOX;
      } else if (sscanf(ptr, FMT_PH_, &original_pH_) == 1) {
        type_ = PH;
      }
      if (type_ != NONE) {
        ptr = infile.NextLine();
        int step_size;
        if (ptr != 0) {
          iscpout = (sscanf(ptr, "Monte Carlo step size: %d", &step_size) == 1);
        }
      }
    }
    infile.CloseFile();
  }
  return iscpout;
}

// DataIO_Cpout::ReadHelp()
void DataIO_Cpout::ReadHelp()
{
  mprintf("\tcpin <file> : CPIN file name.\n");
}

// DataIO_Cpout::processReadArgs()
int DataIO_Cpout::processReadArgs(ArgList& argIn)
{
  cpin_file_ = argIn.GetStringKey("cpin");
  return 0;
}

/** Read CPIN file to get state information for each residue. */
int DataIO_Cpout::ReadCpin(FileName const& fname) {
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) return 1;

  enum VarType { NONE = 0, CHRGDAT, PROTCNT, RESNAME, RESSTATE, NUM_STATES, 
                 NUM_ATOMS, FIRST_STATE, FIRST_CHARGE, TRESCNT };
  VarType vtype = NONE;

  enum NamelistMode { OUTSIDE_NAMELIST = 0, INSIDE_NAMELIST };
  NamelistMode nmode = OUTSIDE_NAMELIST;

  /// One for each residue.
  StateArray States;

  typedef std::vector<std::string> Sarray;
  Sarray resnames;
  std::string system;
  Iarray protcnt;
  int trescnt = 0;
  Darray charges; // TODO may not need these

  int stateinf_ridx = 0;

  std::string token;
  const char* ptr = infile.Line(); // &CNSTPH
  bool inQuote = false;
  while (ptr != 0) {
    for (const char* p = ptr; *p != '\0'; ++p)
    {
      if (nmode == OUTSIDE_NAMELIST) {
        if (isspace(*p)) {
          //mprintf("NAMELIST: %s\n", token.c_str());
          token.clear();
          nmode = INSIDE_NAMELIST;
        } else
          token += *p;
      } else {
        if ( *p == '/' || *p == '&' )
          break;
        else if ( *p == '\'' || *p == '"' ) {
          if (inQuote)
            inQuote = false;
          else
            inQuote = true;
        // ---------------------------------------
        } else if ( *p == '=' ) {
          // Variable. Figure out which one.
          //mprintf("VAR: %s\n", token.c_str());
          if ( token == "CHRGDAT" )
            vtype = CHRGDAT;
          else if (token == "PROTCNT")
            vtype = PROTCNT;
          else if (token == "RESNAME")
            vtype = RESNAME;
          else if (token == "TRESCNT")
            vtype = TRESCNT;
          else if (token.compare(0,8,"STATEINF") == 0) {
            // Expect STATEINF(<res>)%<field>
            ArgList stateinf(token, "()%=");
            if (stateinf.Nargs() != 3) {
              mprintf("Error: Malformed STATEINF: %s\n", token.c_str());
              return 1;
            }
            stateinf_ridx = atoi(stateinf[1].c_str());
            //mprintf("DEBUG Res %i : %s\n", stateinf_ridx, stateinf[2].c_str());
            if (stateinf_ridx >= (int)States.size())
              States.resize(stateinf_ridx+1); // TODO bounds check
            if (stateinf[2] == "NUM_STATES")
              vtype = NUM_STATES;
            else if (stateinf[2] == "NUM_ATOMS")
              vtype = NUM_ATOMS;
            else if (stateinf[2] == "FIRST_STATE")
              vtype = FIRST_STATE;
            else if (stateinf[2] == "FIRST_CHARGE")
              vtype = FIRST_CHARGE;
            else
              vtype = NONE;
          } else
            vtype = NONE;
          token.clear();
        // ---------------------------------------
        } else if (*p == ',') {
          // Value. Assign to appropriate variable
          //mprintf("\tValue: %s mode %i\n", token.c_str(), (int)vtype);
          if (vtype == CHRGDAT)
            charges.push_back( atof( token.c_str() ) );
          else if (vtype == PROTCNT)
            protcnt.push_back( atoi( token.c_str() ) );
          else if (vtype == RESNAME) {
            if (token.compare(0,6,"System")==0)
              system = token.substr(7);
            else
              resnames.push_back( token.substr(8) );
          } else if (vtype == TRESCNT)
            trescnt = atoi( token.c_str() );
          else if (vtype == NUM_STATES)
            States[stateinf_ridx].num_states_ = atoi( token.c_str() );
          else if (vtype == NUM_ATOMS)
            States[stateinf_ridx].num_atoms_  = atoi( token.c_str() );
          else if (vtype == FIRST_STATE)
            States[stateinf_ridx].first_state_  = atoi( token.c_str() );
          else if (vtype == FIRST_CHARGE)
            States[stateinf_ridx].first_charge_ = atoi( token.c_str() );
          token.clear();
        } else if (!isspace(*p) || inQuote)
          token += *p;
      }
    }
    ptr = infile.Line();
  }

  // DEBUG
  if (debug_ > 1) {
    mprintf("%zu charges.\n", charges.size());
    int col = 0;
    for (Darray::const_iterator c = charges.begin(); c != charges.end(); ++c) {
      mprintf(" %12.4f", *c);
      if (++col == 5) {
        mprintf("\n");
        col = 0;
      }
    }
    if (col != 0) mprintf("\n");
  }
  mprintf("\tSystem: %s\n", system.c_str());
  if (debug_ > 0) {
    mprintf("%zu protcnt=", protcnt.size());
    for (Iarray::const_iterator p = protcnt.begin(); p != protcnt.end(); ++p)
      mprintf(" %i", *p);
    mprintf("\n");
    mprintf("trescnt = %i\n", trescnt);
    for (StateArray::const_iterator it = States.begin(); it != States.end(); ++it)
      mprintf("\tnum_states= %i  num_atoms= %i  first_charge= %i  first_state= %i\n",
              it->num_states_, it->num_atoms_, it->first_charge_, it->first_state_);
    for (Sarray::const_iterator it = resnames.begin(); it != resnames.end(); ++it)
      mprintf("\t%s\n", it->c_str());
  }
  // Checks
  if (trescnt != (int)States.size()) {
    mprinterr("Error: Number of states in CPIN (%zu) != TRESCNT in CPIN (%i)\n",
              States.size(), trescnt);
    return 1;
  }
  if (trescnt != (int)resnames.size()) {
    mprinterr("Error: Number of residues in CPIN (%zu) != TRESCNT in CPIN (%i)\n",
              resnames.size(), trescnt);
    return 1;
  }

  // Define residues
  Sarray::const_iterator rname = resnames.begin();
  for (StateArray::const_iterator it = States.begin(); it != States.end(); ++it, ++rname)
  {
    Iarray res_protcnt;
    int max_prots = -1;
    for (int j = 0; j < it->num_states_; j++) {
      res_protcnt.push_back( protcnt[it->first_state_ + j] );
      max_prots = std::max( max_prots, res_protcnt.back() );
    }
    ArgList split(*rname);
    if (split.Nargs() != 2) {
      mprinterr("Error: Malformed residue name/number '%s'\n", rname->c_str());
      return 1;
    }
    Residues_.push_back( 
      CphResidue(split[0], atoi(split[1].c_str()), res_protcnt, max_prots) );
  }

  return 0;
}

// DataIO_Cpout::ReadData()
int DataIO_Cpout::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  // Require a CPIN file. 
  if (cpin_file_.empty()) {
    mprinterr("Error: No CPIN file specified.\n");
    return 1;
  }
  Residues_.clear();
  if (ReadCpin( cpin_file_ )) return 1;

  // Open CPOUT file.
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) return 1;
  const char* ptr = infile.Line();

  // Determine type and number of residues. 
  const char* fmt = 0;
  const char* rFmt = 0;
  if (sscanf(ptr, FMT_REDOX_, &original_pH_) == 1) {
    type_ = REDOX;
    mprintf("\tRedOx output file.\n");
    fmt = FMT_REDOX_;
    rFmt = "Residue %d State: %d E: %f V";
  } else if (sscanf(ptr, FMT_PH_, &original_pH_) == 1) {
    type_ = PH;
    mprintf("\tConstant pH output file.\n");
    fmt = FMT_PH_;
    rFmt = "Residue %d State: %d pH: %f";
  } else {
    mprinterr("Error: Could not determine CPOUT file type.\n");
    return 1;
  }
  ptr = infile.Line(); // Monte Carlo step size
  ptr = infile.Line(); // Current MD time step
  ptr = infile.Line(); // Current MD time
  ptr = infile.Line(); // First residue
  int res, state;
  float pHval; 
  int nscan = sscanf(ptr, rFmt, &res, &state, &pHval);
  if (nscan == 2) {
    mprintf("\tNot from REMD.\n");
  } else if (nscan == 3) {
    mprintf("\tpH values from REMD detected.\n");
  } else {
    mprintf("Got %i values from first Residue line, expected only 2 or 3.\n", nscan);
    return 1;
  }
  infile.CloseFile();
  if (infile.OpenFileRead( fname )) return 1;

  // Allocate DataSets
  int err = 1;
  if (nscan == 2) {
    // Sorted constant pH
    err = ReadSorted(infile, dsl, dsname, fmt, rFmt);
  } else if (nscan == 3) {
    // Unsorted constant pH
    err = ReadUnsorted(infile, dsl, dsname, fmt, rFmt);
  }
  infile.CloseFile();
/*
  if (debug_ > 1) {
    for (DataSet_pH_REMD::const_iterator res = phdata->begin(); res != phdata->end(); ++res) {
      mprintf("DEBUG: Res %u:\n", res-phdata->begin());
      for (CphResidue::const_iterator state = res->begin();
                                               state != res->end(); ++state)
        mprintf(" %i", *state);
      mprintf("\n");
    }
    mprintf("DEBUG: pH values:\n");
    for (DataSet_pH_REMD::ph_iterator ph = phdata->pH_Values().begin();
                                 ph != phdata->pH_Values().end(); ++ph)
      mprintf(" %6.2f", *ph);
    mprintf("\n");
  }

  mprintf("\tTitratable Residues:\n");
  for (DataSet_pH_REMD::const_iterator res = phdata->begin(); res != phdata->end(); ++res)
    res->Print();
  mprintf("\t%u frames\n", nframes);
*/
  return err;
}

int DataIO_Cpout::ReadSorted(BufferedLine& infile, DataSetList& DSL, std::string const& dsname, const char* fmt, const char* rFmt) {
  // Sorted constant pH
  typedef std::vector<DataSet_pH*> Parray;
  Parray ResSets;
  ResSets.reserve( Residues_.size() );
  for (Rarray::iterator res = Residues_.begin(); res != Residues_.end(); ++res)
  {
    MetaData md( dsname, *(res->Name()), res->Num() );
    DataSet* ds = DSL.CheckForSet(md);
    if (ds == 0) {
      // New set
      ds = DSL.AddSet( DataSet::PH, md );
      if (ds == 0) return 1;
      ((DataSet_pH*)ds)->SetResidueInfo( *res );
      ((DataSet_pH*)ds)->Set_Solvent_pH( original_pH_ ); // TODO combine with above?
    } else {
      if (ds->Type() != DataSet::PH) {
        mprinterr("Error: Set '%s' type does not match, cannot append.\n", ds->legend());
        return 1;
      }
      mprintf("\tAppending to set '%s'\n", ds->legend());
      // TODO check # residues etc?
    }
    ResSets.push_back( (DataSet_pH*)ds );
  }
  int maxRes = (int)ResSets.size();

  unsigned int nframes = 0;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    float solvent_pH;
    if (sscanf(ptr, fmt, &solvent_pH) == 1) {
      // Full record
      //mprintf("DEBUG: pH= %f\n", solvent_pH);
      ptr = infile.Line(); // Monte Carlo step size
      ptr = infile.Line(); // Current MD time step
      int step;
      if (sscanf(ptr,"Time step: %d", &step) != 1) {
        mprinterr("Error: Could not get step.\n");
        return 1;
      }
      //mprintf("DEBUG: step= %i\n", step);
      ptr = infile.Line(); // Current time (ps)
      float time;
      if (sscanf(ptr,"Time: %f", &time) != 1) {
        mprinterr("Error: Could not get time.\n");
        return 1;
      }
      //mprintf("DEBUG: time= %f\n", time);
      ptr = infile.Line(); // Residue
    }
    // delta record or full record Residue read
    int res, state;
    while (sscanf(ptr, rFmt, &res, &state) == 2) {
      //mprintf("DEBUG: res= %i state= %i\n", res, state);
      if (res < maxRes)
        ResSets[res]->AddState( state );
      else {
        mprinterr("Error: Res %i in CPOUT > max # res in CPIN (%i)\n", res, maxRes);
        return 1;
      }
      ptr = infile.Line();
    }
    nframes++;
    ptr = infile.Line();
  }
  return 0;
}

int DataIO_Cpout::ReadUnsorted(BufferedLine& infile, DataSetList& DSL, std::string const& dsname, const char* fmt, const char* rFmt) {
  // Unsorted constant pH
  DataSet* ds = DSL.CheckForSet(dsname);
  if (ds == 0) {
    // New set
    ds = DSL.AddSet( DataSet::PH_REMD, dsname, "ph" );
    if (ds == 0) return 1;
    ((DataSet_pH_REMD*)ds)->SetResidueInfo( Residues_ );
  } else {
    if (ds->Type() != DataSet::PH_REMD) {
      mprinterr("Error: Set '%s' is not pH data.\n", ds->legend());
      return 1;
    }
    mprintf("\tAppending to set '%s'\n", ds->legend());
    // TODO check # residues etc?
  }
  DataSet_pH_REMD* phdata = (DataSet_pH_REMD*)ds;
  int maxRes = (int)phdata->Residues().size();
  Iarray resStates( phdata->Residues().size() );

  float solvent_pH = original_pH_;
  unsigned int nframes = 0;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (sscanf(ptr, fmt, &solvent_pH) == 1) {
      // Full record
      //mprintf("DEBUG: pH= %f\n", solvent_pH);
      ptr = infile.Line(); // Monte Carlo step size
      ptr = infile.Line(); // Current MD time step
      int step;
      if (sscanf(ptr,"Time step: %d", &step) != 1) {
        mprinterr("Error: Could not get step.\n");
        return 1;
      }
      //mprintf("DEBUG: step= %i\n", step);
      ptr = infile.Line(); // Current time (ps)
      float time;
      if (sscanf(ptr,"Time: %f", &time) != 1) {
        mprinterr("Error: Could not get time.\n");
        return 1;
      }
      //mprintf("DEBUG: time= %f\n", time);
      ptr = infile.Line(); // Residue
    } 
    // delta record or full record Residue read
    int res, state;
    float pHval;
    while (sscanf(ptr, rFmt, &res, &state, &pHval) >= 2) {
      //mprintf("DEBUG: res= %i state= %i pH= %f\n", res, state, pHval);
      if (res < maxRes)
        resStates[res] = state;
      else {
        mprinterr("Error: Res %i in CPOUT > max # res in CPIN (%i)\n", res, maxRes);
        return 1;
      }
      ptr = infile.Line();
    }
    phdata->AddState(resStates, pHval);
    nframes++;
    ptr = infile.Line();
  }
  return 0;
}

// DataIO_Cpout::WriteHelp()
void DataIO_Cpout::WriteHelp()
{

}

// DataIO_Cpout::processWriteArgs()
int DataIO_Cpout::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Cpout::WriteHeader()
void DataIO_Cpout::WriteHeader(CpptrajFile& outfile, float solventPH, int stepSize,
                               int step, double time) const
{
  outfile.Printf("Solvent pH: %8.5f\n"
                 "Monte Carlo step size: %8i\n"
                 "Time step: %8i\n"
                 "Time: %10.3f\n", solventPH, stepSize, step, time);
}


// DataIO_Cpout::WriteData()
int DataIO_Cpout::WriteData(FileName const& fname, DataSetList const& dsl)
{
  if (dsl.empty()) return 1;

  DataSet::DataType dtype = dsl[0]->Type();
  if (dtype != DataSet::PH && dtype != DataSet::PH_REMD) {
    mprinterr("Internal Error: Set '%s' is not a pH set.\n", dsl[0]->legend() );
    return 1;
  }

  unsigned int maxFrames = dsl[0]->Size();
  for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds) {
    if ((*ds)->Type() != dtype) {
      mprinterr("Error: Cannot mix sorted and unsorted pH sets.\n");
      return 1;
    }
    if (maxFrames != (*ds)->Size()) {
      mprintf("Warning: Set '%s' frames (%zu) != frames in previous set(s) (%u)\n",
              (*ds)->legend(), (*ds)->Size(), maxFrames);
      maxFrames = std::min( maxFrames, (unsigned int)(*ds)->Size() );
    }
  }

  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open %s for writing.\n", fname.full());
    return 1;
  }

  if (dtype == DataSet::PH_REMD) {
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds) {
      DataSet_pH_REMD const& PH = static_cast<DataSet_pH_REMD const&>( *(*ds) );
      unsigned int idx = 0;
      unsigned int maxres = PH.Residues().size();
      for (unsigned int frame = 0; frame != maxFrames; frame++) {
        if (frame == 0)
          WriteHeader(outfile, PH.pH_Values()[0], 100, 100, 0.198);
        for (unsigned int res = 0; res != maxres; res++, idx++)
          outfile.Printf("Residue %4u State: %2i pH: %7.3f\n",
                         res, PH.ResStates()[idx], PH.pH_Values()[frame]);
        outfile.Printf("\n");
      }
    }
  } else {
    for (unsigned int frame = 0; frame != maxFrames; frame++) {
      if (frame == 0) // TODO check all same pH
        WriteHeader(outfile, ((DataSet_pH*)dsl[0])->Solvent_pH(), 100, 100, 0.198);
      for (unsigned int res = 0; res != dsl.size(); res++)
        outfile.Printf("Residue %4u State: %2i\n", res, ((DataSet_pH*)dsl[res])->State(frame));
      outfile.Printf("\n");
    }
  }

  outfile.CloseFile();
  return 0;
}
