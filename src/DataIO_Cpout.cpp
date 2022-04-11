#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include "DataIO_Cpout.h"
#include "CpptrajStdio.h"
#include "DataSet_pH.h"
#include "DataSet_PHREMD_Explicit.h"
#include "DataSet_PHREMD_Implicit.h"

/// CONSTRUCTOR
DataIO_Cpout::DataIO_Cpout() :
  type_(NONE),
  original_pH_(0.0),
  maxRes_(0),
  nframes_(0),
  recType_(0),
  mc_stepsize_(0),
  step_(0),
  s0_(0),
  nRes_(0), // FIXME may not need to be a class var
  solvent_pH_(0.0),
  pHval_(0.0),
  time_(0.0),
  t0_(0.0)
{
  SetValid( DataSet::PH );
  SetValid( DataSet::PH_EXPL );
  SetValid( DataSet::PH_IMPL );
}

const char* DataIO_Cpout::FMT_REDOX_ = "Redox potential: %f V";
//const char* DataIO_Cpout::FMT_REDOX_ = "Redox potential: %f V Temperature: %f K";

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
      Cph::CpRes(split[0], atoi(split[1].c_str()), res_protcnt, max_prots) );
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
  if (ReadCpin( cpin_file_ )) {
    mprinterr("Error: Could not read CPIN file '%s'\n", cpin_file_.full());
    return 1;
  }

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
  // Try to determine the number of residues
  nRes_ = 0;
  while (sscanf(ptr, rFmt, &res, &state, &pHval) >= 2) {
    nRes_++;
    ptr = infile.Line();
  }
  mprintf("\t%i residues in first record.\n", nRes_);
  maxRes_ = nRes_;
  resStates_.resize( maxRes_, 0 );
  // If unsorted data, attempt to determine if from implicit solvent
  // (delta records have info for only one residue per step). Means
  // less data needs to be stored and sorting needs to happen differently.
  bool isImplicit = false;
  if (nscan == 3) {
    // Look ahead 3 records.
    for (int nframe = 1; nframe < 3; nframe++)
    {
      state = ReadRecord(infile, fmt, rFmt);
      if (state == -1) return 1;
      if (state == 1 && recType_ >= 0) {
        // Unsorted implicit
        isImplicit = true;
        break;
      }
    }
  }
  if (isImplicit)
    mprintf("\tUnsorted implicit pH data detected.\n");
  infile.CloseFile();
  if (infile.OpenFileRead( fname )) return 1;

  // Allocate DataSets
  t0_ = -1.0;
  s0_ = -1;
  nframes_ = 0;
  int err = 1;
  if (nscan == 2) {
    // Sorted constant pH
    err = ReadSorted(infile, dsl, dsname, fmt, rFmt);
  } else if (nscan == 3) {
    // Unsorted constant pH
    if (isImplicit)
      err = ReadUnsortedImplicit(infile, dsl, dsname, fmt, rFmt);
    else
      err = ReadUnsortedExplicit(infile, dsl, dsname, fmt, rFmt);
  }
  infile.CloseFile();
/*
  if (debug_ > 1) {
    for (DataSet_pH_REMD::const_iterator res = phdata->begin(); res != phdata->end(); ++res) {
      mprintf("DEBUG: Res %u:\n", res-phdata->begin());
      for (Cph::Res::const_iterator state = res->begin();
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

/** Read a single constant pH data record. For full records sets solvent_pH_,
  * mc_stepsize_, step_, s0_, and time_. For full and delta records sets
  * recType_, resStates_, and pHval_.
  */
int DataIO_Cpout::ReadRecord(BufferedLine& infile, const char* fmt, const char* rFmt) {
  const char* ptr = infile.Line();
  //mprintf("DEBUG: Record: '%s'\n", ptr);
  if (ptr == 0) return 0;
  recType_ = Cph::PARTIAL_RECORD;
  if (sscanf(ptr, fmt, &solvent_pH_) == 1) {
    // Full record
    recType_ = Cph::FULL_RECORD;
    //mprintf("DEBUG: pH= %f\n", solvent_pH_);
    // Monte Carlo step size - should never change
    ptr = infile.Line();
    sscanf(ptr, "Monte Carlo step size: %i", &mc_stepsize_);
    // Current MD time step
    ptr = infile.Line();
    if (sscanf(ptr,"Time step: %d", &step_) != 1) {
      mprinterr("Error: Could not get step.\n");
      return -1;
    }
    if (s0_ < 0) s0_ = step_;
    //mprintf("DEBUG: step= %i\n", step_);
    // Current time (ps)
    ptr = infile.Line();
    if (sscanf(ptr, "Time: %lf", &time_) != 1) {
      mprinterr("Error: Could not get time.\n");
      return -1;
    }
    if (t0_ < 0.0) t0_ = time_;
    //mprintf("DEBUG: time= %f\n", time_);
    ptr = infile.Line(); // Residue
  }
  // delta record or full record Residue read
  int res, state;
  pHval_ = solvent_pH_;
  nRes_ = 0;
  //mprintf("DEBUG: Res: '%s'\n", ptr);
  while (sscanf(ptr, rFmt, &res, &state, &pHval_) >= 2) {
    //mprintf("DEBUG: res= %i state= %i pH= %f\n", res, state, pHval_);
    //mprintf("DEBUG: res= %i state= %i\n", res, state);
    if (res < maxRes_)
      resStates_[res] = state;
    else {
      mprinterr("Error: Res %i in CPOUT > max # res in CPIN (%i)\n", res, maxRes_);
      return -1;
    }
    nRes_++;
    ptr = infile.Line();
    if (ptr == 0) break;
  }
  if (nRes_ == 1)
    recType_ = res;
  else if (nRes_ < maxRes_) {
    mprinterr("Error: Only read %i residues - expected %i\n", nRes_, maxRes_);
    return -1;
  }
  //mprintf("DEBUG: %6.2f", pHval_);
  //for (Iarray::const_iterator it = resStates_.begin(); it != resStates_.end(); ++it)
  //  mprintf(" %2i", *it);
  //mprintf("\n");
  nframes_++;
  return 1;
};

/** Calculate time step based on initial and final time and step count. */
double DataIO_Cpout::CalcTimeStep() const {
  double dt = (time_ - t0_) / ((double)(step_ - s0_));
  mprintf("\tMC step size %i, t0 = %.3f, tf = %.3f, nframes= %u, dt = %.3f\n",
          mc_stepsize_, t0_, time_, nframes_, dt);
  return dt;
}

/** Read sorted pH data. Will generate a single DataSet for each residue. */
int DataIO_Cpout::ReadSorted(BufferedLine& infile, DataSetList& DSL, std::string const& dsname, const char* fmt, const char* rFmt)
{
  // Create data sets for each residue 
  typedef std::vector<DataSet_pH*> Parray;
  Parray ResSets;
  ResSets.reserve( Residues_.size() );
  for (Rarray::iterator res = Residues_.begin(); res != Residues_.end(); ++res)
  {
    MetaData md( dsname, res->Name().Truncated(), res->Num() );
    DataSet* ds = DSL.CheckForSet(md);
    if (ds == 0) {
      // New set
      ds = DSL.AddSet( DataSet::PH, md );
      if (ds == 0) return 1;
      ((DataSet_pH*)ds)->SetResidueInfo( *res );
      ((DataSet_pH*)ds)->Set_Solvent_pH( original_pH_ ); // TODO combine with above?
    } else {
      // TODO may need to skip reading first record on append
      if (ds->Type() != DataSet::PH) {
        mprinterr("Error: Set '%s' type does not match, cannot append.\n", ds->legend());
        return 1;
      }
      mprintf("\tAppending to set '%s'\n", ds->legend());
      // TODO check # residues etc?
    }
    ResSets.push_back( (DataSet_pH*)ds );
  }
  // Loop over records
  while ( ReadRecord(infile, fmt, rFmt) == 1 )
  {
    for (unsigned int idx = 0; idx < resStates_.size(); idx++)
      ResSets[idx]->AddState( resStates_[idx], recType_ );
  }
  double dt = CalcTimeStep();
  Dimension xdim(t0_, dt, "Time");
  for (Parray::iterator p = ResSets.begin(); p != ResSets.end(); ++p)
  {
    (*p)->SetTimeValues(Cph::CpTime(mc_stepsize_, t0_, dt));
    (*p)->SetDim(Dimension::X, xdim);
  }
  return 0;
}

/** Read unsorted explicit pH data (from e.g. replica exchange). Will
  * create a single DataSet containing pH values and states for
  * all residues.
  */
int DataIO_Cpout::ReadUnsortedExplicit(BufferedLine& infile, DataSetList& DSL, std::string const& dsname, const char* fmt, const char* rFmt)
{
  // Unsorted constant pH
  DataSet* ds = DSL.CheckForSet(dsname);
  if (ds == 0) {
    // New set
    ds = DSL.AddSet( DataSet::PH_EXPL, dsname, "ph" );
    if (ds == 0) return 1;
    ((DataSet_PHREMD_Explicit*)ds)->SetResidueInfo( Residues_ );
  } else {
    // TODO may need to skip reading first record on append
    if (ds->Type() != DataSet::PH_EXPL) {
      mprinterr("Error: Set '%s' is not unsorted explicit pH data.\n", ds->legend());
      return 1;
    }
    mprintf("\tAppending to set '%s'\n", ds->legend());
    // TODO check # residues etc?
  }
  DataSet_PHREMD_Explicit* phdata = (DataSet_PHREMD_Explicit*)ds;

  //float solvent_pH = original_pH_;
  while ( ReadRecord(infile, fmt, rFmt) == 1 )
    phdata->AddState(resStates_, pHval_, recType_);
  phdata->SetTimeValues(Cph::CpTime(mc_stepsize_, t0_, CalcTimeStep()));
  return 0;
}

/** Read unsorted implicit pH data (from e.g. replica exchange). Will
  * create a single DataSet containing pH values and states for
  * all residues.
  */
int DataIO_Cpout::ReadUnsortedImplicit(BufferedLine& infile, DataSetList& DSL, std::string const& dsname, const char* fmt, const char* rFmt)
{
  // Unsorted constant pH
  DataSet* ds = DSL.CheckForSet(dsname);
  if (ds == 0) {
    // New set
    ds = DSL.AddSet( DataSet::PH_IMPL, dsname, "ph" );
    if (ds == 0) return 1;
    ((DataSet_PHREMD_Implicit*)ds)->SetResidueInfo( Residues_ );
  } else {
    // TODO may need to skip reading first record on append
    if (ds->Type() != DataSet::PH_IMPL) {
      mprinterr("Error: Set '%s' is not unsorted implicit pH data.\n", ds->legend());
      return 1;
    }
    mprintf("\tAppending to set '%s'\n", ds->legend());
    // TODO check # residues etc?
  }
  DataSet_PHREMD_Implicit* phdata = (DataSet_PHREMD_Implicit*)ds;

  //float solvent_pH = original_pH_;
  while ( ReadRecord(infile, fmt, rFmt) == 1 )
    phdata->AddRecord(DataSet_PHREMD_Implicit::Record(pHval_, recType_, resStates_));
  phdata->SetTimeValues(Cph::CpTime(mc_stepsize_, t0_, CalcTimeStep()));
  return 0;
}

// =============================================================================
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
void DataIO_Cpout::WriteHeader(CpptrajFile& outfile, double time0, double dt, float solventPH, int frame) const
{
  int time_step = (frame+1)*mc_stepsize_;
  double time = time0 + ((double)(time_step - mc_stepsize_) * dt);
  outfile.Printf("Solvent pH: %8.5f\n"
                 "Monte Carlo step size: %8i\n"
                 "Time step: %8i\n"
                 "Time: %10.3f\n", solventPH, mc_stepsize_, time_step, time);
}

// DataIO_Cpout::WriteData()
int DataIO_Cpout::WriteData(FileName const& fname, DataSetList const& dsl)
{
  if (dsl.empty()) return 1;

  DataSet::DataType dtype = dsl[0]->Type();
  if (dtype != DataSet::PH && dsl[0]->Group() != DataSet::PHREMD) {
    mprinterr("Internal Error: Set '%s' is not a pH set.\n", dsl[0]->legend() );
    return 1;
  }

  unsigned int maxFrames = dsl[0]->Size();
  for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds) {
    if ((*ds)->Type() != dtype) {
      mprinterr("Error: Cannot mix different pH set types.\n");
      return 1;
    }
    // For sorted sets need to have same number of frames
    if (dtype == DataSet::PH) {
      if (maxFrames != (*ds)->Size()) {
        mprintf("Warning: Set '%s' frames (%zu) != frames in previous set(s) (%u)\n",
                (*ds)->legend(), (*ds)->Size(), maxFrames);
        maxFrames = std::min( maxFrames, (unsigned int)(*ds)->Size() );
      }
    }
  }

  mprintf("\tWriting %u frames\n", maxFrames);
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open %s for writing.\n", fname.full());
    return 1;
  }

  if (dtype == DataSet::PH_EXPL) {
    // Unsorted explicit pH - all complete residue records.
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds) {
      DataSet_PHREMD_Explicit const& PH = static_cast<DataSet_PHREMD_Explicit const&>( *(*ds) );
      unsigned int idx = 0;
      unsigned int maxres = PH.Residues().size();
      mc_stepsize_ = PH.Time().MonteCarloStepSize();
      for (unsigned int frame = 0; frame != maxFrames; frame++) {
        if (PH.RecordType(frame) == Cph::FULL_RECORD)
          WriteHeader(outfile, PH.Time().InitialTime(), PH.Time().TimeStep(), PH.pH_Values()[frame], frame);
        for (unsigned int res = 0; res != maxres; res++, idx++)
          outfile.Printf("Residue %4u State: %2i pH: %7.3f\n",
                         res, PH.ResStates()[idx], PH.pH_Values()[frame]);
        outfile.Printf("\n");
      }
    }
  } else if (dtype == DataSet::PH_IMPL) {
    // Unsorted implicit pH - mix of complete and single residue records.
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds) {
      DataSet_PHREMD_Implicit const& PH = static_cast<DataSet_PHREMD_Implicit const&>( *(*ds) );
      unsigned int maxres = PH.Residues().size();
      mc_stepsize_ = PH.Time().MonteCarloStepSize();
      for (unsigned int frame = 0; frame != maxFrames; frame++) {
        DataSet_PHREMD_Implicit::Record const& Rec = PH.Records()[frame];
        int rectype = Rec.RecType();
        if ( rectype < 0 ) {
          if (rectype == Cph::FULL_RECORD)
            WriteHeader(outfile, PH.Time().InitialTime(), PH.Time().TimeStep(), Rec.pH(), frame);
          for (unsigned int res = 0; res != maxres; res++)
            outfile.Printf("Residue %4u State: %2i pH: %7.3f\n",
                           res, Rec.ResStates()[res], Rec.pH());
          outfile.Printf("\n");
        } else {
            outfile.Printf("Residue %4u State: %2i pH: %7.3f\n\n",
                           rectype, Rec.ResStates()[0], Rec.pH());
        }
      }
    }
  } else {
    // TODO Check that all are at the same pH and have same time values.
    DataSet_pH* firstSet = ((DataSet_pH*)dsl[0]);
    float solventPH = firstSet->Solvent_pH();
    mc_stepsize_ = firstSet->Time().MonteCarloStepSize();
    for (unsigned int frame = 0; frame != maxFrames; frame++) {
      int rectype = firstSet->RecordType(frame);
      if ( rectype < 0 ) {
        if ( rectype == Cph::FULL_RECORD)
          WriteHeader(outfile, firstSet->Time().InitialTime(), firstSet->Time().TimeStep(), solventPH, frame);
        for (unsigned int res = 0; res != dsl.size(); res++)
          outfile.Printf("Residue %4u State: %2i\n", res, ((DataSet_pH*)dsl[res])->State(frame));
        outfile.Printf("\n");
      } else
        outfile.Printf("Residue %4u State: %2i\n\n",
                       rectype, ((DataSet_pH*)dsl[rectype])->State(frame));
    }
  }

  outfile.CloseFile();
  return 0;
}
