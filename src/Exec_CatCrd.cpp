#include "Exec_CatCrd.h"
#include "CpptrajStdio.h"

// Exec_CatCrd::Help()
void Exec_CatCrd::Help() const
{
  mprintf("\t<set name0> <set name1> [<set nameX> ...] name <name>\n"
          "  Concatenate 2 or more COORDS data sets into a single one.\n");
}

// Exec_CatCrd::Execute()
Exec::RetType Exec_CatCrd::Execute(CpptrajState& State, ArgList& argIn)
{
  // Set up metadata with file name and output set name
  std::string setname = argIn.GetStringKey("name");
  if (setname.empty()) {
    mprinterr("Error: No output COORDS set name specified: 'name <output setname>'\n");
    return CpptrajState::ERR;
  }
  Topology* parm = 0;
  // Get COORDS data sets to concatenate.
  typedef std::vector<DataSet_Coords*> DCarray;
  DCarray inputSets;
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( dsarg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
    {
      if ( (*ds)->Group() != DataSet::COORDINATES )
        mprintf("Warning: Set '%s' is not COORDS, skipping.\n", (*ds)->legend());
      else {
        DataSet_Coords* coordsIn = (DataSet_Coords*)*ds;
        if (parm == 0)
          parm = coordsIn->TopPtr();
        else {
          // Check that topology matches. For now just check # atoms.
          if (parm->Natom() != coordsIn->Top().Natom()) {
            mprinterr("Error: Set '%s' # atoms (%i) differs from first (%i)\n",
                      coordsIn->legend(), coordsIn->Top().Natom(), parm->Natom());
            return CpptrajState::ERR;
          }
        }
        inputSets.push_back( coordsIn );
      }
    }
    dsarg = argIn.GetStringNext();
  }
  if (inputSets.empty()) {
    mprinterr("Error: No input COORDS sets.\n");
    return CpptrajState::ERR;
  }

  // Check if output already present
  DataSet_Coords* coordsOut = 0;
  DataSet* ds = State.DSL().FindSetOfType( setname, DataSet::COORDS );
  if (ds == 0) {
    // Create Set 
    MetaData md( setname );
    coordsOut = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, md);
    if (coordsOut == 0) {
      mprinterr("Error: Could not allocate COORDS data set.\n");
      return CpptrajState::ERR;
    }
    coordsOut->CoordsSetup( *parm, inputSets.front()->CoordsInfo() );
    mprintf("\tNew COORDS set '%s'\n", coordsOut->legend());
  } else {
    // Check that set is actually coords.
    if (ds->Type() != DataSet::COORDS) {
      mprinterr("Error: Set %s present but is not of type COORDS.\n", ds->legend());
      return CpptrajState::ERR;
    }
    coordsOut = (DataSet_Coords*)ds;
    // Check that topology matches. For now just check # atoms.
    if (parm->Natom() != coordsOut->Top().Natom()) {
      mprinterr("Error: # atoms %i does not match COORDS data set '%s' (%i)\n",
                parm->Natom(), coordsOut->legend(), coordsOut->Top().Natom());
      return CpptrajState::ERR;
    }
    mprintf("\tAppending to COORDS data set '%s'\n", coordsOut->legend());
  }

  for (DCarray::const_iterator in = inputSets.begin(); in != inputSets.end(); ++in)
  {
    mprintf("\t'%s'\n", (*in)->legend());
    Frame frameIn = (*in)->AllocateFrame();
    for (unsigned int frm = 0; frm != (*in)->Size(); frm++) {
      (*in)->GetFrame(frm, frameIn);
      coordsOut->AddFrame( frameIn );
    }
  }

  return CpptrajState::OK;
}
