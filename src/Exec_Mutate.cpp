#include "Exec_Mutate.h"
#include "CpptrajStdio.h"
#include "Structure/Creator.h"

// Exec_Mutate::Help()
void Exec_Mutate::Help() const
{
  mprintf("\tcrdset <COORDS set> resmask <mask> [outset <output COORDS>]\n"
          "\t[%s]\n"
          "\t[{%s} ...]\n"
          "\t[{%s} ...]\n",
          Cpptraj::Structure::Creator::other_keywords_,
          Cpptraj::Structure::Creator::template_keywords_,
          Cpptraj::Structure::Creator::parm_keywords_);
}

// Exec_Mutate::Execute()
Exec::RetType Exec_Mutate::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  Cpptraj::Structure::Creator creator( debug_ );
  if (creator.InitCreator(argIn, State.DSL(), debug_)) {
    return CpptrajState::ERR;
  }
  if (!creator.HasTemplates()) {
    mprinterr("Error: No residue templates loaded.\n");
  }

  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: Specify COORDS dataset name with 'crdset'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: No COORDS set with name %s found.\n", setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());

  bool modify_input_coords = true;
  DataSet_Coords* OUT = 0;
  std::string outset = argIn.GetStringKey("outset");
  if (!outset.empty()) {
    OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, outset );
    if (OUT == 0) {
      mprinterr("Error: Could not allocate output set '%s'\n", outset.c_str());
      return CpptrajState::ERR;
    }
    modify_input_coords = false;
    mprintf("\tOutput set: %s\n", OUT->legend());
  } else {
    // Modifying input COORDS set
    if (State.DSL().PopSet( CRD ) == 0) {
      mprinterr("Internal Error: Exec_Mutate::Execute: Could not pop input COORDS set off master DataSetList()\n");
      return CpptrajState::ERR;
    }
    OUT = (DataSet_Coords*)State.DSL().AddSet( DataSet::COORDS, CRD->Meta() );
    mprintf("\tWill modify COORDS set %s\n", CRD->legend());
  }

  CpptrajState::RetType ret = doMutate( State, argIn, CRD, OUT, creator );

  if (modify_input_coords && CRD != 0) delete CRD;

  return ret;
}

/** Actually do the mutation(s) */
CpptrajState::RetType Exec_Mutate::doMutate(CpptrajState& State, ArgList& argIn, DataSet_Coords* CRD, DataSet_Coords* OUT,
                                           Cpptraj::Structure::Creator const& creator)
const
{

  std::string resmask = argIn.GetStringKey("resmask");
  if (resmask.empty()) {
    mprinterr("Error: Specify mask of residues to mutate with 'resmask'\n");
    return CpptrajState::ERR;
  }

  std::string templateName = argIn.GetStringKey("to");
  if (templateName.empty()) {
    mprinterr("Error: Specify template name to mutate to with 'to'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* UNIT = creator.IdTemplateFromName( templateName );
  if (UNIT == 0) {
    mprinterr("Error: Could not get template for '%s'\n", templateName.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tMutate residues selected by '%s' to '%s'\n", resmask.c_str(), UNIT->legend());

  AtomMask mask;
  if (mask.SetMaskString( resmask )) {
    mprinterr("Error: Could not set mask '%s'\n", resmask.c_str());
    return CpptrajState::ERR;
  }
  if (CRD->Top().SetupIntegerMask( mask )) {
    mprinterr("Error: Could not setup mask '%s'\n", mask.MaskString());
    return CpptrajState::ERR;
  }
  if (mask.None()) {
    mprinterr("Error: Nothing selected by mask '%s'\n", mask.MaskString());
    return CpptrajState::ERR;
  }
  //mask.MaskInfo();
  std::vector<int> resnums = CRD->Top().ResnumsSelectedBy( mask );
  mprintf("\t%zu residues selected by '%s'\n", resnums.size(), mask.MaskString());

  AtomMask toKeep;
  toKeep.SetNatoms( CRD->Top().Natom() );
  std::vector<NameType> SourceAtomNames;
  SourceAtomNames.resize( CRD->Top().Natom() );
  std::vector<int>::const_iterator rnum = resnums.begin();
  for (int ires = 0; ires != CRD->Top().Nres(); ires++)
  {
    // Is this a selected res?
    if (rnum != resnums.end() && ires == *rnum) {
      // Selected res. Keep only atoms present in template
      ++rnum;
      int nTgtAtomsMissing = 0;
      std::vector<int> templateToRes = creator.MapAtomsToTemplate( CRD->Top(), ires, UNIT, SourceAtomNames, nTgtAtomsMissing );
      //if (debug_ > 1) {
        mprintf("\tResidue %i Atom map:\n", ires + 1);
        // DEBUG - print map
        for (int iref = 0; iref != UNIT->Top().Natom(); iref++) {
          mprintf("\t\t%6i %6s =>", iref+1, *(UNIT->Top()[iref].Name()));
          if (templateToRes[iref] == -1)
            mprintf(" No match\n");
          else
            mprintf(" %6i %6s\n", templateToRes[iref]+1, *(CRD->Top()[templateToRes[iref]].Name()));
        }
      //}
      // For each template atom, only keep what was mapped
      for (int iref = 0; iref != UNIT->Top().Natom(); iref++)
        if (templateToRes[iref] > -1)
          toKeep.AddSelectedAtom( templateToRes[iref] );
    } else {
      // Not a selected res. Keep all atoms
      Residue const& currentRes = CRD->Top().Res( ires );
      for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); ++at)
        toKeep.AddSelectedAtom( at );
    }
  }

/*
  for (std::vector<int>::const_iterator rnum = resnums.begin(); rnum != resnums.end(); ++rnum)
  {
    Residue const& currentRes = CRD->Top().Res( *rnum );
    int atomsToRemove = 0;
    for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); at++)
    {
      // Does this atom exist in the template?
      int idx = UNIT->Top().FindAtomInResidue(0, CRD->Top()[at].Name());
      if (idx < 0) {
        toRemove.AddSelectedAtom( at );
        atomsToRemove++;
      }
      if (idx > -1)
        mprintf("DEBUG: Found atom %s in template.\n", CRD->Top().AtomMaskName(at).c_str());
      else
        mprintf("DEBUG: Atom %s not in template.\n", CRD->Top().AtomMaskName(at).c_str());
    }
    mprintf("DEBUG: Removing %i atoms from residue %i\n", atomsToRemove, *rnum + 1 );
    if (currentRes.NumAtoms() - atomsToRemove < 1) {
      mprinterr("Error: Number of atoms to remove %i >= number of atoms in residue %s (%i)\n",
                atomsToRemove, CRD->Top().TruncResNameNum(*rnum).c_str(), currentRes.NumAtoms());
      return CpptrajState::ERR;
    }
  }
  toRemove.InvertMask();
  Topology* newTop = CRD->Top().modifyStateByMask( toRemove );
*/
  Topology* newTop = CRD->Top().modifyStateByMask( toKeep );
  if (newTop == 0) {
    mprinterr("Error: Could not remove atoms from '%s'\n", CRD->legend());
    return CpptrajState::ERR;
  }
  newTop->Summary();
  // Change the residue names
  for (std::vector<int>::const_iterator rnum = resnums.begin(); rnum != resnums.end(); ++rnum)
    newTop->SetRes( *rnum ).SetName( templateName );

  // Set up output coords
  OUT->CoordsSetup( *newTop, CRD->CoordsInfo() );
  Frame newFrame;
  newFrame.SetupFrameV(newTop->Atoms(), CRD->CoordsInfo());

  // Strip all input coords frames
  Frame inputFrame = CRD->AllocateFrame();
  for (unsigned int frm = 0; frm != CRD->Size(); ++frm)
  {
    CRD->GetFrame(frm, inputFrame);

    newFrame.SetFrame(inputFrame, toKeep);
    OUT->AddFrame(newFrame);
  }

  if (newTop != 0) delete newTop;

  return CpptrajState::OK;
}
