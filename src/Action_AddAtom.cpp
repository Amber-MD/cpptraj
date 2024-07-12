#include "Action_AddAtom.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Action_AddAtom::Action_AddAtom() :
  newParm_(0)
{}

/** DESTRUCTOR */
Action_AddAtom::~Action_AddAtom() {
  if (newParm_ != 0) delete newParm_;
}

// Action_AddAtom::Help()
void Action_AddAtom::Help() const {

}

// Action_AddAtom::Init()
Action::RetType Action_AddAtom::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get output stripped parm filename
  topWriter_.InitTopWriter(actionArgs, "stripped", debugIn);

  std::string aname = actionArgs.GetStringKey("aname");
  if (aname.empty()) {
    mprinterr("Error: Must specify atom name with 'aname'.\n");
    return Action::ERR;
  }
  NameType atomName( aname );

  std::string elt = actionArgs.GetStringKey("elt");
  if (elt.empty())
    elt.assign("H");
  if (elt.size() > 2) {
    mprinterr("Error: Element name '%s' is too big; should be 2 characters max.\n", elt.c_str());
    return Action::ERR;
  }

  std::string rname = actionArgs.GetStringKey("rname");
  if (rname.empty())
    rname.assign("TMP");
  residueName_ = NameType( rname );

  newAtom_ = Atom(atomName, elt.c_str());

  mprintf("    ADDATOM: Adding atom named '%s', element %s, residue name '%s'\n",
          *atomName, elt.c_str(), *residueName_);
  topWriter_.PrintOptions();

  return Action::OK;
}

// Action_AddAtom::Setup()
Action::RetType Action_AddAtom::Setup(ActionSetup& setup)
{

}

// Action_AddAtom::DoAction()
Action::RetType Action_AddAtom::DoAction(int frameNum, ActionFrame& frm)
{

}
