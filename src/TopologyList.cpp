// ParmList
#include "TopologyList.h"
#include "CpptrajStdio.h"
#include "AtomMask.h"
#include "ParmFile.h"

// CONSTRUCTOR 
TopologyList::TopologyList() : 
  hasCopies_(false)
{}

// DESTRUCTOR
TopologyList::~TopologyList() {
  if (!hasCopies_) {
    for (std::vector<Topology*>::iterator top = TopList_.begin();
                                          top != TopList_.end(); top++)
      delete *top;
  }
}

void TopologyList::Help_Parm() {
  mprintf("parm <filename> [<tag>] [nobondsearch | bondsearch [<offset>]]\n");
  mprintf("\tAdd <filename> to parm list\n");
}

void TopologyList::Help_ParmInfo() {
  mprintf("parminfo [<parmindex>] [<mask>]:\n");
  mprintf("\tPrint information on parm <parmindex> (0 by default). If <mask> is given\n");
  mprintf("print info on atoms in mask. If no mask given print overall information.\n");
}

void TopologyList::Help_ParmWrite() {
  mprintf("parmwrite out <filename> [<parmindex>]\n");
  mprintf("\tWrite parm <parmindex> to <filename>\n");
}

void TopologyList::Help_ParmStrip() {
  mprintf("parmstrip <mask> [<parmindex>]\n");
  mprintf("\tStrip atoms in mask from parm\n");
}

void TopologyList::Help_ParmBox() {
  mprintf("[parm]box [<parmindex>] [x <xval>] [y <yval>] [z <zval>]");
  mprintf(" [alpha <a>] [beta <b>] [gamma <g>] [nobox]\n");
  mprintf("\tSet the given parm box info to what is specified. If nobox, remove box info.\n");
}

void TopologyList::Help_Solvent() {
  mprintf("solvent [<parmindex>] <mask>\n");
  mprintf("\tSet solvent for the given parm (default 0) based on <mask>\n");
}

void TopologyList::Help_BondInfo() {
  mprintf("parmbondinfo [<parmindex>]\n");
  mprintf("\tPrint bond information for parm <parmindex> (0 by default).\n");
}

void TopologyList::Help_ResInfo() {
  mprintf("parmresinfo [<parmindex>]\n");
  mprintf("\tPrint residue information for parm <parmindex> (0 by default).\n");
}

void TopologyList::Help_MolInfo() {
  mprintf("parmmolinfo [<parmindex>]\n");
  mprintf("\tPrint molecule information for parm <parmindex> (0 by default).\n");
}

enum ParmCmdTypes { LOADPARM=0, PARMINFO, PARMWRITE, PARMSTRIP, PARMBOX,
                    SOLVENT, BONDINFO, RESINFO, MOLINFO, DEPRECATED,
                    SELECT };
// TODO: Make deprecated a separate list
const DispatchObject::Token TopologyList::ParmCmds[] = {
  { DispatchObject::PARM, "box", 0, Help_ParmBox, PARMBOX },
  { DispatchObject::PARM, "parm", 0, Help_Parm, LOADPARM },
  { DispatchObject::PARM, "parmbondinfo", 0, Help_BondInfo, BONDINFO },
  { DispatchObject::PARM, "parmbox", 0, Help_ParmBox, PARMBOX },
  { DispatchObject::PARM, "parminfo", 0, Help_ParmInfo, PARMINFO },
  { DispatchObject::PARM, "parmmolinfo", 0, Help_MolInfo, MOLINFO },
  { DispatchObject::PARM, "parmresinfo", 0, Help_ResInfo, RESINFO },
  { DispatchObject::PARM, "parmstrip", 0, Help_ParmStrip, PARMSTRIP },
  { DispatchObject::PARM, "parmwrite", 0, Help_ParmWrite, PARMWRITE },
  { DispatchObject::PARM, "select", 0, 0, SELECT },
  { DispatchObject::PARM, "solvent", 0, Help_Solvent, SOLVENT },
  { DispatchObject::PARM, "molsearch", 0, 0, DEPRECATED },
  { DispatchObject::PARM, "nomolsearch", 0, 0, DEPRECATED },
  { DispatchObject::PARM, "bondsearch", 0, 0, DEPRECATED },
  { DispatchObject::PARM, "nobondsearch", 0, 0, DEPRECATED },
  { DispatchObject::NONE,                  0, 0,                 0, 0 }
};
  
int TopologyList::LoadParm(ArgList& argIn) {
  std::string parmtag = argIn.getNextTag();
  bool bondsearch = !argIn.hasKey("nobondsearch");
  double offset = argIn.getKeyDouble("bondsearch", -1.0);
  return AddParmFile(argIn.GetStringNext(), parmtag, bondsearch, offset);
}

int TopologyList::ParmInfo(ArgList& argIn) {
  int pindex = argIn.getNextInteger(0);
  if (pindex>=0 && pindex<(int)TopList_.size()) {
    std::string maskarg = argIn.GetMaskNext();
    if (!maskarg.empty()) 
      TopList_[pindex]->PrintAtomInfo( maskarg );
    else 
      TopList_[pindex]->Summary();
  } else {
    mprinterr("Error: parminfo: parm index %i not loaded.\n",pindex);
    return 1;
  }
  return 0;
}

int TopologyList::ParmWrite(ArgList& argIn) {
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: parmwrite: No output filename specified (use 'out <filename>').\n");
    return 1;
  }
  int pindex = argIn.getNextInteger(0);
  if (pindex < 0 || pindex >= (int)TopList_.size()) {
    mprinterr("Error: parmwrite: parm index %i out of bounds.\n",pindex);
    return 1;
  }
  mprintf("\tWriting parm %i (%s) to Amber parm %s\n",pindex,
          TopList_[pindex]->c_str(), outfilename.c_str());
  ParmFile pfile;
  pfile.Write( *TopList_[pindex], outfilename, ParmFile::AMBERPARM, debug_ );
  return 0;
}

int TopologyList::ParmStrip(ArgList& argIn) {
  int pindex = argIn.getNextInteger(0);
  if (pindex < 0 || pindex >= (int)TopList_.size()) {
    mprinterr("Error: parmstrip: parm index %i out of bounds.\n",pindex);
    return 1;
  }
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMask();
  TopList_[pindex]->SetupIntegerMask( tempMask );
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(), 
           TopList_[pindex]->Natom() - tempMask.Nselected(), TopList_[pindex]->c_str());
  Topology* tempParm = TopList_[pindex]->modifyStateByMask(tempMask);
  if (tempParm==0) { 
    mprinterr("Error: parmstrip: Could not strip parm.\n");
    return 1;
  } else {
    // Replace parm with stripped version
    tempParm->ParmInfo();
    if (!hasCopies_) delete TopList_[pindex];
    TopList_[pindex] = tempParm;
  }
  return 0;
}

int TopologyList::ParmBox(ArgList& argIn) {
  Box pbox;
  pbox.SetX( argIn.getKeyDouble("x",0) );
  pbox.SetY( argIn.getKeyDouble("y",0) );
  pbox.SetZ( argIn.getKeyDouble("z",0) );
  pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
  pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
  pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
  bool nobox = argIn.hasKey("nobox"); 
  // Get parm index
  int pindex = argIn.getNextInteger(0);
  if (pindex < 0 || pindex >= (int)TopList_.size()) {
    mprinterr("Error: box: parm index %i out of bounds.\n",pindex);
    return 1;
  }
  if (nobox)
    TopList_[pindex]->SetBox( Box() );
  else {
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( TopList_[pindex]->ParmBox() );
    TopList_[pindex]->SetBox( pbox );
  }
  return 0;
}

int TopologyList::ParmSolvent(ArgList& argIn) {
  std::string maskexpr = argIn.GetMaskNext();
  if ( maskexpr.empty() ) {
    mprinterr("Error: solvent: No mask specified.\n");
    return 1;
  }
  // Get parm index
  int pindex = argIn.getNextInteger(0);
  if (pindex < 0 || pindex >= (int)TopList_.size()) {
    mprinterr("Error: solvent: parm index %i out of bounds.\n",pindex);
    return 1;
  }
  TopList_[pindex]->SetSolvent( maskexpr );
  return 0;
}

int TopologyList::Select(ArgList& argIn) {
  AtomMask tempMask( argIn.GetMaskNext() );
  int pindex = argIn.getNextInteger(0);
  if (pindex < 0 || pindex >= (int)TopList_.size()) {
    mprinterr("Error: select: parm index %i out of bounds.\n",pindex);
    return 1;
  }
  TopList_[pindex]->SetupIntegerMask( tempMask );
  tempMask.PrintMaskAtoms("Selected");
  return 0;
}
 
// TopologyList::CheckCommand()
/** Check if the command in the arglist pertains to topology files.
  * \return 0 if command was recognized, 1 if not.
  */
int TopologyList::CheckCommand(int cmdidxIn, ArgList& argIn) {
  int err = 0;
  ParmCmdTypes cmdidx = (ParmCmdTypes) cmdidxIn;
  switch ( cmdidx ) {
    case LOADPARM: err = LoadParm( argIn ); break;
    case PARMINFO: err = ParmInfo( argIn ); break;
    case PARMWRITE: err = ParmWrite( argIn ); break;
    case PARMSTRIP: err = ParmStrip( argIn ); break;
    case PARMBOX: err = ParmBox( argIn ); break;
    case SOLVENT: err = ParmSolvent(argIn); break;
    case SELECT: err = Select(argIn); break;
    case DEPRECATED:
      mprintf("Warning: %s is deprecated.\n", argIn.Command());
      break;
    default: err = -1; // One of the parmXinfo commands
  }
  if (err == -1) {
    int pindex = argIn.getNextInteger(0);
    if (pindex < 0 || pindex >= (int)TopList_.size()) {
      mprinterr("Error: %s: parm %i not loaded.\n",argIn.Command(), pindex);
      return 1;
    }
    err = 0;
    switch ( cmdidx ) {
      case BONDINFO: TopList_[pindex]->PrintBondInfo(argIn.GetMaskNext()); break;
      case RESINFO : TopList_[pindex]->PrintResidueInfo(); break;
      case MOLINFO : TopList_[pindex]->PrintMoleculeInfo(argIn.GetMaskNext()); break;
      default: err = 1; // Should never get here
    }
  }
  return err;
}

// TopologyList::GetParm()
/** Return the parm structure with index num. */
Topology *TopologyList::GetParm(int num) {
  if (num>=(int)TopList_.size() || num<0) return NULL;
  return TopList_[num];
}

// TopologyList::GetParm()
/** Return the parm structure based on arguments in the given arg list. 
  *   parm <parm name>
  *   parmindex <parm index>
  * \param argIn argument list that contains parm-related keyword
  * \return parm specified by 'parm' or 'parmindex', or the first parm. NULL on error.
  */
Topology *TopologyList::GetParm(ArgList &argIn) {
  // Get any parm keywords if present
  std::string parmfilename = argIn.GetStringKey("parm");
  int pindex = argIn.getKeyInt("parmindex",0);
  // Associate trajectory with parameter file. Associate with default parm if none specified
  if (!parmfilename.empty())
    pindex = FindName(parmfilename);
  Topology *ParmOut = GetParm(pindex);
  if (ParmOut==NULL) {
    mprintf("Warning: Could not get parameter file:\n");
    mprintf("Warning: parmname=%s, pindex=%i\n",parmfilename.c_str(),pindex);
    return NULL;
  }

  return ParmOut;
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list. */
int TopologyList::AddParmFile(std::string const& filename) {
  std::string emptystring;
  return AddParmFile(filename, emptystring, true, -1.0);
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list with optional tag. */
int TopologyList::AddParmFile(std::string const& filename, std::string const& ParmTag,
                              bool bondsearch, double offset) 
{
  // Dont let a list that has copies add a new file
  if (hasCopies_) {
    mprintf("Warning: Attempting to add parm %s to a list that already\n",filename.c_str());
    mprintf("Warning: has copies of parm files. This should not occur.\n");
    mprintf("Warning: Skipping.\n");
    return 0;
  }

  // Check if this file has already been loaded
  if (FindName(filename)!=-1) {
    mprintf("Warning: Parm %s already loaded, skipping.\n",filename.c_str());
    return 0;
  }

  // If tag specified, check if tag already in use
  if (FindName(ParmTag)!=-1) {
    mprinterr("Error: Parm tag [%s] already in use.\n",ParmTag.c_str());
    return 1;
  }

  Topology *parm = new Topology();
  parm->SetDebug( debug_ );
  parm->SetOffset( offset );
  ParmFile pfile;
  int err = pfile.Read(*parm, filename, bondsearch, debug_);
  if (err!=0) {
    mprinterr("Error: Could not open parm %s\n",filename.c_str());
    delete parm;
    return 1;
  }

  if (debug_>0) 
    mprintf("    PARAMETER FILE %zu: %s\n",TopList_.size(),filename.c_str());
  // pindex is used for quick identification of the parm file
  parm->SetPindex( TopList_.size() );
  TopList_.push_back(parm);
  AddNameWithTag( filename, pfile.BaseName(), ParmTag);
  return 0;
}

// TopologyList::AddParm()
/** Add an existing AmberParm to parm file list. Currently used to keep track
  * of parm files corresponding to frames in the reference frame list.
  */
int TopologyList::AddParm(Topology *ParmIn) {
  if (ParmIn==NULL) return 1;
  if (!hasCopies_ && !TopList_.empty()) {
    mprinterr("Error: Attempting to add copy of parm to list with non-copies!\n");
    return 1;
  }
  // Set the hasCopies flag so we know not to try and delete these parms
  hasCopies_=true;
  //P->pindex=Nparm; // pindex should already be set
  TopList_.push_back(ParmIn);
  return 0;
}

// TopologyList::List()
/** Print list of loaded parameter files */
void TopologyList::List() {
  mprintf("\nPARAMETER FILES:\n");
  if (TopList_.empty()) {
    mprintf("  No parameter files defined.\n");
    return;
  }
  for (std::vector<Topology*>::iterator top = TopList_.begin(); top != TopList_.end(); top++)
    (*top)->ParmInfo();

    //mprintf("  %i: %s, %i atoms (%i trajectory frames associated)\n",
    //        i,ParmList[i]->File.filename, ParmList[i]->natom, ParmList[i]->parmFrames);
}
