// ParmList
#include "TopologyList.h"
#include "CpptrajStdio.h"
#include "AtomMask.h"
#include "ParmFile.h"

// CONSTRUCTOR 
TopologyList::TopologyList() : 
  hasCopies_(false),
  bondsearch_(true),
  molsearch_(true)
{}

// DESTRUCTOR
TopologyList::~TopologyList() {
  if (!hasCopies_) {
    for (std::vector<Topology*>::iterator top = TopList_.begin();
                                          top != TopList_.end(); top++)
      delete *top;
  }
}

// TopologyList::CheckCommand()
/** Check if the command in the arglist pertains to topology files.
  * \return 0 if command was recognized, 1 if not.
  */
int TopologyList::CheckCommand(ArgList &argIn) {
  int pindex;
  // parm <filename> [<tag>]: Add <filename> to parm list
  if (argIn.CommandIs("parm")) {
    std::string parmtag = argIn.getNextTag();
    this->AddParmFile(argIn.getNextString(), parmtag);
    return 0;
  }
  // parmlist: Print list of loaded parm files
  if (argIn.CommandIs("parmlist")) {
    this->Print();
    return 0;
  }
  // parminfo [<parmindex>] [<mask>]: Print information on parm <parmindex> 
  //     (0 by default). If <mask> is given print info on atoms in mask. If
  //     no mask given print overall information.
  if (argIn.CommandIs("parminfo")) {
    pindex = argIn.getNextInteger(0);
    if (pindex>=0 && pindex<(int)TopList_.size()) {
      char *maskarg = argIn.getNextMask();
      if (maskarg!=NULL) 
        TopList_[pindex]->PrintAtomInfo( maskarg );
      else 
        TopList_[pindex]->Summary();
    } else
      mprinterr("Error: parminfo: parm index %i not loaded.\n",pindex);
    return 0;
  }
  // parmwrite out <filename> [<parmindex>]: Write parm <parmindex> to <filename>
  if (argIn.CommandIs("parmwrite")) {
    char *outfilename = argIn.getKeyString("out",NULL);
    if (outfilename==NULL) {
      mprinterr("Error: parmwrite: No output filename specified (use 'out <filename>').\n");
      return 0;
    }
    pindex = argIn.getNextInteger(0);
    if (pindex < 0 || pindex >= (int)TopList_.size()) {
      mprinterr("Error: parmwrite: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    mprintf("\tWriting parm %i (%s) to Amber parm %s\n",pindex,
            TopList_[pindex]->c_str(), outfilename);
    ParmFile pfile;
    pfile.SetDebug( debug_ );
    pfile.Write( *TopList_[pindex], outfilename, ParmFile::AMBERPARM );
    return 0;
  }
  // parmstrip <mask> [<parmindex>]: Strip atoms int mask from parm
  if (argIn.CommandIs("parmstrip")) {
    char *mask0 = argIn.getNextMask();
    pindex = argIn.getNextInteger(0);
    if (pindex < 0 || pindex >= (int)TopList_.size()) {
      mprinterr("Error: parmstrip: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    AtomMask tempMask(mask0);
    // Since want to keep atoms outside mask, invert selection
    tempMask.InvertMask();
    TopList_[pindex]->SetupIntegerMask( tempMask );
    mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(), 
             TopList_[pindex]->Natom() - tempMask.Nselected(), TopList_[pindex]->c_str());
    Topology *tempParm = TopList_[pindex]->modifyStateByMask(tempMask);
    if (tempParm==NULL) 
      mprinterr("Error: parmstrip: Could not strip parm.\n");
    else {
      // Replace parm with stripped version
      tempParm->ParmInfo();
      if (!hasCopies_) delete TopList_[pindex];
      TopList_[pindex] = tempParm;
    }
    return 0;
  }
  // [parm]box [<parmindex>] [x <xval>] [y <yval>] [z <zval>] [alpha <a>] [beta <b>] [gamma <g>]
  //           [nobox]
  // Set the given parm box info to what is specified. If nobox, remove box info.
  if (argIn.CommandIs("box") || argIn.CommandIs("parmbox")) {
    Box pbox;
    pbox.SetX( argIn.getKeyDouble("x",0) );
    pbox.SetY( argIn.getKeyDouble("y",0) );
    pbox.SetZ( argIn.getKeyDouble("z",0) );
    pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
    pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
    pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
    bool nobox = argIn.hasKey("nobox"); 
    // Get parm index
    pindex = argIn.getNextInteger(0);
    if (pindex < 0 || pindex >= (int)TopList_.size()) {
      mprinterr("Error: box: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    if (nobox)
      TopList_[pindex]->ParmBox().SetNoBox();
    else {
      // Fill in missing parm box information from specified parm
      pbox.SetMissingInfo( TopList_[pindex]->ParmBox() );
      TopList_[pindex]->ParmBox() = pbox;
    }
    return 0;
  }
  // solvent [<parmindex>] <mask>
  // Set solvent for the given parm (default 0) based on <mask>
  if (argIn.CommandIs("solvent")) {
    char* maskexpr = argIn.getNextMask();
    if ( maskexpr == NULL ) {
      mprinterr("Error: solvent: No mask specified.\n");
      return 0;
    }
    // Get parm index
    pindex = argIn.getNextInteger(0);
    if (pindex < 0 || pindex >= (int)TopList_.size()) {
      mprinterr("Error: solvent: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    TopList_[pindex]->SetSolvent( maskexpr );
    return 0;
  }
  // parmbondinfo [<parmindex>]: Print bond information for parm <parmindex>
  //     (0 by default).
  if (argIn.CommandIs("parmbondinfo")) {
    pindex = argIn.getNextInteger(0);
    if (pindex>=0 && pindex<(int)TopList_.size()) 
      TopList_[pindex]->PrintBondInfo();
    else
      mprinterr("Error: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmmolinfo [<parmindex>]: Print molecule information for parm
  //     <parmindex> (0 by default).
  if (argIn.CommandIs("parmmolinfo")) {
    pindex = argIn.getNextInteger(0);
    if (pindex>=0 && pindex<(int)TopList_.size())
      TopList_[pindex]->PrintMoleculeInfo();
    else
      mprinterr("Error: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmresinfo [<parmindex>]: Print residue information for parm
  //     <parmindex> (0 by default).
  if (argIn.CommandIs("parmresinfo")) {
    pindex = argIn.getNextInteger(0);
    if (pindex>=0 && pindex<(int)TopList_.size())
      TopList_[pindex]->PrintResidueInfo();
    else
      mprinterr("Error: parm %i not loaded.\n",pindex);
    return 0;
  }
  // bondsearch: Indicate that if bond information not found in topology
  //     it should be determined by distance search.
  if (argIn.CommandIs("bondsearch")) {
    mprintf("\tInfo: Bond info will be determined from distance search if not present.\n");
    bondsearch_=true;
    return 0;
  }
  // molsearch: Indicate that if molecule information not found in 
  //     topology file it should be determined by bonding information.
  if (argIn.CommandIs("molsearch")) {
    mprintf("\tInfo: Molecule info will be determined from bonds if not present.\n");
    molsearch_=true;
    return 0;
  }
  // nobondsearch: Turn off bond search.
  if (argIn.CommandIs("nobondsearch")) {
    mprintf("\tInfo: Bond search is off.\n");
    bondsearch_=false;
    return 0;
  }
  // nomolsearch: Turn off molecule search.
  if (argIn.CommandIs("nomolsearch")) {
    mprintf("\tInfo: Molecule search is off.\n");
    molsearch_=false;
    return 0;
  }
  // Unrecognized parm command
  return 1;
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
    mprinterr("    Error: Could not get parameter file:\n");
    mprinterr("           parmname=%s, pindex=%i\n",parmfilename.c_str(),pindex);
    return NULL;
  }

  return ParmOut;
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list. */
int TopologyList::AddParmFile(std::string const& filename) {
  std::string emptystring;
  return AddParmFile(filename, emptystring);
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list with optional tag. */
int TopologyList::AddParmFile(std::string const& filename, std::string const& ParmTag) 
{
  // Dont let a list that has copies add a new file
  if (hasCopies_) {
    mprintf("    Warning: Attempting to add parm %s to a list that already\n",filename.c_str());
    mprintf("             has copies of parm files. This should not occur.\n");
    mprintf("             Skipping.\n");
    return 0;
  }

  // Check if this file has already been loaded
  if (FindName(filename)!=-1) {
    mprintf("    Warning: Parm %s already loaded, skipping.\n",filename.c_str());
    return 1;
  }

  // If tag specified, check if tag already in use
  if (FindName(ParmTag)!=-1) {
    mprintf("    Warning: Parm tag [%s] already in use.\n",ParmTag.c_str());
    return 1;
  }

  Topology *parm = new Topology();
  parm->SetDebug( debug_ );
  ParmFile pfile;
  pfile.SetDebug( debug_ );
  int err = pfile.Read(*parm, filename.c_str(), bondsearch_, molsearch_);
  if (err!=0) {
    mprinterr("Error: Could not open parm %s\n",filename.c_str());
    delete parm;
    return 1;
  }

  // pindex is used for quick identification of the parm file
  if (debug_>0) 
    mprintf("    PARAMETER FILE %zu: %s\n",TopList_.size(),filename.c_str());
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

// TopologyList::Print()
/// Print list of loaded parameter files
void TopologyList::Print() {
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
