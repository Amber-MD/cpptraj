// ParmList
#include "TopologyList.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR 
TopologyList::TopologyList() {}

// DESTRUCTOR
TopologyList::~TopologyList() {
  Clear();
}

void TopologyList::Clear() {
  for (std::vector<Topology*>::iterator top = TopList_.begin();
                                        top != TopList_.end(); top++)
    delete *top;
  TopList_.clear();
  FileList::Clear();
}

// TopologyList::GetParm()
/** Return the parm structure with index num. */
Topology* TopologyList::GetParm(int num) const {
  if (num>=(int)TopList_.size() || num<0) return 0;
  return TopList_[num];
}

// TopologyList::GetParmByIndex()
Topology* TopologyList::GetParmByIndex(ArgList& argIn) const {
  int pindex = argIn.getNextInteger(0);
  Topology* parm = GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parm index %i not loaded.\n",pindex);
    return 0;
  }
  return parm;
}

// TopologyList::GetParm()
/** Return the parm structure based on arguments in the given arg list. 
  *   parm <parm name>
  *   parmindex <parm index>
  * \param argIn argument list that contains parm-related keyword
  * \return parm specified by 'parm' or 'parmindex', or the first parm. null on error.
  */
Topology* TopologyList::GetParm(ArgList &argIn) const {
  // Get any parm keywords if present
  std::string parmfilename = argIn.GetStringKey("parm");
  int pindex = argIn.getKeyInt("parmindex",0);
  // Associate trajectory with parameter file. Associate with default parm if none specified
  if (!parmfilename.empty())
    pindex = FindName(parmfilename);
  Topology *ParmOut = GetParm(pindex);
  if (ParmOut==0) {
    mprintf("Warning: Could not get parameter file:\n");
    mprintf("Warning: parmname=%s, pindex=%i\n",parmfilename.c_str(),pindex);
    return 0;
  }

  return ParmOut;
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list. */
int TopologyList::AddParmFile(std::string const& filename) {
  return AddParmFile(filename, std::string(), true, -1.0);
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list with optional tag. */
int TopologyList::AddParmFile(std::string const& filename, std::string const& ParmTag,
                              bool bondsearch, double offset) 
{
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
  AddNameWithTag( pfile.ParmFilename(), ParmTag );
  return 0;
}

// TopologyList::WriteParm()
int TopologyList::WriteParm(ArgList& argIn) const {
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: %s: No output filename specified (use 'out <filename>').\n", argIn.Command());
    return 1;
  }
  Topology* parm = GetParmByIndex( argIn );
  if (parm == 0) return 1;
  mprintf("\tWriting parm %i (%s) to Amber parm %s\n",parm->Pindex(),
          parm->c_str(), outfilename.c_str());
  ParmFile pfile;
  pfile.Write( *parm, outfilename, ParmFile::AMBERPARM, debug_ );
  return 0;
}

// TopologyList::ReplaceParm()
// TODO: Use mem address instead of pindex?
void TopologyList::ReplaceParm(int pindex, Topology* newParm) {
  if (pindex < 0 || pindex >= (int)TopList_.size()) {
    mprinterr("Error: ReplaceParm: parm index %i out of bounds.\n",pindex);
    return;
  }
  delete TopList_[pindex];
  TopList_[pindex] = newParm;
}

// TopologyList::List()
/** Print list of loaded parameter files */
void TopologyList::List() const {
  mprintf("\nPARAMETER FILES:\n");
  if (TopList_.empty()) {
    mprintf("  No parameter files defined.\n");
    return;
  }
  for (std::vector<Topology*>::const_iterator top = TopList_.begin(); top != TopList_.end(); top++)
  {
    mprintf(" %i:", (*top)->Pindex());
    (*top)->ParmInfo();
    if ((*top)->Nframes() > 0)
      mprintf(", %i frames", (*top)->Nframes());
    mprintf("\n");
  }
}
