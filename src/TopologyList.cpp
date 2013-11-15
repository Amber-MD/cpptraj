// ParmList
#include "TopologyList.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR 
TopologyList::TopologyList() : debug_(0) {}

// DESTRUCTOR
TopologyList::~TopologyList() {
  Clear();
}

void TopologyList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("TopologyList debug level set to %i\n", debug_);
}

void TopologyList::Clear() {
  for (std::vector<Topology*>::iterator top = TopList_.begin();
                                        top != TopList_.end(); top++)
    delete *top;
  TopList_.clear();
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
  Topology* ParmOut = 0;
  // Get any parm keywords if present
  int pindex = argIn.getKeyInt("parmindex",0);
  std::string parmfilename = argIn.GetStringKey("parm");
  if (!parmfilename.empty()) {
    for (std::vector<Topology*>::const_iterator tf = TopList_.begin();
                                              tf != TopList_.end(); ++tf)
    {
      if ( (*tf)->Tag()==parmfilename || 
           (*tf)->OriginalFilename().MatchFullOrBase( parmfilename ) )
      {
        ParmOut = *tf;
        break;
      }
    }
  } else {
    ParmOut = GetParm(pindex);
  }
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
  // Check if filename/parmtag already in use
  for (std::vector<Topology*>::const_iterator tf = TopList_.begin();
                                              tf != TopList_.end(); ++tf)
  {
    if ( (*tf)->OriginalFilename().Full() == filename ) {
      mprintf("Warning: Parm '%s' already loaded, skipping.\n",filename.c_str());
      return 0;
    }
    if ( !ParmTag.empty() && (*tf)->Tag() == ParmTag ) {
      mprinterr("Error: Parm tag '%s' already in use.\n",ParmTag.c_str());
      return 1;
    }
  }

  Topology *parm = new Topology();
  parm->SetDebug( debug_ );
  parm->SetOffset( offset );
  ParmFile pfile;
  // TODO: Pass in FileName
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
  parm->SetTag( ParmTag );
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

// TopologyList::List()
/** Print list of loaded parameter files */
void TopologyList::List() const {
  if (!TopList_.empty()) {
    mprintf("\nPARAMETER FILES:\n");
    for (std::vector<Topology*>::const_iterator top = TopList_.begin();
                                                top != TopList_.end(); top++)
    {
      mprintf(" %i:", (*top)->Pindex());
      (*top)->Brief();
      if ((*top)->Nframes() > 0)
        mprintf(", %i frames", (*top)->Nframes());
      mprintf("\n");
    }
  }
}
