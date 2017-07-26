#include <algorithm> // std::min, std::max
#include "Exec_CombineCoords.h"
#include "CpptrajStdio.h"

void Exec_CombineCoords::Help() const {
  mprintf("\t<crd1> <crd2> ... [parmname <topname>] [crdname <crdname>] [nobox]\n"
          "  Combine two or more COORDS data sets.\n");
}

// Exec_CombineCoords::Execute()
Exec::RetType Exec_CombineCoords::Execute(CpptrajState& State, ArgList& argIn) {
# ifdef MPI
  if (Parallel::World().Size() > 1) {
    mprinterr("Error: '%s' cannot be run in parallel.\n", argIn.Command());
    return CpptrajState::ERR;
  }
# endif
  std::string parmname = argIn.GetStringKey("parmname");
  std::string crdname  = argIn.GetStringKey("crdname");
  bool nobox = argIn.hasKey("nobox");
  // Get COORDS DataSets.
  std::vector<DataSet_Coords*> CRD;
  std::string setname = argIn.GetStringNext();
  while (!setname.empty()) {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL().FindCoordsSet( setname );
    if (ds == 0) {
      mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
      return CpptrajState::ERR;
    }
    CRD.push_back( ds );
    setname = argIn.GetStringNext();
  }
  if (CRD.size() < 2) {
    mprinterr("Error: %s: Must specify at least 2 COORDS data sets\n", argIn.Command());
    return CpptrajState::ERR;
  }
  // Only add the topology to the list if parmname specified
  bool addTop = true;
  Topology CombinedTop;
  CombinedTop.SetDebug( State.Debug() );
  if (parmname.empty()) {
    parmname = CRD[0]->Top().ParmName() + "_" + CRD[1]->Top().ParmName();
    addTop = false;
  }
  CombinedTop.SetParmName( parmname, FileName() );
  // TODO: Check Parm box info.
  enum BoxStatus {NOT_SET=0, SET, INVALID};
  BoxStatus boxStatus = NOT_SET;
  if (nobox) boxStatus = INVALID;
  Box combinedBox;
  size_t minSize = CRD[0]->Size();
  for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum) {
    if (CRD[setnum]->Size() < minSize)
      minSize = CRD[setnum]->Size();
    if (CRD[setnum]->CoordsInfo().HasBox()) {
      if (boxStatus == NOT_SET) {
        combinedBox = CRD[setnum]->CoordsInfo().TrajBox();
        boxStatus = SET;
      } else if (boxStatus == SET) {
        // Make sure it is the same type of box. TODO Check angles
        if (combinedBox.Type() != CRD[setnum]->CoordsInfo().TrajBox().Type())
        {
          mprintf("Warning: COORDS '%s' box type '%s' differs from other COORDS. Disabling box.\n",
                  CRD[setnum]->legend(), CRD[setnum]->CoordsInfo().TrajBox().TypeName());
          combinedBox.SetNoBox();
          boxStatus = INVALID;
        }
      }
    }
    CombinedTop.AppendTop( CRD[setnum]->Top() );
  }
  CombinedTop.SetParmBox( combinedBox );
  CombinedTop.Brief("Combined parm:");
  if (addTop) {
    if (State.AddTopology( CombinedTop, parmname )) return CpptrajState::ERR;
  }

  // Combine coordinates from all sets for all frames.
  if (crdname.empty())
    crdname = CRD[0]->Meta().Legend() + "_" + CRD[1]->Meta().Legend();
  mprintf("\tCombining %zu frames from each set into %s\n", minSize, crdname.c_str());
  DataSet_Coords* CombinedCrd = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS,crdname,"CRD");
  if (CombinedCrd == 0) {
    mprinterr("Error: Could not create COORDS data set.\n");
    return CpptrajState::ERR;
  }
  // FIXME: Only copying coords+box for now
  CoordinateInfo combinedInfo(combinedBox, false, false, false);
  CombinedCrd->CoordsSetup( CombinedTop, combinedInfo );
  Frame CombinedFrame = CombinedCrd->AllocateFrame();
  // Allocate space for input frames from each COORDS set.
  std::vector<Frame> InputFrames;
  for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    InputFrames.push_back( CRD[setnum]->AllocateFrame() );
  // Loop over all frames from each COORDS set.
  for (size_t nf = 0; nf != minSize; ++nf) {
    double *output = CombinedFrame.xAddress();
    // Coords
    for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    {
      Frame& input = InputFrames[setnum];
      CRD[setnum]->GetFrame( nf, input );
      std::copy(input.xAddress(), input.xAddress()+input.size(), output);
      output += input.size();
    }
    // Box info
    if (combinedBox.Type() != Box::NOBOX) {
      double* cBox = CombinedFrame.bAddress();
      // Only use angles from first coords set.
      std::copy(InputFrames[0].bAddress(), InputFrames[0].bAddress()+6, cBox);
      for (unsigned int setnum = 1; setnum < CRD.size(); ++setnum) {
        // Use max X/Y/Z among coords
        cBox[0] = std::max(cBox[0], InputFrames[setnum].BoxCrd().BoxX());
        cBox[1] = std::max(cBox[1], InputFrames[setnum].BoxCrd().BoxY());
        cBox[2] = std::max(cBox[2], InputFrames[setnum].BoxCrd().BoxZ());
      }
    }
    CombinedCrd->AddFrame( CombinedFrame );
  }
/* FIXME: This code is fast but only works for DataSet_Coords_CRD
  Frame::CRDtype CombinedFrame( CombinedTop->Natom() * 3 );
  for (size_t nf = 0; nf != minSize; ++nf) {
    size_t offset = 0;
    for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    {
      size_t crd_offset = (size_t)CRD[setnum]->Top().Natom() * 3;
      std::copy( CRD[setnum]->CRD(nf).begin(), CRD[setnum]->CRD(nf).begin() + crd_offset,
                 CombinedFrame.begin() + offset );
      offset += crd_offset;
    }
    CombinedCrd->AddCRD( CombinedFrame );
  }
*/
  return CpptrajState::OK;
}
