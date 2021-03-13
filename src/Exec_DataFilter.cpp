#include "Exec_DataFilter.h"
#include "DataFilter.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

void Exec_DataFilter::Help() const {
  mprintf("\t{<dataset arg> min <min> max <max> ...} [out <file>] [name <setname>]\n"
          "\t[{multi | filterset <set> [newset <newname>]}]\n"
          "  Create a data set (optionally named <setname>) containing 1 for\n"
          "  data within given <min> and <max> criteria for each specified\n"
          "  data set. There must be at least one <min> and <max> argument,\n"
          "  and can be as many as there are specified data sets.\n"
          "  If 'multi' is specified then only filter data sets will be created for each\n"
          "  data set instead.\n"
          "  If 'filterset' is specified, the specified <set> will be modified\n"
          "  to only contain '1' frames; cannot be used with 'multi'. If 'newset'\n"
          "  is also specified, a new set will be created containing the '1' frames instead.\n"
          "  The 'filterset' functionality only works for 1D scalar sets.\n");
}

// Exec_DataFilter::Execute()
Exec::RetType Exec_DataFilter::Execute(CpptrajState& State, ArgList& argIn) {
  DataFilter dataFilter;

  if (dataFilter.InitFilter(argIn, State.DSL(), State.DFL(), State.Debug()))
    return CpptrajState::ERR;
/*
  // Get args specific to this exec.
  std::string filterSetName = argIn.GetStringKey("filterset");
  std::string newSetName;
  if (!filterSetName.empty()) {
    // Not intended to work with 'multi'
    if (argIn.hasKey("multi")) {
      mprinterr("Error: 'filterset' can not be used with 'multi' keyword.\n");
      return CpptrajState::ERR;
    }
    newSetName = argIn.GetStringKey("newset");
  }
  // Init and set up Action_FilterByData
  Action_FilterByData filterAction;
  ActionInit state(State.DSL(), State.DFL());
  if (filterAction.Init(argIn, state, State.Debug()) != Action::OK)
    return CpptrajState::ERR;
  size_t nframes = filterAction.DetermineFrames();
  if (nframes < 1) {
    mprinterr("Error: No data to filter. All sets must contain some data.\n");
    return CpptrajState::ERR;
  }
  // Get the set to be filtered if needed.
  DataSet_1D* SetToBeFiltered = 0;
  DataSet_1D* FilteredSet = 0;
  DataSet_double Xvals;
  DataSet_integer* FilterSet = 0;
  if (!filterSetName.empty()) {
    DataSet* ds = filterAction.FilterSet();
    if (ds == 0 || ds->Type() != DataSet::INTEGER) {
      mprinterr("Internal Error: FilterSet was not created.\n");
      return CpptrajState::ERR;
    }
    FilterSet = (DataSet_integer*)ds;
    ds = State.DSL().GetDataSet( filterSetName );
    if (ds == 0) {
      mprinterr("Error: Set to be filtered '%s' not found.\n", filterSetName.c_str());
      return CpptrajState::ERR;
    }
    if (ds->Group() != DataSet::SCALAR_1D)
    {
      mprinterr("Error: '%s' is not 1D scalar.\n", ds->legend());
      return CpptrajState::ERR;
    }
    if (ds->Size() < nframes) {
      mprinterr("Error: Set to be filtered size (%zu) is too small (%zu)\n",
                ds->Size(), nframes);
      return CpptrajState::ERR;
    }
    SetToBeFiltered = (DataSet_1D*)ds;
    // Create new filter set
    
    FilteredSet = (DataSet_1D*)DataSetList::Allocate(SetToBeFiltered->Type());
    MetaData md(SetToBeFiltered->Meta());
    if (!newSetName.empty()) {
      md.SetName(newSetName);
      ds = State.DSL().CheckForSet( md );
      if (ds != 0) {
        mprinterr("Error: New set name '%s' already exists.\n", newSetName.c_str());
        return CpptrajState::ERR;
      }
      mprintf("\tA new filtered set will be created from set '%s'\n", SetToBeFiltered->legend());
    } else
      mprintf("\tFiltering set '%s'\n", SetToBeFiltered->legend());
    FilteredSet->SetMeta( md );
    FilteredSet->SetDim(0, SetToBeFiltered->Dim(0));
    // Do not allocate here to save memory.
  }
 */

  size_t nelements = dataFilter.MinNumElements();
  if (nelements < 1) {
    mprinterr("Error: No data to filter. All sets must contain some data.\n");
    return CpptrajState::ERR;
  }
 
  ProgressBar progress( nelements );
  //int newidx = 0;
  for (size_t idx = 0; idx != nelements; idx++) {
    progress.Update( idx );
    dataFilter.FilterIndex( idx );
/*    // Filter does not need frame but does need trajout num.
    ActionFrame frm(0, frame);
    filterAction.DoAction(frame, frm);
    if (SetToBeFiltered != 0) {
      if ((*FilterSet)[frame] == 1) {
        Xvals.AddElement(SetToBeFiltered->Xcrd(frame));
        FilteredSet->Add(newidx++, SetToBeFiltered->VoidPtr(frame));
      }
    }*/
  }
/*
  // Add/replace the filtered set if necessary
  if (SetToBeFiltered != 0) {
    if (newSetName.empty())
      State.DSL().RemoveSet( SetToBeFiltered );
    DataSetList::DataListType inputSets(1);
    inputSets[0] = FilteredSet;
    State.DSL().AddOrAppendSets( FilteredSet->Dim(0).Label(), Xvals.Data(), inputSets );
  } */
  if (dataFilter.Finalize())
    return CpptrajState::ERR;
  // Trigger master datafile write just in case
  State.MasterDataFileWrite();
  return CpptrajState::OK;
}
