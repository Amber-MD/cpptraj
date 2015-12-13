#include "Exec_DataSetCmd.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

void Exec_DataSetCmd::Help() const {
  mprintf("\t{ legend <legend> <set> | makexy <Xset> <Yset> [name <name>] |\n"
          "\t  cat <set0> <set1> ... [name <name>] [nooffset]\n"
          "\t  [mode <mode>] [type <type>] <set arg1> [<set arg 2> ...] }\n"
          "\t<mode>: ");
  for (int i = 0; i != (int)MetaData::M_MATRIX; i++) // TODO: Allow matrix?
    mprintf(" %s", MetaData::ModeString((MetaData::scalarMode)i));
  mprintf("\n\t<type>: ");
  for (int i = 0; i != (int)MetaData::DIST; i++)
    mprintf(" %s", MetaData::TypeString((MetaData::scalarType)i));
  mprintf("\n\tOptions for 'type noe':\n"
          "\t  %s\n"
          "  legend: Set the legend for a single data set\n"
          "  makexy: Create new data set with X values from one set and Y values from another.\n"
          "  cat   : Concatenate 2 or more data sets.\n"
          "  Otherwise, change the mode/type for one or more data sets.\n",
          AssociatedData_NOE::HelpText);
}

Exec::RetType Exec_DataSetCmd::Execute(CpptrajState& State, ArgList& argIn) {
  if (argIn.Contains("legend")) { // Set legend for one data set
    std::string legend = argIn.GetStringKey("legend");
    DataSet* ds = State.DSL().GetDataSet( argIn.GetStringNext() );
    if (ds == 0) return CpptrajState::ERR;
    mprintf("\tChanging legend '%s' to '%s'\n", ds->legend(), legend.c_str());
    ds->SetLegend( legend );
  // ---------------------------------------------
  } else if (argIn.hasKey("makexy")) { // Combine values from two sets into 1
    std::string name = argIn.GetStringKey("name");
    DataSet* ds1 = State.DSL().GetDataSet( argIn.GetStringNext() );
    DataSet* ds2 = State.DSL().GetDataSet( argIn.GetStringNext() );
    if (ds1 == 0 || ds2 == 0) return CpptrajState::ERR;
    if (ds1->Ndim() != 1 || ds2->Ndim() != 1) {
      mprinterr("Error: makexy only works for 1D data sets.\n");
      return CpptrajState::ERR;
    }
    DataSet* ds3 = State.DSL().AddSet( DataSet::XYMESH, name, "XY" );
    if (ds3 == 0) return CpptrajState::ERR;
    mprintf("\tUsing values from '%s' as X, values from '%s' as Y, output set '%s'\n",
            ds1->legend(), ds2->legend(), ds3->legend());
    DataSet_1D const& ds_x = static_cast<DataSet_1D const&>( *ds1 );
    DataSet_1D const& ds_y = static_cast<DataSet_1D const&>( *ds2 );
    DataSet_1D&       out  = static_cast<DataSet_1D&>( *ds3 );
    size_t nframes = std::min( ds_x.Size(), ds_y.Size() );
    if (ds_x.Size() != ds_y.Size())
      mprintf("Warning: Data sets do not have equal sizes, only using %zu frames.\n", nframes);
    double XY[2];
    for (size_t i = 0; i != nframes; i++) {
      XY[0] = ds_x.Dval(i);
      XY[1] = ds_y.Dval(i);
      out.Add( i, XY );
    }
  // ---------------------------------------------
  } else if (argIn.hasKey("cat")) { // Concatenate two or more data sets
    std::string name = argIn.GetStringKey("name");
    bool use_offset = !argIn.hasKey("nooffset");
    DataSet* ds3 = State.DSL().AddSet( DataSet::XYMESH, name, "CAT" );
    if (ds3 == 0) return CpptrajState::ERR;
    DataSet_1D& out = static_cast<DataSet_1D&>( *ds3 );
    mprintf("\tConcatenating sets into '%s'\n", out.legend());
    if (use_offset)
      mprintf("\tX values will be offset.\n");
    else
      mprintf("\tX values will not be offset.\n");
    std::string dsarg = argIn.GetStringNext();
    double offset = 0.0;
    while (!dsarg.empty()) {
      DataSetList dsl = State.DSL().GetMultipleSets( dsarg );
      double XY[2];
      for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
      {
        if ( (*ds)->Type() != DataSet::INTEGER &&
             (*ds)->Type() != DataSet::DOUBLE &&
             (*ds)->Type() != DataSet::FLOAT &&
             (*ds)->Type() != DataSet::XYMESH )
        {
          mprintf("Warning: '%s': Concatenation only supported for 1D scalar data sets.\n",
                  (*ds)->legend());
        } else {
          DataSet_1D const& set = static_cast<DataSet_1D const&>( *(*ds) );
          mprintf("\t\t'%s'\n", set.legend());
          for (size_t i = 0; i != set.Size(); i++) {
            XY[0] = set.Xcrd( i ) + offset;
            XY[1] = set.Dval( i );
            out.Add( i, XY ); // NOTE: value of i does not matter for mesh
          }
          if (use_offset) offset = XY[0];
        }
      }
      dsarg = argIn.GetStringNext();
    }
  // ---------------------------------------------
  } else { // Default: change mode/type for one or more sets.
    std::string modeKey = argIn.GetStringKey("mode");
    std::string typeKey = argIn.GetStringKey("type");
    if (modeKey.empty() && typeKey.empty()) {
      mprinterr("Error: No valid keywords specified.\n");
      return CpptrajState::ERR;
    }
    // First determine mode if specified.
    MetaData::scalarMode dmode = MetaData::UNKNOWN_MODE;
    if (!modeKey.empty()) {
      dmode = MetaData::ModeFromKeyword( modeKey );
      if (dmode == MetaData::UNKNOWN_MODE) {
        mprinterr("Error: Invalid mode keyword '%s'\n", modeKey.c_str());
        return CpptrajState::ERR;
      }
    }
    // Next determine type if specified.
    MetaData::scalarType dtype = MetaData::UNDEFINED;
    if (!typeKey.empty()) {
      dtype = MetaData::TypeFromKeyword( typeKey, dmode );
      if (dtype == MetaData::UNDEFINED) {
        mprinterr("Error: Invalid type keyword '%s'\n", typeKey.c_str());
        return CpptrajState::ERR;
      }
    }
    // Additional options for type 'noe'
    AssociatedData_NOE noeData;
    if (dtype == MetaData::NOE) {
      if (noeData.NOE_Args(argIn))
        return CpptrajState::ERR;
    }
    if (dmode != MetaData::UNKNOWN_MODE)
      mprintf("\tDataSet mode = %s\n", MetaData::ModeString(dmode));
    if (dtype != MetaData::UNDEFINED)
      mprintf("\tDataSet type = %s\n", MetaData::TypeString(dtype));
    // Loop over all DataSet arguments 
    std::string ds_arg = argIn.GetStringNext();
    while (!ds_arg.empty()) {
      DataSetList dsl = State.DSL().GetMultipleSets( ds_arg );
      for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
      {
        if ( (*ds)->Ndim() != 1 ) // TODO remove restriction
          mprintf("Warning:\t\t'%s': Can only set mode/type for 1D data sets.\n",
                  (*ds)->legend());
        else {
          if ( dtype == MetaData::NOE ) (*ds)->AssociateData( &noeData );
          mprintf("\t\t'%s'\n", (*ds)->legend());
          MetaData md = (*ds)->Meta();
          md.SetScalarMode( dmode );
          md.SetScalarType( dtype );
          (*ds)->SetMeta( md );
        }
      }
      ds_arg = argIn.GetStringNext();
    }
  }
  return CpptrajState::OK;
}
