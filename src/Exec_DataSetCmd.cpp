#include <algorithm> // std::min, std::max
#include "Exec_DataSetCmd.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_Vector.h"
#include "DataSet_string.h"
#include "DataSet_Mesh.h"
#include "DataSet_Mat3x3.h"
#include "StringRoutines.h"

// Exec_DataSetCmd::Help()
void Exec_DataSetCmd::Help() const {
  mprintf("\t{legend|makexy|vectorcoord|cat|make2d|droppoints|keeppoints|remove|\n"
          "\t dim|outformat|invert|shift|mode|type} <options>\n");
  mprintf("  Type 'help dataset <cmd>' for detailed subcommand help.\n");
}

// Exec_DataSetCmd::Help()
void Exec_DataSetCmd::Help(ArgList& argIn) const {
  if (argIn.hasKey("legend")) {
    mprintf("  legend <legend> <set>\n"
            "    Set the legend for a single data set.\n");
  } else if (argIn.hasKey("makexy")) {
    mprintf("  makexy <Xset> <Yset> [name <name>]\n"
            "    Create new data set with X values from one set and Y values from another.\n");
  } else if (argIn.hasKey("vectorcoord")) {
    mprintf("  vectorcoord {X|Y|Z} [name <name>]\n"
            "    Extract X, Y, or Z component of vector data into new set.\n");
  } else if (argIn.hasKey("cat")) {
    mprintf("  cat <set0> <set1> ... [name <name>] [nooffset]\n"
            "    Concatenate 2 or more data sets.\n");
  } else if (argIn.hasKey("make2d")) {
    mprintf("  make2d <1D set> cols <ncols> rows <nrows> [name <name>]\n"
            "    Create new 2D data set from 1D data set, assumes row-major ordering.\n");
  } else if (argIn.hasKey("droppoints") || argIn.hasKey("keeppoints")) {
    Help_ModifyPoints();
  } else if (argIn.hasKey("remove")) {
    mprintf("  remove <criterion> <select> <value> [and <value2>] [<set selection>]\n"
            "      <criterion>: ");
    for (int i = 1; i < (int)N_C; i++)
      mprintf(" '%s'", CriterionKeys[i]);
    mprintf("\n      <select>   : ");
    for (SelectPairType const* ptr = SelectKeys; ptr->key_ != 0; ptr++)
      mprintf(" '%s'", ptr->key_);
    mprintf("\n    Remove data sets according to specified criterion and selection.\n");
  } else if (argIn.hasKey("dim")) {
    Help_ChangeDim();
  } else if (argIn.hasKey("outformat")) {
    mprintf("  outformat {double|scientific|general} <set arg1> [<set arg 2> ...]\n"
            "    Change output format of double-precision data:\n"
            "      double     - \"Normal\" output, e.g. 0.4032\n"
            "      scientific - Scientific \"E\" notation output, e.g. 4.032E-1\n"
            "      general    - Use 'double' or 'scientific', whichever is shortest.\n");
  } else if (argIn.hasKey("invert")) {
    Help_InvertSets();
  } else if (argIn.hasKey("shift")) {
    Help_Shift();
  } else if (argIn.hasKey("mode")) {
    mprintf("  [mode <mode>] [type <type>] <set arg1> [<set arg 2> ...]\n");
    mprintf("      <mode>: ");
    for (int i = 0; i != (int)MetaData::UNKNOWN_MODE; i++)
      mprintf(" '%s'", MetaData::ModeString((MetaData::scalarMode)i));
    mprintf("    Change the mode for one or more data sets.\n");
  } else if (argIn.hasKey("type")) {
    mprintf("\n      <type>: ");
    for (int i = 0; i != (int)MetaData::UNDEFINED; i++)
      mprintf(" '%s'", MetaData::TypeString((MetaData::scalarType)i));
    mprintf("\n    Options for 'type noe':\n"
            "      %s\n", AssociatedData_NOE::HelpText);
    mprintf("    Change the type for one or more data sets.\n");
  } else
    Help();
}

// Exec_DataSetCmd::Execute()
Exec::RetType Exec_DataSetCmd::Execute(CpptrajState& State, ArgList& argIn) {
  modifiedSets_.clear();
  RetType err = CpptrajState::OK;
  if (argIn.Contains("legend")) {         // Set legend for one data set
    std::string legend = argIn.GetStringKey("legend");
    DataSet* ds = State.DSL().GetDataSet( argIn.GetStringNext() );
    if (ds == 0) return CpptrajState::ERR;
    mprintf("\tChanging legend '%s' to '%s'\n", ds->legend(), legend.c_str());
    ds->SetLegend( legend );
  // ---------------------------------------------
  } else if (argIn.hasKey("outformat")) {   // Change double precision set output format
    err = ChangeOutputFormat(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("remove")) {      // Remove data sets by various criteria
    err = Remove(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("makexy")) {      // Combine values from two sets into 1
    err = MakeXY(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("make2d")) {      // Create 2D matrix from 1D set
    err = Make2D(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("vectorcoord")) { // Extract vector X/Y/Z coord as new set
    err = VectorCoord(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("filter")) {      // Filter points in data set to make new data set
    err = Filter(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("cat")) {         // Concatenate two or more data sets
    err = Concatenate(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("droppoints")) {  // Drop points from set
    err = ModifyPoints(State, argIn, true);
  // ---------------------------------------------
  } else if (argIn.hasKey("keeppoints")) {  // Keep points in set
    err = ModifyPoints(State, argIn, false);
  // ---------------------------------------------
  } else if (argIn.hasKey("dim")) {         // Modify dimension of set(s)
    err = ChangeDim(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("invert")) {      // Invert set(s) X/Y, create new sets
    err = InvertSets(State, argIn);
  // ---------------------------------------------
  } else if (argIn.hasKey("shift")) {       // Shift data in set(s) that match criteria by offset
    err = ShiftData(State, argIn);
  // ---------------------------------------------
  } else {                                  // Default: change mode/type for one or more sets.
    err = ChangeModeType(State, argIn);
  }
  // If any sets were modified, write out data again.
  if (err == CpptrajState::OK && !modifiedSets_.empty()) {
    State.DFL().ResetWriteStatIfContain( modifiedSets_ );
    State.DFL().WriteAllDF();
  }
  return err;
}

const char* Exec_DataSetCmd::CriterionKeys[] = { 0, "ifaverage", "ifsize", "ifmode", "iftype" };

Exec_DataSetCmd::SelectPairType Exec_DataSetCmd::SelectKeys[] = {
  {EQUAL,        "equal"},
  {EQUAL,        "=="},
  {NOT_EQUAL,    "notequal"},
  {NOT_EQUAL,    "!="},
  {LESS_THAN,    "lessthan"},
  {LESS_THAN,    "<"},
  {GREATER_THAN, "greaterthan"},
  {GREATER_THAN, ">"},
  {BETWEEN,      "between"},
  {OUTSIDE,      "outside"},
  {UNKNOWN_S,    0}
};

void Exec_DataSetCmd::Help_ModifyPoints() {
  mprintf("  drop|keep}points {range <range arg> | [start <#>] [stop <#>] [offset <#>]}\n"
          "                   [name <output set>] <set arg1> ...\n"
          "    Drop specified points from or keep specified points in data set(s).\n"
          "    This command currently only sowrks for 1D scalar and 3x3 matrix sets.\n");
}

/** Add the X and Y values from set in at idx to set out. */
static void KeepPoint_1D(const DataSet* in, DataSet* out, int idx) {
  DataSet_1D const& set1d = static_cast<DataSet_1D const&>( *in );
  DataSet_Mesh& mesh = static_cast<DataSet_Mesh&>( *out );
  mesh.AddXY( set1d.Xcrd(idx), set1d.Dval(idx) );
}

/** Add matrix from set in at idx to set out. */
static void KeepPoint_Mat3x3(const DataSet* in, DataSet* out, int idx) {
  DataSet_Mat3x3 const& setmat = static_cast<DataSet_Mat3x3 const&>( *in );
  DataSet_Mat3x3& mat = static_cast<DataSet_Mat3x3&>( *out );
  mat.AddMat3x3( setmat[idx] );
}

// Exec_DataSetCmd::ModifyPoints()
Exec::RetType Exec_DataSetCmd::ModifyPoints(CpptrajState& State, ArgList& argIn, bool drop) {
  const char* mode;
  if (drop)
    mode = "Drop";
  else
    mode = "Kee";
  // Hold pointer to function used to Keep points for data set type.
  typedef void (*KeepFxnType)(const DataSet*, DataSet*, int);
  // Keywords
  std::string name = argIn.GetStringKey("name");
  int start = argIn.getKeyInt("start", 0) - 1;
  int stop = argIn.getKeyInt("stop", -1);
  int offset = argIn.getKeyInt("offset", -1);
  Range points;
  if (start < 0 && stop < 0 && offset < 0) {
    std::string rangearg = argIn.GetStringKey("range");
    if (rangearg.empty()) {
      mprinterr("Error: Must specify range or start/stop/offset.\n");
      return CpptrajState::ERR;
    }
    points.SetRange( rangearg );
    if (points.Empty()) {
      mprinterr("Error: Range '%s' is empty.\n", rangearg.c_str());
      return CpptrajState::ERR;
    }
    mprintf("\t%sping points in range %s\n", mode, rangearg.c_str());
    // User args start from 1
    points.ShiftBy(-1);
  }
  // Get data set to drop/keep points from
  // Loop over all DataSet arguments 
  std::string ds_arg = argIn.GetStringNext();
  while (!ds_arg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( ds_arg );
    for (DataSetList::const_iterator it = dsl.begin(); it != dsl.end(); ++it)
    {
      DataSet* DS = *it;
      if (DS->Size() < 1) {
        mprinterr("Error: Set '%s' is empty.\n", DS->legend());
        return CpptrajState::ERR;
      }
      // Restrict to 1D sets for now TODO more types
      KeepFxnType keepFxn;
      DataSet::DataType outType;
      if (DS->Group() == DataSet::SCALAR_1D) {
        keepFxn = KeepPoint_1D;
        outType = DataSet::XYMESH;
      } else if (DS->Type() == DataSet::MAT3X3) {
        keepFxn = KeepPoint_Mat3x3;
        outType = DataSet::MAT3X3;
      } else {
        mprinterr("Error: Currently only works for 1D scalar or 3x3 matrix data sets.\n");
        return CpptrajState::ERR;
      }
      // Output data set.
      // NOTE: We want to preserve the original X values. The easiest way
      //       to do this is to always make the output set an XY mesh
      //       regardless of the input set type. Ideally we would detect
      //       if a set is monotonic in X and allocate appropriately, but
      //       since currently there are no non-monotonic (i.e. sparse)
      //       integer/float DataSet_1D classes, XYMESH is the easiest
      //       solution.
      DataSet* out = 0;
      if (name.empty()) {
        // Modifying this set. Create new temporary set.
        out = State.DSL().Allocate( outType );
        if (out == 0) return CpptrajState::ERR;
        // This gives the new set same MetaData, format etc as the original.
        *out = *DS;
        mprintf("\tOverwriting set '%s'\n", DS->legend());
      } else {
        // Write to new set
        MetaData md = DS->Meta();
        md.SetName( name );
        out = State.DSL().AddSet(outType, md);
        if (out == 0) return CpptrajState::ERR;
        mprintf("\tNew set is '%s'\n", out->Meta().PrintName().c_str());
      }
      // NOTE: Allocate() must be valid for the set type!
      out->Allocate(DataSet::SizeArray(1, DS->Size()));
      if (points.Empty()) {
        // Drop by start/stop/offset. Set defaults if needed
        if (start < 0)  start = 0;
        if (stop < 0)   stop = DS->Size();
        if (offset < 0) offset = 1;
        mprintf("\t%sping points from %i to %i, step %i\n", mode, start+1, stop, offset);
        for (int idx = start; idx < stop; idx += offset)
          points.AddToRange( idx );
      } // TODO check that range values are valid?
      if (State.Debug() > 0) mprintf("DEBUG: Keeping points:");
      Range::const_iterator pt = points.begin();
      int idx = 0;
      if (drop) {
        // Drop points TODO should idx be unsigned?
        for (; idx < (int)DS->Size(); idx++) {
          if (pt == points.end()) break;
          if (*pt != idx) {
            if (State.Debug() > 0) mprintf(" %i", idx + 1);
            keepFxn(DS, out, idx);
          } else
            ++pt;
        }
        // Keep all remaining points
        for (; idx < (int)DS->Size(); idx++) {
          if (State.Debug() > 0) mprintf(" %i", idx + 1);
          keepFxn(DS, out, idx);
        }
      } else {
        // Keep points
        for (; pt != points.end(); pt++) {
          if (*pt >= (int)DS->Size()) break;
          if (State.Debug() > 0) mprintf(" %i", *pt + 1);
          keepFxn(DS, out, *pt);
        }
      }
      if (State.Debug() > 0) mprintf("\n");
      if (name.empty()) {
        // Replace old set with new set
        State.DSL().RemoveSet( DS );
        State.DSL().AddSet( out );
      }
    } // END loop over sets
    ds_arg = argIn.GetStringNext();
  } // END loop over set args
  return CpptrajState::OK;
}

// Exec_DataSetCmd::VectorCoord()
Exec::RetType Exec_DataSetCmd::VectorCoord(CpptrajState& State, ArgList& argIn) {
  // Keywords
  std::string name = argIn.GetStringKey("name");
  int tgtIdx;
  if (argIn.hasKey("X"))
    tgtIdx = 0;
  else if (argIn.hasKey("Y"))
    tgtIdx = 1;
  else if (argIn.hasKey("Z"))
    tgtIdx = 2;
  else {
    mprinterr("Error: 'vectorcoord' requires specifying X, Y, or Z.\n");
    return CpptrajState::ERR;
  }
  // Data set(s)
  typedef std::vector<DataSet_Vector*> DVarray;
  DVarray inputSets;
  std::string vecsetarg = argIn.GetStringNext();
  if (vecsetarg.empty())
    mprintf("Warning: No data set arguments specified.\n");
  while (!vecsetarg.empty()) {
    DataSetList dsl1 = State.DSL().GetMultipleSets( vecsetarg );
    for (DataSetList::const_iterator it = dsl1.begin(); it != dsl1.end(); ++it)
    {
      if ( (*it)->Group() != DataSet::VECTOR_1D) {
        mprintf("Warning: '%s' 'vectorcoord' only works with vector data sets.\n", (*it)->legend());
      } else if ( (*it)->Size() < 1) {
        mprintf("Warning: '%s' is empty.\n", (*it)->legend());
      } else {
        inputSets.push_back( static_cast<DataSet_Vector*>( *it ) );
      }
    }
    vecsetarg = argIn.GetStringNext();
  }
  if (inputSets.empty()) {
    mprinterr("Error: 'vectorcoord': No data sets selected.\n");
    return CpptrajState::ERR;
  }
  mprintf("\t%zu sets.\n", inputSets.size());
  // Default name
  if (name.empty())
    name = State.DSL().GenerateDefaultName( "COORD" );
  // Loop over input sets
  int idx = -1;
  if (inputSets.size() > 1)
    idx = 0;
  for (DVarray::const_iterator it = inputSets.begin(); it != inputSets.end(); ++it)
  {
    // Create output set.
    MetaData md( name );
    if (idx > -1) 
      md.SetIdx( idx++ );
    static const char* XYZchar[3] = { "X", "Y", "Z" };
    DataSet* out = State.DSL().AddSet( DataSet::DOUBLE, md );
    if (out == 0) return CpptrajState::ERR;
    // Extract data
    mprintf("\tExtracting %s coordinate from vector %s to %s\n",
            XYZchar[tgtIdx], (*it)->legend(), out->Meta().PrintName().c_str());
    DataSet_Vector const& vec = static_cast<DataSet_Vector const&>( *(*it) );
    for (unsigned int n = 0; n != vec.Size(); n++) {
      double d = vec.VXYZ(n)[tgtIdx];
      out->Add( n, &d );
    }
  }
  return CpptrajState::OK;
}

// Exec_DataSetCmd::ChangeOutputFormat()
Exec::RetType Exec_DataSetCmd::ChangeOutputFormat(CpptrajState const& State, ArgList& argIn)
{
  TextFormat::FmtType fmt;
  if (argIn.hasKey("double"))
    fmt = TextFormat::DOUBLE;
  else if (argIn.hasKey("scientific"))
    fmt = TextFormat::SCIENTIFIC;
  else if (argIn.hasKey("general"))
    fmt = TextFormat::GDOUBLE;
  else {
    mprinterr("Error: Expected either 'double', 'scientific', or 'general'\n");
    return CpptrajState::ERR;
  }
  // Loop over all DataSet arguments 
  std::string ds_arg = argIn.GetStringNext();
  while (!ds_arg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( ds_arg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
      if ((*ds)->SetupFormat().SetFormatType(fmt))
        mprintf("\tSet '%s' output format changed to '%s'\n",
                (*ds)->legend(), TextFormat::typeDescription(fmt));
    ds_arg = argIn.GetStringNext();
  }
  return CpptrajState::OK;
}

// Exec_DataSetCmd::Remove()
Exec::RetType Exec_DataSetCmd::Remove(CpptrajState& State, ArgList& argIn) {
  std::string status;
  // Get criterion type
  CriterionType criterion = UNKNOWN_C;
  for (int i = 1; i < (int)N_C; i++)
    if (argIn.hasKey( CriterionKeys[i] )) {
      criterion = (CriterionType)i;
      status.assign( CriterionKeys[i] );
      break;
    }
  if (criterion == UNKNOWN_C) {
    mprinterr("Error: No criterion specified for 'remove'.\n");
    return CpptrajState::ERR;
  }
  // Get select type
  SelectType select = UNKNOWN_S;
  std::string val1, val2;
  for (const SelectPairType* ptr = SelectKeys; ptr->key_ != 0; ptr++)
    if (argIn.Contains( ptr->key_ )) {
      select = ptr->type_;
      val1 = argIn.GetStringKey( ptr->key_ );
      status.append( " " + std::string(ptr->key_) + " " + val1 );
      // Get 'and' value for between/outside. TODO put nargs in SelectPairType?
      if (select == BETWEEN || select == OUTSIDE) {
        val2 = argIn.GetStringKey("and");
        if (val2.empty()) {
          mprinterr("Error: Missing 'and' value for selection '%s'\n", ptr->key_);
          return CpptrajState::ERR;
        }
        status.append(" and " + val2);
      }
      break;
    }
  if (select == UNKNOWN_S || val1.empty()) {
    mprinterr("Error: No selection specified for 'remove'.\n");
    return CpptrajState::ERR;
  }
  if ( (criterion == SMODE || criterion == STYPE) &&
       (select != EQUAL && select != NOT_EQUAL) )
  {
    mprinterr("Error: Specified select not valid for criterion '%s'\n", CriterionKeys[criterion]);
    return CpptrajState::ERR;
  }
  mprintf("\tRemoving data sets");
  std::string setSelectArg = argIn.GetStringNext();
  if (setSelectArg.empty())
    setSelectArg.assign("*");
  else
    mprintf(" within selection '%s'", setSelectArg.c_str());
  mprintf(" %s\n", status.c_str());
  DataSetList tempDSL = State.DSL().GetMultipleSets( setSelectArg );
  if (tempDSL.empty()) {
    mprinterr("Error: No data sets selected.\n");
    return CpptrajState::ERR;
  }
  // Remove sets
  unsigned int Nremoved = 0;
  if ( criterion == AVERAGE ) {
    if (!validDouble( val1 )) {
      mprinterr("Error: '%s' is not a valid number\n", val1.c_str());
      return CpptrajState::ERR;
    }
    double d_val1 = convertToDouble( val1 );
    double d_val2  = d_val1;
    if (!val2.empty()) {
      if (!validDouble( val2 )) {
        mprinterr("Error: '%s' is not a valid number\n", val2.c_str());
        return CpptrajState::ERR;
      }
      d_val2 = convertToDouble( val2 );
    }
    for (DataSetList::const_iterator ds = tempDSL.begin(); ds != tempDSL.end(); ++ds)
    {
      if ( (*ds)->Group() != DataSet::SCALAR_1D )
        mprintf("Warning: '%s' is not a valid data set for 'average' criterion.\n", (*ds)->legend());
      else {
        DataSet_1D const& ds1 = static_cast<DataSet_1D const&>( *(*ds) );
        double avg = ds1.Avg();
        bool remove = false;
        switch (select) {
          case EQUAL        : remove = (avg == d_val1); break;
          case NOT_EQUAL    : remove = (avg != d_val1); break;
          case LESS_THAN    : remove = (avg < d_val1); break;
          case GREATER_THAN : remove = (avg > d_val1); break;
          case BETWEEN      : remove = (avg > d_val1 && avg < d_val2); break;
          case OUTSIDE      : remove = (avg < d_val1 || avg > d_val2); break;
          case UNKNOWN_S:
          case N_S      : return CpptrajState::ERR; // Sanity check
        }
        if (remove) {
          mprintf("\t  Removing set '%s' (avg is %g)\n", (*ds)->legend(), avg);
          State.RemoveDataSet( *ds );
          ++Nremoved;
        }
      }
    }
  } else if ( criterion == SIZE ) {
    if (!validInteger( val1 )) {
      mprinterr("Error: '%s' is not a valid number\n", val1.c_str());
      return CpptrajState::ERR;
    }
    unsigned int i_val1 = (unsigned int)convertToInteger( val1 );
    unsigned int i_val2 = i_val1;
    if (!val2.empty()) {
      if (!validInteger( val2 )) {
        mprinterr("Error: '%s' is not a valid number\n", val2.c_str());
        return CpptrajState::ERR;
      }
      i_val2 = convertToInteger( val2 );
    }
    for (DataSetList::const_iterator ds = tempDSL.begin(); ds != tempDSL.end(); ++ds)
    {
      unsigned int size = (*ds)->Size();
      bool remove = false;
      switch ( select ) {
        case EQUAL        : remove = (size == i_val1); break;
        case NOT_EQUAL    : remove = (size != i_val1); break;
        case LESS_THAN    : remove = (size < i_val1); break;
        case GREATER_THAN : remove = (size > i_val1); break;
        case BETWEEN      : remove = (size > i_val1 && size < i_val2); break;
        case OUTSIDE      : remove = (size < i_val1 || size > i_val2); break;
        case UNKNOWN_S:
        case N_S      : return CpptrajState::ERR; // Sanity check
      }
      if (remove) {
        mprintf("\t  Removing set '%s' (size is %u)\n", (*ds)->legend(), size);
        State.RemoveDataSet( *ds );
        ++Nremoved;
      }
    }
  } else if ( criterion == SMODE ) {
    MetaData::scalarMode mode_val = MetaData::ModeFromKeyword( val1 );
    if (mode_val == MetaData::UNKNOWN_MODE) {
      mprinterr("Error: '%s' is not a valid mode.\n", val1.c_str());
      return CpptrajState::ERR;
    }
    for (DataSetList::const_iterator ds = tempDSL.begin(); ds != tempDSL.end(); ++ds)
    {
      bool remove = false;
      MetaData::scalarMode mode = (*ds)->Meta().ScalarMode();
      if      (select == EQUAL    ) remove = ( mode == mode_val );
      else if (select == NOT_EQUAL) remove = ( mode != mode_val );
      else return CpptrajState::ERR; // Sanity check
      if (remove) {
        mprintf("\t  Removing set '%s' (mode is '%s')\n", (*ds)->legend(), MetaData::ModeString(mode));
        State.RemoveDataSet( *ds );
        ++Nremoved;
      }
    }
  } else if ( criterion == STYPE ) {
    MetaData::scalarType type_val = MetaData::TypeFromKeyword( val1, MetaData::UNKNOWN_MODE );
    if (type_val == MetaData::UNDEFINED) {
      mprinterr("Error: '%s' is not a valid type.\n", val1.c_str());
      return CpptrajState::ERR;
    }
    for (DataSetList::const_iterator ds = tempDSL.begin(); ds != tempDSL.end(); ++ds)
    {
      bool remove = false;
      MetaData::scalarType type = (*ds)->Meta().ScalarType();
      if      (select == EQUAL    ) remove = ( type == type_val );
      else if (select == NOT_EQUAL) remove = ( type != type_val );
      else return CpptrajState::ERR; // Sanity check
      if (remove) {
        mprintf("\t  Removing set '%s' (typeis '%s')\n", (*ds)->legend(), MetaData::TypeString(type));
        State.RemoveDataSet( *ds );
        ++Nremoved;
      }
    }
  } else {
    mprinterr("Internal Error: Criterion not yet implemented.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tRemoved %u of %zu sets.\n", Nremoved, tempDSL.size());
  return CpptrajState::OK;
}

// Exec_DataSetCmd::MakeXY()
Exec::RetType Exec_DataSetCmd::MakeXY(CpptrajState& State, ArgList& argIn) {
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
  return CpptrajState::OK;
}

// Exec_DataSetCmd::Make2D()
Exec::RetType Exec_DataSetCmd::Make2D(CpptrajState& State, ArgList& argIn) {
  std::string name = argIn.GetStringKey("name");
  int ncols = argIn.getKeyInt("ncols", 0);
  int nrows = argIn.getKeyInt("nrows", 0);
  if (ncols <= 0 || nrows <= 0) {
    mprinterr("Error: Must specify both ncols and nrows\n");
    return CpptrajState::ERR;
  }
  DataSet* ds1 = State.DSL().GetDataSet( argIn.GetStringNext() );
  if (ds1 == 0) return CpptrajState::ERR;
  if (ds1->Ndim() != 1) {
    mprinterr("Error: make2d only works for 1D data sets.\n");
    return CpptrajState::ERR;
  }
  if (nrows * ncols != (int)ds1->Size()) {
    mprinterr("Error: Size of '%s' (%zu) != nrows X ncols.\n", ds1->legend(), ds1->Size());
    return CpptrajState::ERR;
  }
  if (name.empty())
    name = State.DSL().GenerateDefaultName("make2d");
  MetaData md(name, MetaData::M_MATRIX);
  DataSet* ds3 = State.DSL().AddSet( DataSet::MATRIX_DBL, md );
  
  if (ds3 == 0) return CpptrajState::ERR;
  mprintf("\tConverting values from 1D set '%s' to 2D matrix '%s' with %i cols, %i rows.\n",
          ds1->legend(), ds3->legend(), ncols, nrows);
  DataSet_1D const& data = static_cast<DataSet_1D const&>( *ds1 );
  DataSet_MatrixDbl& matrix = static_cast<DataSet_MatrixDbl&>( *ds3 );
  if (matrix.Allocate2D( ncols, nrows )) return CpptrajState::ERR;
  for (unsigned int idx = 0; idx != data.Size(); idx++)
    matrix.AddElement( data.Dval(idx) );
  return CpptrajState::OK;
}

// Exec_DataSetCmd::Filter()
Exec::RetType Exec_DataSetCmd::Filter(CpptrajState& State, ArgList& argIn) {
  std::string name = argIn.GetStringKey("name");
  int rowmin = argIn.getKeyInt("rowmin", -1);
  int rowmax = argIn.getKeyInt("rowmax", -1);
  int colmin = argIn.getKeyInt("colmin", -1);
  int colmax = argIn.getKeyInt("colmax", -1);
  
  DataSet* ds1 = State.DSL().GetDataSet( argIn.GetStringNext() );
  if (ds1 == 0) return CpptrajState::ERR;
  if ( ds1->Group() == DataSet::SCALAR_1D ) {
    mprinterr("Error: Not yet set up for 1D sets.\n");
    return CpptrajState::ERR;
  } else if (ds1->Group() == DataSet::MATRIX_2D) {
    DataSet_2D const& matrixIn = static_cast<DataSet_2D const&>( *ds1 );
    if (rowmin < 0) rowmin = 0;
    if (rowmax < 0) rowmax = matrixIn.Nrows();
    int nrows = rowmax - rowmin;
    if (nrows < 1) {
      mprinterr("Error: Keeping less than 1 row.\n");
      return CpptrajState::ERR;
    } else if (nrows > (int)matrixIn.Nrows())
      nrows = matrixIn.Nrows();
    if (colmin < 0) colmin = 0;
    if (colmax < 0) colmax = matrixIn.Ncols();
    int ncols = colmax - colmin;
    if (ncols < 1) {
      mprinterr("Error: Keeping less than 1 column.\n");
      return CpptrajState::ERR;
    } else if (ncols > (int)matrixIn.Ncols())
      ncols = matrixIn.Ncols();
    mprintf("\tMatrix to filter: %s\n", ds1->legend());
    mprintf("\tKeeping rows >= %i and < %i\n", rowmin, rowmax);
    mprintf("\tKeeping cols >= %i and < %i\n", colmin, colmax);
    mprintf("\tCreating new matrix with %i rows and %i columns.\n", nrows, ncols);
    DataSet* ds3 = State.DSL().AddSet( DataSet::MATRIX_DBL, name, "make2d" );
    if (ds3 == 0) return CpptrajState::ERR;
    DataSet_MatrixDbl& matrixOut = static_cast<DataSet_MatrixDbl&>( *ds3 );
    matrixOut.Allocate2D(ncols, nrows);
    matrixOut.SetDim( Dimension::X, Dimension(matrixIn.Dim(0).Coord(colmin),
                                              matrixIn.Dim(0).Step(),
                                              matrixIn.Dim(0).Label()) );
    matrixOut.SetDim( Dimension::Y, Dimension(matrixIn.Dim(1).Coord(rowmin),
                                              matrixIn.Dim(1).Step(),
                                              matrixIn.Dim(1).Label()) );
    for (int row = 0; row < (int)matrixIn.Nrows(); row++)
    {
      if (row >= rowmin && row < rowmax)
      {
        for (int col = 0; col < (int)matrixIn.Ncols(); col++)
        {
          if (col >= colmin && col < colmax)
          {
            double val = matrixIn.GetElement(col, row);
            matrixOut.SetElement( col-colmin, row-rowmin, val );
            //mprintf("DEBUG:\tmatrixIn(%i, %i) = %f  to matrixOut(%i, %i)\n",
            //        col, row, val, col-colmin, row-rowmin);
          }
        }
      }
    }
  }
  return CpptrajState::OK;
}

// Exec_DataSetCmd::Concatenate()
Exec::RetType Exec_DataSetCmd::Concatenate(CpptrajState& State, ArgList& argIn) {
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
      if ( (*ds)->Group() != DataSet::SCALAR_1D )
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
  return CpptrajState::OK;
}

void Exec_DataSetCmd::Help_ChangeDim() {
  mprintf("  dim {xdim|ydim|zdim|ndim <#>} [label <label>] [min <min>] [step <step>]\n"
          "    Change specified dimension in set(s).\n");
}

// Exec_DataSetCmd::ChangeDim()
Exec::RetType Exec_DataSetCmd::ChangeDim(CpptrajState const& State, ArgList& argIn) {
  int ndim = -1;
  if (argIn.hasKey("xdim"))
    ndim = 0;
  else if (argIn.hasKey("ydim"))
    ndim = 1;
  else if (argIn.hasKey("zdim"))
    ndim = 2;
  else
    ndim = argIn.getKeyInt("ndim", -1);
  if (ndim < 0) {
    mprinterr("Error: Specify xdim/ydim/zdim or dimension number with ndim.\n");
    return CpptrajState::ERR;
  }
  if (ndim < 3) {
    static const char DIMSTR[3] = { 'X', 'Y', 'Z' };
    mprintf("\tChanging the following in the %c dimension:\n", DIMSTR[ndim]);
  } else
    mprintf("\tChanging the following in dimension %i\n", ndim);

  bool changeLabel, changeMin, changeStep;
  std::string label;
  double min = 0.0;
  double step = 0.0;
  if (argIn.Contains("label")) {
    label = argIn.GetStringKey("label");
    changeLabel = true;
    mprintf("\tNew Label: %s\n", label.c_str());
  } else
    changeLabel = false;
  if (argIn.Contains("step")) {
    step = argIn.getKeyDouble("step", 0.0);
    changeStep = true;
    mprintf("\tNew step: %g\n", step);
  } else
    changeStep = false;
  if (argIn.Contains("min")) {
    min = argIn.getKeyDouble("min", 0.0);
    changeMin = true;
    mprintf("\tNew min: %g\n", min);
  } else
    changeMin = false;
  // Loop over all DataSet arguments 
  std::string ds_arg = argIn.GetStringNext();
  while (!ds_arg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( ds_arg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
    {
      if (ndim < (int)(*ds)->Ndim()) {
        mprintf("\t%s\n", (*ds)->legend());
        Dimension dim = (*ds)->Dim(ndim);
        if (changeLabel) dim.SetLabel( label );
        if (changeMin)   dim.ChangeMin( min );
        if (changeStep)  dim.ChangeStep( step );
        (*ds)->SetDim(ndim, dim);
      } else
        mprintf("Warning: Set '%s' has fewer then %i dimensions - skipping.\n",
                (*ds)->legend(), ndim);
    }
    ds_arg = argIn.GetStringNext();
  }
  return CpptrajState::OK;
}

// Exec_DataSetCmd::ChangeModeType()
Exec::RetType Exec_DataSetCmd::ChangeModeType(CpptrajState const& State, ArgList& argIn) {
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
      if (dmode != MetaData::UNKNOWN_MODE) {
        // Warn if mode does not seem appropriate for the data set type.
        if ( dmode >= MetaData::M_DISTANCE &&
             dmode <= MetaData::M_RMS &&
             (*ds)->Group() != DataSet::SCALAR_1D )
          mprintf("Warning: '%s': Expected scalar 1D data set type for mode '%s'\n",
                  (*ds)->legend(), MetaData::ModeString(dmode));
        else if ( dmode == MetaData::M_VECTOR &&
                  (*ds)->Group() != DataSet::VECTOR_1D )
          mprintf("Warning: '%s': Expected vector data set type for mode '%s'\n",
                  (*ds)->legend(), MetaData::ModeString(dmode));
        else if ( dmode == MetaData::M_MATRIX &&
                  (*ds)->Group() != DataSet::MATRIX_2D )
          mprintf("Warning: '%s': Expected 2D matrix data set type for mode '%s'\n",
                  (*ds)->legend(), MetaData::ModeString(dmode));
      }
      if ( dtype == MetaData::NOE ) (*ds)->AssociateData( &noeData );
      mprintf("\t\t'%s'\n", (*ds)->legend());
      MetaData md = (*ds)->Meta();
      md.SetScalarMode( dmode );
      md.SetScalarType( dtype );
      (*ds)->SetMeta( md );
    }
    ds_arg = argIn.GetStringNext();
  }
  return CpptrajState::OK;
}

// Exec_DataSetCmd::Help_InvertSets()
void Exec_DataSetCmd::Help_InvertSets() {
  mprintf("  invert <set arg0> ... name <new name> [legendset <set>]\n"
          "    Given M input sets of length N, create N new sets of length M by\n"
          "    inverting the input sets.\n");
}

/** Syntax: dataset invert <set arg0> ... name <new name> */
Exec::RetType Exec_DataSetCmd::InvertSets(CpptrajState& State, ArgList& argIn) {
  DataSetList& DSL = State.DSL();
  // Get keywords
  DataSet* inputNames = 0;
  std::string dsname = argIn.GetStringKey("legendset");
  if (!dsname.empty()) {
    inputNames = DSL.GetDataSet( dsname );
    if (inputNames == 0) {
      mprinterr("Error: Name set '%s' not found.\n", dsname.c_str());
      return CpptrajState::ERR;
    }
    if (inputNames->Type() != DataSet::STRING) {
      mprinterr("Error: Set '%s' does not contain strings.\n", inputNames->legend());
      return CpptrajState::ERR;
    }
    mprintf("\tUsing names from set '%s' as legends for inverted sets.\n", inputNames->legend());
  }
  dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: 'invert' requires that 'name <new set name>' be specified.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tNew sets will be named '%s'\n", dsname.c_str());
  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );
  if (outfile != 0)
    mprintf("\tNew sets will be output to '%s'\n", outfile->DataFilename().full());
  // TODO determine type some other way
  DataSet::DataType outtype = DataSet::DOUBLE;
  // Get input DataSets
  std::vector<DataSet_1D*> input_sets; 
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    DataSetList sets = DSL.GetMultipleSets( dsarg );
    for (DataSetList::const_iterator ds = sets.begin(); ds != sets.end(); ++ds)
    {
      if ( (*ds)->Group() != DataSet::SCALAR_1D ) {
        mprintf("Warning: '%s': Inversion only supported for 1D scalar data sets.\n",
                (*ds)->legend());
      } else {
        if (!input_sets.empty()) {
          if ( (*ds)->Size() != input_sets.back()->Size() ) {
            mprinterr("Error: Set '%s' has different size (%zu) than previous set (%zu)\n",
                      (*ds)->legend(), (*ds)->Size(), input_sets.back()->Size());
            return CpptrajState::ERR;
          }
        }
        input_sets.push_back( (DataSet_1D*)*ds );
      }
    }
    dsarg = argIn.GetStringNext();
  }
  if (input_sets.empty()) {
    mprinterr("Error: No sets selected.\n");
    return CpptrajState::ERR;
  }
  if (inputNames != 0 && inputNames->Size() != input_sets.front()->Size()) {
    mprinterr("Error: Name set '%s' size (%zu) differs from # data points (%zu).\n",
              inputNames->legend(), inputNames->Size(), input_sets.front()->Size());
    return CpptrajState::ERR;
  }
  mprintf("\t%zu input sets; creating %zu output sets.\n",
          input_sets.size(), input_sets.front()->Size());
  // Need an output data set for each point in input sets
  std::vector<DataSet*> output_sets;
  int column = 1;
  for (int idx = 0; idx != (int)input_sets[0]->Size(); idx++, column++) {
    DataSet* ds = 0;
    ds = DSL.AddSet(outtype, MetaData(dsname, column));
    if (ds == 0) return CpptrajState::ERR;
    if (inputNames != 0)
      ds->SetLegend( (*((DataSet_string*)inputNames))[idx] );
    output_sets.push_back( ds );
    if (outfile != 0) outfile->AddDataSet( ds );
  }
  // Create a data set containing names of each input data set
  DataSet* nameset = DSL.AddSet(DataSet::STRING, MetaData(dsname, column));
  if (nameset == 0) return CpptrajState::ERR;
  if (inputNames != 0)
    nameset->SetLegend("Names");
  if (outfile != 0) outfile->AddDataSet( nameset );
  // Populate output data sets
  for (int jdx = 0; jdx != (int)input_sets.size(); jdx++)
  {
    DataSet_1D const& INP = static_cast<DataSet_1D const&>( *(input_sets[jdx]) );
    nameset->Add( jdx, INP.legend() );
    for (unsigned int idx = 0; idx != INP.Size(); idx++)
    {
      double dval = INP.Dval( idx );
      output_sets[idx]->Add( jdx, &dval );
    }
  }

  return CpptrajState::OK;
}

// Exec_DataSetCmd::Help_Shift()
void Exec_DataSetCmd::Help_Shift() {
  mprintf("  shift [above <value> by <offset>] [below <value> by <offset>] <set arg0> ...\n"
          "    Shift the values in given set(s) by an offset if the points meet\n"
          "    certain criteria.\n");
}

/** Apply an offset to data elements when a certain criterion is met. */
Exec::RetType Exec_DataSetCmd::ShiftData(CpptrajState& State, ArgList& argIn) {
  std::vector<SelectType> criteria;
  std::vector<double> vals;
  std::vector<double> offsets;
  if (argIn.Contains("below")) {
    criteria.push_back(LESS_THAN);
    vals.push_back( argIn.getKeyDouble("below", 0) );
    if (!argIn.Contains("by")) {
      mprinterr("Error: 'by <offset>' argument missing for 'below'\n");
      return CpptrajState::ERR;
    }
    offsets.push_back( argIn.getKeyDouble("by", 0) );
  }
  if (argIn.Contains("above")) {
    criteria.push_back(GREATER_THAN);
    vals.push_back( argIn.getKeyDouble("above", 0) );
    if (!argIn.Contains("by")) {
      mprinterr("Error: 'by <offset>' argument missing for 'above'\n");
      return CpptrajState::ERR;
    }
    offsets.push_back( argIn.getKeyDouble("by", 0) );
  }
  if (criteria.empty()) {
    mprinterr("Error: shift requires 'below <value> by <offset>' or 'above <value> by <offset>\n");
    return CpptrajState::ERR;
  }
  // Sanity checks
  if (criteria.size() != vals.size()) {
    mprinterr("Error: Missing values for above/below; have %zu, expected %zu.\n",
              vals.size(), criteria.size());
    return CpptrajState::ERR;
  }
  if (criteria.size() != offsets.size()) {
    mprinterr("Error: Missing offsets for above/below; have %zu, expected %zu.\n",
              offsets.size(), criteria.size());
    return CpptrajState::ERR;
  }
  for (unsigned int ii = 0; ii < criteria.size(); ii++)
  {
    if (criteria[ii] == LESS_THAN)
      mprintf("\tValues below %g will be shifted by %g\n", vals[ii], offsets[ii]);
    else if (criteria[ii] == GREATER_THAN)
      mprintf("\tValues above %g will be shifted by %g\n", vals[ii], offsets[ii]);
  }
  // Loop over all DataSet arguments 
  std::string ds_arg = argIn.GetStringNext();
  if (ds_arg.empty()) {
    mprinterr("Error: No sets specfied.\n");
    return CpptrajState::ERR;
  }
  while (!ds_arg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( ds_arg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
    {
      if ((*ds)->Group() != DataSet::SCALAR_1D) {
        mprintf("Warning: Set '%s' is not scalar 1D, skipping.\n", (*ds)->legend());
      } else {
        DataSet_1D& set = static_cast<DataSet_1D&>( *(*ds) );
        mprintf("\tModifying set: %s\n", set.legend());
        // Check value against all criteria, modify by offset if necessary
        for (unsigned int idx = 0; idx < set.Size(); idx++)
        {
          double dval = set.Dval(idx);
          //mprintf("DBG: %u %g", idx, dval);
          for (unsigned int ii = 0; ii < criteria.size(); ii++)
          {
            if (criteria[ii] == LESS_THAN && dval < vals[ii])
              dval += offsets[ii];
            else if (criteria[ii] == GREATER_THAN && dval > vals[ii])
              dval += offsets[ii];
          }
          //mprintf(" %g\n", dval);
          set.SetY( idx, dval );
        } // END loop over set values
        modifiedSets_.push_back( *ds );
      }
    } // END loop over sets
    ds_arg = argIn.GetStringNext();
  } // END loop over data set args

  return CpptrajState::OK;
}
