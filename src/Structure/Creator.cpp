#include "Creator.h"
#include "GenerateConnectivityArrays.h" // For setting atom scan direction
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h" // TODO new coords type
#include "../DataSet_NameMap.h"
#include "../DataSet_Parameters.h"
#include "../DataSet_PdbResMap.h"
#include "../DataSetList.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Creator::Creator() :
  mainParmSet_(0),
  pdbResidueMap_(0),
  debug_(0),
  free_parmset_mem_(false)
{}

/** CONSTRUCTOR */
Creator::Creator(int d) :
  mainParmSet_(0),
  debug_(d),
  free_parmset_mem_(false)
{}

/** DESTRUCTOR */
Creator::~Creator() {
  if (mainParmSet_ != 0 && free_parmset_mem_)
    delete mainParmSet_;
}

const char* Creator::parm_keywords_ = "parmset <parameter setname>";

const char* Creator::template_keywords_ = "lib <template setname>";

const char* Creator::other_keywords_ = "atomscandir {f|b}";

/** Initialize */
int Creator::InitCreator(ArgList& argIn, DataSetList const& DSL, int debugIn)
{
  t_total_.Start();
  debug_ = debugIn;

  // PDB residue map TODO handle multiple maps?
  pdbResidueMap_ = (DataSet_PdbResMap*)DSL.FindSetOfType( "*", DataSet::PDBRESMAP );
  if (pdbResidueMap_ != 0) {
    if (debug_ > 0) {
      mprintf("DEBUG: PDB residue map data set: %s\n", pdbResidueMap_->legend());
      pdbResidueMap_->PrintPdbResMap();
    }
  }

  // Atom scan direction
  std::string atomscandir = argIn.GetStringKey("atomscandir");
  if (!atomscandir.empty()) {
    if (atomscandir == "f")
      Cpptraj::Structure::SetAtomScanDirection(Cpptraj::Structure::SCAN_ATOMS_FORWARDS);
    else if (atomscandir == "b")
      Cpptraj::Structure::SetAtomScanDirection(Cpptraj::Structure::SCAN_ATOMS_BACKWARDS);
    else {
      mprinterr("Error: Unrecognized keyword for 'atomscandir' : %s\n", atomscandir.c_str());
      return 1;
    }
  }

  t_get_templates_.Start();
  if (getTemplates(argIn, DSL)) return 1;
  t_get_templates_.Stop();
  t_get_parameters_.Start();
  if (getParameterSets(argIn, DSL)) return 1;
  t_get_parameters_.Stop();
  UpdateTemplateElements();

  // Get any atom name maps
  DataSetList nameMapSets = DSL.GetSetsOfType("*", DataSet::NAMEMAP);
  for (DataSetList::const_iterator ds = nameMapSets.begin();
                                   ds != nameMapSets.end(); ++ds)
  {
    NameMaps_.push_back( static_cast<DataSet_NameMap*>( *ds ) );
    if (debug_ > 0)
      mprintf("DEBUG: Atom name map: %s\n", NameMaps_.back()->legend());
  }
  t_total_.Stop();
  return 0;
}

/** Write timing info to stdout. */
void Creator::TimingInfo(double total, int indent) const {
  t_total_.WriteTiming(indent, "Get templates/parms:", total);
  t_get_templates_.WriteTiming (indent+1, "Get templates  :", t_total_.Total());
  t_get_parameters_.WriteTiming(indent+1, "Get parameters :", t_total_.Total());
}

/** Update template atom elements from atom types in parameter set. */
void Creator::UpdateTemplateElements() const {
  if (mainParmSet_ == 0) return;

  for (Carray::const_iterator cit = Templates_.begin(); cit != Templates_.end(); ++cit)
  {
    DataSet_Coords* crd = cit->second;
    if (debug_ > 0)
      mprintf("DEBUG: Updating atom elements in '%s'\n", crd->legend());
    // Loop over template atoms
    Topology& templateTop = *(crd->TopPtr());
    for (int at = 0; at != templateTop.Natom(); at++) {
      Cpptraj::Parm::ParmHolder<AtomType>::const_iterator it;
      if (templateTop[at].HasType()) {
        it = mainParmSet_->AT().GetParam( TypeNameHolder(templateTop[at].Type()) );
        if (it != mainParmSet_->AT().end()) {
          if (debug_ > 1)
            mprintf("DEBUG:\t\tSetting atom %s element to %s\n", *(templateTop[at].Name()), it->second.EltStr());
          templateTop.SetAtom(at).SetElementFromSymbol( it->second.EltStr()[0],
                                                        it->second.EltStr()[1] );
        }
      }
    } // END loop over template atoms
  } // END loop over templates
}

/** Try to identify residue template DataSet from the given name;
  * could be data set name or residue name.
  */
DataSet_Coords* Creator::IdTemplateFromName(std::string const& nameIn)
const
{
  Carray::const_iterator it = Templates_.find( nameIn );
  if (it == Templates_.end()) {
    mprintf("Warning: No template found named '%s'\n", nameIn.c_str());
    return 0;
  }
  return it->second;
  //MetaData::SearchString search( nameIn );
  //for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
  //  if ((*it)->Meta().Match_WildCard(search)) {
  //    return *it;
  //  }
  //}
  // No aspect. Convert to NameType TODO check for truncation
  //NameType rname( nameIn );
  //return IdTemplateFromResname( NameType(nameIn), Cpptraj::Structure::NON_TERMINAL );
}

/** Try to identify residue template DataSet from the given residue
  * name (from e.g. the PDB/Mol2/etc file).
  */
DataSet_Coords* Creator::IdTemplateFromResname(NameType const& rname,
                                              TerminalType termType)
const
{
  std::vector<DataSet_Coords*> Out;

  // See if a PDB residue map exists.
  std::string targetUnitName;
  if (pdbResidueMap_ != 0) {
    targetUnitName = pdbResidueMap_->FindUnitName(rname, termType);
//    if (!targetUnitName.empty())
//       mprintf("DEBUG: Found mapped name for '%s' (%s) -> '%s'\n", *rname, Structure::terminalStr(termType), targetUnitName.c_str());
  }
  if (targetUnitName.empty()) {
    targetUnitName = rname.Truncated();
//    mprintf("DEBUG: Target unit name: %s\n", targetUnitName.c_str());
  }
  // Most residue templates have name in aspect currently.
  // Residue templates loaded separately (via a mol2) may just
  // have name.
//  for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
//    if ((*it)->Meta().Aspect().empty()) {
//      if ((*it)->Meta().Name() == targetUnitName)
//        Out.push_back( *it );
//    } else if ((*it)->Meta().Aspect() == targetUnitName)
//      Out.push_back( *it );
//  }
  Carray::const_iterator it = Templates_.find( targetUnitName );
  if (it == Templates_.end()) {
    mprintf("Warning: No template found named '%s'\n", targetUnitName.c_str());
    return 0;
  }
  return it->second;

/*
  //DataSet_Coords* out = 0;
  if (termType != Cpptraj::Structure::NON_TERMINAL) {
    // Looking for a terminal residue. Need to get sets with AssociatedData_ResId
    for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
      AssociatedData* ad = (*it)->GetAssociatedData( AssociatedData::RESID );
      if (ad != 0) {
        AssociatedData_ResId const& resid = static_cast<AssociatedData_ResId const&>( *ad );
        if (rname == resid.ResName() && termType == resid.TermType()) {
          //out = *it;
          //break;
          Out.push_back( *it );
        }
      }
    }
  }
  if (Out.empty()) {
  //if (out == 0) {
    // Terminal residue with alias not found or non-terminal residue.
    if (debug_ > 0 && termType != Cpptraj::Structure::NON_TERMINAL)
      mprintf("DEBUG: No aliased terminal residue found for '%s'\n", *rname);
    // Assume Coords set aspect is what we need
    for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
      if ( (*it)->Meta().Aspect().size() < NameType::max() ) {
        if ( rname == NameType( (*it)->Meta().Aspect() ) ) {
          //out = *it;
          //break;
          Out.push_back( *it );
        }
      }
    }
  }
  if (Out.empty()) {
  //if (out == 0) {
    // As a final attempt, just look for the name
    for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
      if ( (*it)->Meta().Name().size() < NameType::max() ) {
        if ( rname == NameType( (*it)->Meta().Name() ) ) {
          //out = *it;
          //break;
          Out.push_back( *it );
        }
      }
    }
  }
*/
/*  if (Out.empty()) return 0;
  if (Out.size() > 1) {
    mprintf("Warning: Multiple templates match '%s':", *rname);
    for (std::vector<DataSet_Coords*>::const_iterator it = Out.begin(); it != Out.end(); ++it)
      mprintf(" %s", (*it)->legend());
    mprintf("\n");
    mprintf("Warning: Using the last template loaded.\n");
  }
  //return out;
  return Out.back();*/
}

/// \return Aspect, or Name if Aspect is empty
static inline std::string getTemplateName(DataSet* ds)
{
  if (ds->Meta().Aspect().empty())
    return ds->Meta().Name();
  else
    return ds->Meta().Aspect();
}

/** Add coords set as a template */
void Creator::addCoordsAsTemplate(DataSet_Coords* ds) {
  std::string templateName = getTemplateName( ds );
  Carray::iterator unit = Templates_.lower_bound( templateName );
  if (unit == Templates_.end() || unit->first != templateName) {
    unit = Templates_.insert( unit, Cpair( templateName, ds ) );
  } else {
    mprintf("Warning: Replacing template %s with %s\n", unit->second->legend(), ds->legend());
    unit->second = ds;
  }
}

/** Get templates */
int Creator::getTemplates(ArgList& argIn, DataSetList const& DSL) {
  // Clear existing templates
  Templates_.clear();
  std::string lib = argIn.GetStringKey("lib");
  if (lib.empty()) {
    mprintf("\tNo template(s) specified with 'lib'; using any loaded templates.\n");
    DataSetList sets = DSL.SelectGroupSets( "*", DataSet::COORDINATES ); // TODO specific set type for units?
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
    {
      DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
      // Should only be a single residue FIXME need new set type
      if ( ds->Top().Nres() == 1 ) {
        addCoordsAsTemplate( ds );
        //Templates_.push_back( (DataSet_Coords*)(*it) );
      }
    }
  } else {
    while (!lib.empty()) {
      DataSetList sets = DSL.SelectGroupSets( lib, DataSet::COORDINATES ); // TODO specific set type for units?
      if (sets.empty()) {
        mprintf("Warning: No sets corresponding to '%s'\n", lib.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
        {
          // Should only be a single residue FIXME need new set type
          //DataSet_Coords const& ds = static_cast<DataSet_Coords const&>( *(*it) );
          //if ( ds.Top().Nres() == 1 )
          //  Templates_.push_back( (DataSet_Coords*)(*it) );
          addCoordsAsTemplate( static_cast<DataSet_Coords*>( *it ) );
        }
      }
      lib = argIn.GetStringKey("lib");
    }
  }
  if (!Templates_.empty()) {
    mprintf("\t%zu residue templates found:\n", Templates_.size());
    if (debug_ > 0) {
      for (Carray::const_iterator it = Templates_.begin(); it != Templates_.end(); ++it) {
        mprintf("\t%s", it->second->legend());
//        AssociatedData* ad = (*it)->GetAssociatedData( AssociatedData::RESID );
//        if (ad != 0) {
//          AssociatedData_ResId const& resid = static_cast<AssociatedData_ResId const&>( *ad );
//          resid.Ainfo();
//        }
        mprintf("\n");
      }
    }
  }

  return 0;
}

/** Get parameter sets. */
int Creator::getParameterSets(ArgList& argIn, DataSetList const& DSL) {
  // Clear any existing set
  if (mainParmSet_ != 0) {
    if (free_parmset_mem_) delete mainParmSet_;
  }
  mainParmSet_ = 0;
  free_parmset_mem_ = false;
  // Look for parmset args
  typedef std::vector<DataSet_Parameters*> Parray;
  Parray ParamSets;
  std::string parmset = argIn.GetStringKey("parmset");
  if (parmset.empty()) {
    mprintf("\tNo parameter set(s) specified with 'parmset'; using any loaded sets.\n");
    // See if there are any parameter sets.
    DataSetList sets = DSL.GetSetsOfType( "*", DataSet::PARAMETERS );
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
      ParamSets.push_back( (DataSet_Parameters*)(*it) );
  } else {
    while (!parmset.empty()) {
      DataSetList sets = DSL.GetSetsOfType( parmset, DataSet::PARAMETERS );
      if (sets.empty()) {
        mprintf("Warning: No parameter sets corresponding to '%s'\n", parmset.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
          ParamSets.push_back( (DataSet_Parameters*)(*it) );
      }
      parmset = argIn.GetStringKey("parmset");
    }
  }
  //if (ParamSets.empty()) {
  //  mprinterr("Error: No parameter sets.\n");
  //  return CpptrajState::ERR;
  //}
  if (!ParamSets.empty()) {
    mprintf("\tParameter sets:\n");
    for (Parray::const_iterator it = ParamSets.begin(); it != ParamSets.end(); ++it)
      mprintf("\t  %s\n", (*it)->legend());

    // Combine parameters if needed

    if (ParamSets.size() == 1)
      mainParmSet_ = ParamSets.front();
    else {
      free_parmset_mem_ = true;
      mprintf("\tCombining parameter sets.\n");
      Parray::const_iterator it = ParamSets.begin();
      mprintf("\t  Initial parameter set: %s\n", (*it)->legend());
      mainParmSet_ = new DataSet_Parameters( *(*it) );
      ++it;
      Cpptraj::Parm::ParameterSet::UpdateCount UC;
      for (; it != ParamSets.end(); ++it) {
        mprintf("\t  Adding parameter set: %s\n", (*it)->legend());
        mainParmSet_->UpdateParamSet( *(*it), UC, debug_, debug_+1 ); // Make it so verbosity is at least 1 to report overwritten params
      }
    }
  }
  return 0;
}

/** Get alias if present */
bool Creator::GetAlias(NameType& newName, NameType const& oldName)
const
{
  for (Narray::const_iterator it = NameMaps_.begin(); it != NameMaps_.end(); ++it)
  {
    if ((*it)->GetName( newName, oldName ))
      return true;
  }
  return false;
}

/** Count missing atoms from template */
int Creator::CountAtomsMissingFromTemplate(Topology const& topIn,
                                           int rnum,
                                           DataSet_Coords* resTemplate)
const
{
  int nTgtAtomsMissing = 0;
  Residue const& resIn = topIn.Res(rnum);
  // For each atom in topIn, find a template atom
  for (int itgt = resIn.FirstAtom(); itgt != resIn.LastAtom(); itgt++)
  {
    NameType const& tgtName = topIn[itgt].Name();
    //mprintf("DEBUG: Search for atom %s\n", *tgtName);
    bool found = false;
    // Check if this atom has an alias.
    NameType alias;
    bool has_alias = GetAlias( alias, tgtName );
//    if (creator.GetAlias( alias, tgtName )) {
//      mprintf("DEBUG: Atom %s alias is %s\n", *tgtName, *alias);
//    }
    // See if tgtName matches a reference (template) atom name.
    for (int iref = 0; iref != resTemplate->Top().Natom(); iref++)
    {
      NameType const& refName = resTemplate->Top()[iref].Name();
      if (refName == tgtName) {
        found = true;
        break;
      }
    }
    if (!found && has_alias) {
      // See if alias matches a reference (template) atom name.
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        NameType const& refName = resTemplate->Top()[iref].Name();
        if (refName == alias) {
          found = true;
          break;
        }
      } // END search template for alias
    } // END do alias search
    if (!found) {
      //mprintf("Warning: Input atom %s was not mapped to a template atom.\n",
      //        topIn.TruncAtomResNameOnumId( itgt ).c_str());
      nTgtAtomsMissing++;
    }
  }
  return nTgtAtomsMissing;
}

/** Map atoms in residue to template. */
std::vector<int> Creator::MapAtomsToTemplate(Topology const& topIn,
                                                int rnum,
                                                DataSet_Coords* resTemplate,
                                                std::vector<NameType>& sourceAtomNames,
                                                int& nTgtAtomsMissing)
const
{
  nTgtAtomsMissing = 0;
  std::vector<int> mapOut(resTemplate->Top().Natom(), -1);
  mapOut.reserve( resTemplate->Top().Natom() );
  Residue const& resIn = topIn.Res(rnum);
  // For each atom in topIn, find a template atom
  for (int itgt = resIn.FirstAtom(); itgt != resIn.LastAtom(); itgt++)
  {
    NameType const& tgtName = topIn[itgt].Name();
    //mprintf("DEBUG: Search for atom %s\n", *tgtName);
    bool found = false;
    // Check if this atom has an alias.
    NameType alias;
    bool has_alias = GetAlias( alias, tgtName );
//    if (creator.GetAlias( alias, tgtName )) {
//      mprintf("DEBUG: Atom %s alias is %s\n", *tgtName, *alias);
//    }
    // See if tgtName matches a reference (template) atom name.
    for (int iref = 0; iref != resTemplate->Top().Natom(); iref++)
    {
      NameType const& refName = resTemplate->Top()[iref].Name();
      if (refName == tgtName) {
        sourceAtomNames[itgt] = tgtName;
        mapOut[iref] = itgt;
        found = true;
        break;
      }
    }
    if (!found && has_alias) {
      // See if alias matches a reference (template) atom name.
      for (int iref = 0; iref != resTemplate->Top().Natom(); iref++) {
        NameType const& refName = resTemplate->Top()[iref].Name();
        if (refName == alias) {
          sourceAtomNames[itgt] = alias;
          mapOut[iref] = itgt;
          found = true;
          break;
        }
      } // END search template for alias
    } // END do alias search
    if (!found) {
      mprintf("Warning: Input atom %s was not mapped to a template atom.\n",
              topIn.TruncAtomResNameOnumId( itgt ).c_str());
      nTgtAtomsMissing++;
    }
  }
/*
  for (int iref = 0; iref != resTemplate->Top().Natom(); iref++)
  {
    // Find this atom name in topIn
    NameType const& refName = resTemplate->Top()[iref].Name();
    int iat = -1;
    for (int itgt = resIn.FirstAtom(); itgt != resIn.LastAtom(); itgt++) {
      if ( refName == topIn[itgt].Name() ) {
        iat = itgt;
        break;
      }
    }
    
    mapOut.push_back( iat );
  }*/
  return mapOut;
}
