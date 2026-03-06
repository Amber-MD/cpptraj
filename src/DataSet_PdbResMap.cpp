#include "DataSet_PdbResMap.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;
/// CONSTRUCTOR
DataSet_PdbResMap::DataSet_PdbResMap() :
  // 0 dim indicates DataSet-specific write
  DataSet(PDBRESMAP, GENERIC, TextFormat(TextFormat::STRING, 12, 0), 0)
{

}

/** Add a PDB residue to unit name mapping. */
int DataSet_PdbResMap::AddPdbResMap(PdbResMapType const& prmIn) {
  MapType::iterator it = pdbResMap_.lower_bound( prmIn.PdbName() );
  if (it == pdbResMap_.end() || it->first != prmIn.PdbName())
  {
    // New PDB residue name.
    PairType pdbUnitPair( prmIn.PdbName(), Sarray(3, "") );
    pdbUnitPair.second[(int)prmIn.TermType()] = prmIn.UnitName();
    it = pdbResMap_.insert( it, pdbUnitPair );
  } else {
    // Is this a terminal type that is not yet filled?
    if ( it->second[(int)prmIn.TermType()].empty() )
      it->second[(int)prmIn.TermType()] = prmIn.UnitName();
    else {
      mprintf("Warning: PDB residue %s is already mapped to unit %s for terminal type %s\n",
              *(it->first), it->second[(int)prmIn.TermType()].c_str(), Structure::terminalStr(prmIn.TermType()));
    }
  }

  return 0;
}

/** \return Unit name based on PDB residue name and terminal type. */
std::string DataSet_PdbResMap::FindUnitName(NameType const& pdbName, Cpptraj::Structure::TerminalType termType)
const
{
  MapType::const_iterator it = pdbResMap_.find( pdbName );
  if (it == pdbResMap_.end())
    return std::string("");
  else
    return it->second[(int)termType];
}

/** Print PDB residue to unit mapping to STDOUT. */
void DataSet_PdbResMap::PrintPdbResMap() const {
  mprintf("PDB to unit residue mapping:\n");
  mprintf("\t%8s : %8s %8s %8s\n", "PDB", "Begin", "Nonterm", "End");
  for (MapType::const_iterator it = pdbResMap_.begin(); it != pdbResMap_.end(); ++it)
  {
    mprintf("\t%8s : %8s %8s %8s\n", *(it->first),
            it->second[0].c_str(), it->second[1].c_str(), it->second[2].c_str());
  }
}
