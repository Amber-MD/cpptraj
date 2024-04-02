#include "InteractionData.h"
#include "CpptrajStdio.h"
#include "DataFile.h"
#include "DataSetList.h"

using namespace Cpptraj;

/** CONSTRUCTOR */
InteractionData::InteractionData()
{}

/** Add set tracking interaction between the two given indices. */
DataSet* InteractionData::AddInteractionSet(DataSetList& dsl, DataSet::DataType typeIn,
                                            MetaData const& metaIn, int i0, int i1,
                                            DataFile* outfile)
{
  // Check for existing interaction
  Ipair interaction;
  if (i1 < i0) {
    interaction.first = i1;
    interaction.second = i0;
  } else {
    interaction.first = i0;
    interaction.second = i1;
  }
  MapType::const_iterator it = setMap_.lower_bound(interaction);
  if (it == setMap_.end() || it->first != interaction)
  {
    // New interaction
    DataSet* ds = dsl.AddSet(typeIn, metaIn);
    if (ds == 0) {
      mprinterr("Error: Could not allocate data set %s\n", metaIn.PrintName().c_str());
      return 0;
    }
    if (outfile != 0) outfile->AddDataSet( ds );
    it = setMap_.insert(it, PairType(interaction, ds));
  }
  // TODO check existing set for issues?
  return it->second;
}
