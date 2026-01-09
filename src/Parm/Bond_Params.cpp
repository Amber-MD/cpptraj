#include "Bond_Params.h"
#include "../CpptrajStdio.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
#include <algorithm> // find

/** Create parameters for given bond based on element types. */
void Cpptraj::Parm::GenerateBondParam(BondParmArray& bondparm, BondType& bnd, BP_mapType& bpMap, std::vector<Atom> const& atoms)
{
  unsigned int bp_idx;
  Atom::AtomicElementType a1Elt = atoms[bnd.A1()].Element();
  Atom::AtomicElementType a2Elt = atoms[bnd.A2()].Element();
  std::set<Atom::AtomicElementType> Eset;
  Eset.insert( a1Elt );
  Eset.insert( a2Elt );
  // Has this bond parameter been defined?
  BP_mapType::iterator bp = std::find(bpMap.begin(), bpMap.end(), Eset);
  if (bp == bpMap.end()) { // Bond parameter Not defined
    bp_idx = bondparm.size();
    bpMap.push_back( Eset );
    bondparm.push_back( BondParmType(0.0, Atom::GetBondLength(a1Elt, a2Elt)) );
  } else
    bp_idx = bp - bpMap.begin();
  //mprintf("DEBUG:\t\t%i:[%s] -- %i:[%s] Cut=%f BPidx=%u\n",
  //        bnd.A1()+1, atoms_[bnd.A1()].c_str(), bnd.A2()+1, atoms_[bnd.A2()].c_str(),
  //        bondparm_[bp_idx].Req(), bp_idx);
  bnd.SetIdx( bp_idx );
}

/** Fill in bond parameters based on atomic element types. */
void Cpptraj::Parm::GenerateBondParameters(Topology& topOut)
{
  mprintf("\tDetermining bond length parameters from element types for '%s'.\n", topOut.c_str());
  BondParmArray& bondparm = topOut.ModifyBondParm();
  bondparm.clear();
  // Hold indices into bondparm for unique element pairs
  BP_mapType bpMap;
  for (BondArray::iterator bnd = topOut.ModifyBondsH().begin(); bnd != topOut.ModifyBondsH().end(); ++bnd)
    GenerateBondParam( bondparm, *bnd, bpMap, topOut.Atoms() ); 
  for (BondArray::iterator bnd = topOut.ModifyBonds().begin();  bnd != topOut.ModifyBonds().end(); ++bnd)
    GenerateBondParam( bondparm, *bnd, bpMap, topOut.Atoms() );
}
