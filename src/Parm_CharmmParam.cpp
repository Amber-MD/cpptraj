#include "Parm_CharmmParam.h"
#include "ParameterSet.h"
#include "CharmmParamFile.h"

static inline void AddBonds(ParameterSet& pout, std::vector<bool>& typeUsed,
                            std::vector<Atom> const& atoms,
                            BondParmArray const& BP, BondArray const& bonds)
{
  for (BondArray::const_iterator it = bonds.begin(); it != bonds.end(); ++it)
  {
    if (it->Idx() != -1) {
      if (!typeUsed[it->Idx()]) {
        AtomTypeHolder types(2);
        types.AddName( atoms[it->A1()].Type() );
        types.AddName( atoms[it->A2()].Type() );
        pout.BP().AddParm( types, BP[it->Idx()], false );
        typeUsed[it->Idx()] = true;
      }
    }
  }
}

static inline void AddAngles(ParameterSet& pout, std::vector<bool>& typeUsed,
                            std::vector<Atom> const& atoms,
                            AngleParmArray const& AP, AngleArray const& angles)
{
  for (AngleArray::const_iterator it = angles.begin(); it != angles.end(); ++it)
  {
    if (it->Idx() != -1) {
      if (!typeUsed[it->Idx()]) {
        AtomTypeHolder types(3);
        types.AddName( atoms[it->A1()].Type() );
        types.AddName( atoms[it->A2()].Type() );
        types.AddName( atoms[it->A3()].Type() );
        pout.AP().AddParm( types, AP[it->Idx()], false );
        typeUsed[it->Idx()] = true;
      }
    }
  }
}

static inline void AddDihedrals(ParameterSet& pout, std::vector<bool>& typeUsed,
                            std::vector<Atom> const& atoms,
                            DihedralParmArray const& DP, DihedralArray const& dihedrals)
{
  for (DihedralArray::const_iterator it = dihedrals.begin(); it != dihedrals.end(); ++it)
  {
    if (it->Idx() != -1) {
      if (!typeUsed[it->Idx()]) {
        AtomTypeHolder types(4);
        types.AddName( atoms[it->A1()].Type() );
        types.AddName( atoms[it->A2()].Type() );
        types.AddName( atoms[it->A3()].Type() );
        types.AddName( atoms[it->A4()].Type() );
        pout.DP().AddParm( types, DP[it->Idx()], false );
        typeUsed[it->Idx()] = true;
      }
    }
  }
}

int Parm_CharmmParam::WriteParm(FileName const& fname, Topology const& topOut) {
  // Convert topology to parameter set
  ParameterSet pout;

  for (Topology::atom_iterator it = topOut.begin(); it != topOut.end(); ++it)
  {
    pout.AT().AddAtomType( it->Type(), AtomType( it->Mass() ) );
  }

  std::vector<bool> typeUsed( topOut.BondParm().size(), false );
  AddBonds( pout, typeUsed, topOut.Atoms(), topOut.BondParm(), topOut.Bonds() );
  AddBonds( pout, typeUsed, topOut.Atoms(), topOut.BondParm(), topOut.BondsH() );

  typeUsed.assign( topOut.AngleParm().size(), false );
  AddAngles( pout, typeUsed, topOut.Atoms(), topOut.AngleParm(), topOut.Angles() );
  AddAngles( pout, typeUsed, topOut.Atoms(), topOut.AngleParm(), topOut.AnglesH() );

  typeUsed.assign( topOut.DihedralParm().size(), false );
  AddDihedrals( pout, typeUsed, topOut.Atoms(), topOut.DihedralParm(), topOut.Dihedrals() );
  AddDihedrals( pout, typeUsed, topOut.Atoms(), topOut.DihedralParm(), topOut.DihedralsH() );

  // Write out
  CharmmParamFile outfile;
  return outfile.WriteParams(pout, fname, debug_);
}
