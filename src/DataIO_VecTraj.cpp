#include "DataIO_VecTraj.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector.h"
#include "ParmFile.h"
#include "Trajout_Single.h"

// CONSTRUCTOR
DataIO_VecTraj::DataIO_VecTraj() :
  trajoutFmt_(TrajectoryFile::UNKNOWN_TRAJ),
  includeOrigin_(true)
{
  SetValid( DataSet::VECTOR );
}

void DataIO_VecTraj::WriteHelp() {
  mprintf("\t[trajfmt <format>] [parmout <file>] [noorigin]\n");
}

int DataIO_VecTraj::processWriteArgs(ArgList& argIn) {
  trajoutFmt_ = TrajectoryFile::WriteFormatFromString( argIn.GetStringKey("trajfmt"),
                                                       TrajectoryFile::AMBERTRAJ );
  parmoutName_ = argIn.GetStringKey("parmout");
  includeOrigin_ = !argIn.hasKey("noorigin");
  return 0;
}

int DataIO_VecTraj::WriteData(FileName const& fname, DataSetList const& SetList) {
  if (SetList.empty()) return 1;
  // Check input data sets
  int vec_size = -1;
  unsigned int num_no_origins = 0;
  typedef std::vector<DataSet_Vector*> Varray;
  Varray VecSets;
  VecSets.reserve( SetList.size() );
  for (DataSetList::const_iterator set = SetList.begin(); set != SetList.end(); ++set) {
    if ((*set)->Type() != DataSet::VECTOR)
      mprintf("Warning: Set '%s' is not a vector, skipping.\n", (*set)->legend());
    else {
      DataSet_Vector const& Vec = static_cast<DataSet_Vector const&>( *(*set) );
      if (vec_size == -1)
        vec_size = (int)Vec.Size();
      else if (vec_size != (int)Vec.Size()) {
        mprinterr("Error: Vector '%s' size %zu != first vector size %zu.\n"
                  "Error:   All vectors must have same size.\n",
                  Vec.legend(), Vec.Size(), vec_size);
        return 1;
      }
      if (!Vec.HasOrigins()) num_no_origins++;
      VecSets.push_back( (DataSet_Vector*)*set );
    }
  }
  if (VecSets.empty()) {
    mprinterr("Error: No vector data sets.\n");
    return 1;
  }
  if (num_no_origins == VecSets.size())
    includeOrigin_ = false;
  // Set up pseudo topology for all vectors
  Topology pseudo;
  BondArray bonds;
  // 1 pseudo bond type, Rk = 0.0, Req = 1.0 Ang.
  pseudo.AddBondParm( BondParmType(0.0, 1.0) );
  int natom = 0;
  for (unsigned int nres = 1; nres <= VecSets.size(); nres++) {
    Residue vec_res("VEC", nres, ' ', ' ');
    if (includeOrigin_)
      pseudo.AddTopAtom(Atom("OXYZ", 0), vec_res);
    pseudo.AddTopAtom(Atom("VXYZ", 0), vec_res);
    if (includeOrigin_) {
      pseudo.AddBond(natom, natom+1, 0); // Bond parm index 0
      natom += 2;
    } else
      natom++;
  }
  pseudo.CommonSetup();
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology( pseudo, parmoutName_, ParmFile::UNKNOWN_PARM, 0 )) {
      mprinterr("Error: Could not write pseudo topology to '%s'\n", parmoutName_.c_str());
      return 1;
    }
  }
  // Write out vectors
  Trajout_Single out;
  if (out.PrepareTrajWrite(fname, ArgList(), &pseudo, CoordinateInfo(),
                           vec_size, trajoutFmt_) == 0)
  {
    Frame outFrame(pseudo.Natom());
    for (int i = 0; i != vec_size; ++i) {
      outFrame.ClearAtoms();
      for (Varray::const_iterator set = VecSets.begin(); set != VecSets.end(); ++set) {
        DataSet_Vector const& Vec = static_cast<DataSet_Vector const&>( *(*set) );
        if (includeOrigin_) {
          Vec3 const& ovec = Vec.OXYZ(i);
          outFrame.AddVec3( ovec );
          outFrame.AddVec3( Vec[i] + ovec );
        } else
          outFrame.AddVec3( Vec[i] );
      }
      if (out.WriteSingle(i, outFrame)) return 1;
    }
    out.EndTraj();
  } else {
    mprinterr("Error: Could not set up '%s' for write.\n", fname.full());
    return 1;
  }
  return 0;
}
