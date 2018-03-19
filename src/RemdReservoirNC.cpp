#include "RemdReservoirNC.h"
#include "CpptrajStdio.h"
#include "NC_Routines.h"
#ifdef BINTRAJ
#include <netcdf.h>
#endif

/** Initialize NetCDF structure reservoir. */
int RemdReservoirNC::InitReservoir(FileName const& fnameIn, std::string const& titleIn,
                                   CoordinateInfo const& cinfoIn,
                                   int natomsIn, bool hasBins, double tempIn, int iseed)
{
# ifdef BINTRAJ
  CoordinateInfo cInfo(cinfoIn.TrajBox(), true, cinfoIn.HasVel(), cinfoIn.HasForce(), false);
  if (NC_create( fnameIn.Full(), NC_AMBERTRAJ, natomsIn, cInfo,
                 "Cpptraj generated structure reservoir", debug_ ))
    return 1;
  if (NC_createReservoir(hasBins, tempIn, iseed, eptotVID_, binsVID_))
    return 1;
  Coord_.resize( Ncatom3() );
  return 0;
# else
  mprinterr("Error: NetCDF reservoir requires NetCDF support. Recompile with -DBINTRAJ.\n");
  return 1;
# endif
}

/** Write a frame to structure reservoir. */
int RemdReservoirNC::WriteReservoir(unsigned int set, Frame const& frame, double energy, int bin)
{
# ifdef BINTRAJ
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  // Coords
  DoubleToFloat(&Coord_[0], frame.xAddress());
  if (NC::CheckErr(nc_put_vara_float(ncid_,coordVID_,start_,count_,&Coord_[0])) ) {
    mprinterr("Error: Netcdf writing reservoir coords %u\n",set);
    return 1;
  }
  // Velo
  if (velocityVID_ != -1) {
    if (frame.vAddress() == 0) { // TODO: Make it so this can NEVER happen.
      mprinterr("Error: Reservoir expects velocities, but no velocities in frame.\n");
      return 1;
    }
    DoubleToFloat(&Coord_[0], frame.vAddress());
    if (NC::CheckErr(nc_put_vara_float(ncid_,velocityVID_,start_,count_,&Coord_[0])) ) {
      mprinterr("Error: Netcdf writing reservoir velocities %i\n",set);
      return 1;
    }
  }
  // Eptot, bins
  if ( NC::CheckErr( nc_put_vara_double(ncid_,eptotVID_,start_,count_,&energy)) ) {
    mprinterr("Error: Writing eptot.\n");
    return 1;
  }
  if (binsVID_ != -1) {
    if ( NC::CheckErr( nc_put_vara_int(ncid_,binsVID_,start_,count_,&bin)) ) {
      mprinterr("Error: Writing bins.\n");
      return 1;
    }
  }
  // Write box
  if (cellLengthVID_ != -1) {
    count_[1] = 3;
    count_[2] = 0;
    if (NC::CheckErr(nc_put_vara_double(ncid_,cellLengthVID_,start_,count_,frame.bAddress())) ) {
      mprinterr("Error: Writing cell lengths.\n");
      return 1;
    }
    if (NC::CheckErr(nc_put_vara_double(ncid_,cellAngleVID_,start_,count_, frame.bAddress()+3)) ) {
      mprinterr("Error: Writing cell angles.\n");
      return 1;
    }
  }
  nc_sync(ncid_); // Necessary after every write??
  return 0;
# else
  return 1;
# endif
}
