// CharmmDcd
#include "Traj_CharmmDcd.h"
#include "CpptrajStdio.h"
#include <cstddef>

// CONSTRUCTOR
CharmmDcd::CharmmDcd() {
  dcdatom=0;
  dcdExtraBlock=false;
  dcd4D=false;
  istart=0;
  nsavc=0;
  namnf=0;
  nfreat=0;
  freeat=NULL;
  timestep=0;
}

// DESTRUCTOR
CharmmDcd::~CharmmDcd() {
  if (freeat!=NULL) delete[] freeat;
}

/* CharmmDcd::openTraj()
 */
int CharmmDcd::openTraj() {
  int err;

  err=0;
  switch (tfile->access) {
    case READ : err = tfile->OpenFile(); break;
    case APPEND :
      mprintf("Error: Append not supported for dcd files.\n");
      err=1;
      break;
    case WRITE :
      // Set up title
      if (title==NULL)
        this->SetTitle((char*)"Cpptraj generated mol2 file.\0");
      err = tfile->OpenFile();
      break;
  }
      
  return err;
}

/* CharmmDcd::closeTraj() {
 */
void CharmmDcd::closeTraj() {
  tfile->CloseFile();
}

// Convert 4 bytes to integer
/*static int byte_to_int(unsigned char byte[4]) {
  int val,temp;
  val=0;
  temp = (int) byte[0];
  temp <<= 24;
  val = val | temp;
  temp = (int) byte[1];
  temp <<= 16;
  val = val | temp;
  temp = (int) byte[2];
  temp <<= 8;
  val = val | temp;
  temp = (int) byte[3];
  val = val | temp;
  return val;
}*/

inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}

/* CharmmDcd::setupRead()
 */
int CharmmDcd::setupRead(AmberParm *trajParm) {
  int Frames=0;
  int dcdkey;
  char *dcdkey_char = (char*) &dcdkey;
  doublebyte LEbyte;
  doublebyte BEbyte;
  bool isBigEndian = false;
  bool is64bit = false;
  unsigned int readSize;
  unsigned char buffer[80];
  int ntitle;
  char dcdtitle[81];
  dcdtitle[80]='\0';

  if ( openTraj() ) return -1;

  // Step 1 - Determine endianness.
  // Read first 8 bytes - first number in dcd header should be 84.
  // If 32 bit the number is in the first 4 bytes, 64 bit is first
  // 8 bytes.
  tfile->IO->Read(LEbyte.c, sizeof(unsigned char), 8);

  dcdkey_char[0]='C';
  dcdkey_char[1]='O';
  dcdkey_char[2]='R';
  dcdkey_char[3]='D';

  // If the sum of the first 2 4-byte blocks is 84 this is a 
  // 64bit little endian
  if ( (LEbyte.i[0]+LEbyte.i[1]) == 84 ) {
    isBigEndian = false;
    is64bit = true;
  // If first 4-byte block is 84 and second is CORD this is
  // 32bit little endian
  } else if ( LEbyte.i[0]==84 && LEbyte.i[1]==dcdkey) {
    isBigEndian = false;
    is64bit = false;
  } else {
    // Flip bytes to convert to big endian
    BEbyte.c[0] = LEbyte.c[3];
    BEbyte.c[1] = LEbyte.c[2];
    BEbyte.c[2] = LEbyte.c[1];
    BEbyte.c[3] = LEbyte.c[0];
    BEbyte.c[4] = LEbyte.c[7];
    BEbyte.c[5] = LEbyte.c[6];
    BEbyte.c[6] = LEbyte.c[5];
    BEbyte.c[7] = LEbyte.c[4];
    // If the sum of the first 2 4-byte blocks is 84 this is a 
    // 64bit big endian
    if ( (BEbyte.i[0] + BEbyte.i[1]) == 84 ) {
      isBigEndian = true;
      is64bit = true;
    // If first 4-byte block is 84 and second is CORD this is
    // 32bit big endian
    } else if ( BEbyte.i[0]==84 && LEbyte.i[1]==dcdkey) {
      isBigEndian = true;
      is64bit = false;
    // Otherwise who know what the heck this is
    } else {
      mprinterr("Error: Unrecognized DCD header [%s].\n",tfile->filename);
      return -1;
    }
  }
  // If 64 bit check next 4 bytes for the dcd key
  if (is64bit) {
    tfile->IO->Read(LEbyte.c, sizeof(unsigned char), 4);
    if ( LEbyte.i[0] != dcdkey ) {
      mprinterr("Error: DCD key not found in 64 bit Charmm DCD file.\n");
      return -1;
    }
  }

  mprintf("\tDCD header:");
  if (isBigEndian)
    mprintf(" Big Endian");
  else
    mprintf(" Little Endian");
  if (is64bit) {
    mprintf(" 64 bit.\n");
    readSize = 8;
  } else {
    mprintf(" 32 bit.\n");
    readSize = 4;
  }

  // Buffer the rest of the header
  if (tfile->IO->Read(buffer, sizeof(unsigned char), 80) < 1) {
    mprinterr("Error: Could not buffer DCD header.\n");
    return -1;
  }

  // load ICNTRL variables and print for debugging
  for (int i=0; i < 80; i+=4) {
    mprintf("\ticntrl[%i]= %i\n",i,*( (int*) (buffer + i))  );
  }

  // Make sure this is Charmm format; last integer in the header should not
  // be zero.
  if ( *( (int*) (buffer + 76) ) != 0 ) {
    mprintf("\tCharmm DCD\n");
    // Check for Charmm-specific flags
    if ( *( (int*) (buffer + 40) ) != 0 ) dcdExtraBlock = true;
    if ( *( (int*) (buffer + 44) ) != 0 ) dcd4D = true;
  } else {
    mprinterr("\tNon-charmm DCD - currently unsupported.\n");
    return -1;
  }

  // Number of sets
  Frames = *( (int*) (buffer     ) );
  // Starting timestep
  istart = *( (int*) (buffer + 4 ) );
  // Number of steps between frames
  nsavc  = *( (int*) (buffer + 8 ) );
  // Number of fixed atoms
  namnf  = *( (int*) (buffer + 32) );
  // Timestep - float
  timestep = *( (float*) (buffer + 36) );
  mprintf("\tTimestep is %f\n",timestep);

  // Read end size of first block
  LEbyte.i[1] = 0;
  if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
    mprinterr("Error: Could not read second 84 from DCD.\n");
    return -1;
  } 
  if ( (LEbyte.i[0] + LEbyte.i[1]) != 84 ) {
    mprinterr("Error: Bad size at end of first block (%i) of DCD.\n",LEbyte.i[0] + LEbyte.i[1]);
    return -1;
  }

  // Read title block size
  LEbyte.i[1] = 0;
  if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
    mprinterr("Error: Could not read size of title block from DCD.\n");
    return -1;
  }  
  // Read titles    
  mprintf("\tTitle block size %i\n",LEbyte.i[0] + LEbyte.i[1]);
  if ( (((LEbyte.i[0] + LEbyte.i[1]) - 4) % 80) == 0 ) {
    // Read ntitle
    if (tfile->IO->Read(&ntitle,sizeof(int),1) < 1) {
      mprintf("Error: DCD Reading ntitle.\n");
      return -1;
    }
    mprintf("\tNtitle %i\n",ntitle);
    for (int i=0; i < ntitle; i++) {
      tfile->IO->Read(dcdtitle,sizeof(char),80);
      mprintf("\tTitle%i: [%s]\n",i+1,dcdtitle);
    }
  }
  // Read title end block size
  LEbyte.i[1] = 0;
  if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
    mprinterr("Error: Could not read size of title end block from DCD.\n");
    return -1;
  }

  // Read in next block size, should be 4
  LEbyte.i[1] = 0;
  if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
    mprinterr("Error: Could not read size of 4 from DCD after title.\n");
    return -1;
  }
  if ( (LEbyte.i[0] + LEbyte.i[1]) != 4 ) {
    mprinterr("Error: Expected to read block size of 4 after title block.\n");
    return -1;
  }
  // Read in number of atoms
  if (tfile->IO->Read(&dcdatom,sizeof(int),1) < 1) {
    mprintf("Error: DCD reading natom.\n");
    return -1;
  }
  mprintf("\tNatom %i\n",dcdatom);

  // Read in next block size, should also be 4
  LEbyte.i[1] = 0;
  if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
    mprinterr("Error: Could not read size of 4 from DCD after natom.\n");
    return -1;
  }  
  if ( (LEbyte.i[0] + LEbyte.i[1]) != 4 ) {
    mprinterr("Error: Expected to read block size of 4 after natom block.\n");
    return -1;
  }

  // If number of fixed atoms not 0, need to read list of free atoms.
  if (namnf!=0) {
    // Set nfreat, natom - namnf
    nfreat = dcdatom - namnf;
    mprintf("\tNfreat %i\n",nfreat);
    // Allocate space for nfreat atom indices
    freeat = new int[ nfreat ];
    // Read index array size
    LEbyte.i[1] = 0;
    if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
      mprinterr("Error: Could not read index array size from DCD.\n");
      return -1;
    }
    if ( (LEbyte.i[0] + LEbyte.i[1]) != (nfreat * 4) ) {
      mprinterr("Error: DCD: Expected index array size %i, got %i\n",(nfreat * 4),
                (LEbyte.i[0] + LEbyte.i[1]));
      return -1;
    }
    // Read index array
    if (tfile->IO->Read( freeat, sizeof(int), nfreat) < 1) {
      mprinterr("Error reading DCD index array.\n");
      return -1;
    }
    // Read end index array size
    LEbyte.i[1] = 0;
    if (tfile->IO->Read(LEbyte.c, sizeof(unsigned char), readSize) < 1) {
      mprinterr("Error: Could not read end index array size from DCD.\n");
      return -1;
    }
    if ( (LEbyte.i[0] + LEbyte.i[1]) != (nfreat * 4) ) {
      mprinterr("Error: DCD: Expected end index array size %i, got %i\n",(nfreat * 4),
              (LEbyte.i[0] + LEbyte.i[1]));
      return -1;
    }
  }


  closeTraj();

  return Frames;
}

/* CharmmDcd::readFrame()
 */
int CharmmDcd::readFrame(int set,double *X, double *V,double *box, double *T) {

  return 1;
}

/* CharmmDcd::processWriteArgs()
 */
int CharmmDcd::processWriteArgs(ArgList *argIn) {
  return 0;
}

/* CharmmDcd::setupWrite()
 * Set parm information required for write, and check write mode against
 * number of frames to be written.
 */
int CharmmDcd::setupWrite(AmberParm *trajParm) {
  return 1;
}

/* CharmmDcd::writeFrame()
 */
int CharmmDcd::writeFrame(int set, double *X, double *V,double *box, double T) {

  return 1;
}
 
/* CharmmDcd::info()
 */
void CharmmDcd::info() {
  mprintf("is a CHARMM DCD file");
}
