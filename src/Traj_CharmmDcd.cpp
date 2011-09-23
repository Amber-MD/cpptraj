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
  timestep=0;
}

// DESTRUCTOR
CharmmDcd::~CharmmDcd() {
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
  byte B;
  bool isBigEndian = false;
  bool is64bit = false;
  unsigned char buffer[80];

  if ( openTraj() ) return -1;

  // Read first 8 bytes - first number in dcd header should be 84
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
    if ( (BEbyte.i[0] + BEbyte.i[1]) == 84 ) {
      isBigEndian = true;
      is64bit = true;
    } else if ( BEbyte.i[0]==84 && LEbyte.i[1]==dcdkey) {
      isBigEndian = true;
      is64bit = false;
    } else {
      mprinterr("Error: Unrecognized DCD header [%s].\n",tfile->filename);
      return -1;
    }
  }
  // If 64 bit check next 4 bytes for the dcd key
  if (is64bit) {
    tfile->IO->Read(B.c, sizeof(unsigned char), 4);
    if ( B.i != dcdkey ) {
      mprinterr("Error: DCD key not found in 64 bit Charmm DCD file.\n");
      return -1;
    }
  }

  mprintf("\tDCD header:");
  if (isBigEndian)
    mprintf(" Big Endian");
  else
    mprintf(" Little Endian");
  if (is64bit)
    mprintf(" 64 bit.\n");
  else
    mprintf(" 32 bit.\n");

  // Buffer the rest of the header
  if (tfile->IO->Read(buffer, sizeof(unsigned char), 80) < 1) {
    mprinterr("Error: Could not buffer DCD header.\n");
    return -1;
  }

  //binaryByteOrder(84,&byteorder);

  //mprintf("CHARMMDCD: Magic %c %c %c %c\n",magic.c[0],magic.c[1],magic.c[2],magic.c[3]);

  // load ICNTRL variables and print for debugging
  for (int i=0; i < 80; i+=4) {
    B.c[0] = buffer[i  ];
    B.c[1] = buffer[i+1];
    B.c[2] = buffer[i+2];
    B.c[3] = buffer[i+3];
    mprintf("\ticntrl[%i]= %i\n",i,B.i);
  }

  // Make sure this is Charmm format; last integer in the header should not
  // be zero.
  if ( *( (int*) (buffer + 76) ) != 0 ) {
    mprintf("\tCharmm DCD\n");
    // Check for Charmm-specific flags
    //if ( *( (int*) (buffer + 40) ) != 0 ) dcdExtraBlock = true;
    //if ( *( (int*) (buffer + 44) ) != 0 ) dcd4D = true;
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
