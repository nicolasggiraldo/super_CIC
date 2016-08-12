////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTURES FOR THE PROGRAM //
////////////////////////////////////////////////////



/////////////////////////////
// PREPROCESSOR DIRECTIVES //
/////////////////////////////
#define X 0
#define Y 1
#define Z 2

#define GAS   0
#define HALO  1
#define DISK  2
#define BULGE 3
#define STARS 4
#define BNDRY 5

#define LENCHAR 300

#define SIZEDOUBLE sizeof(double)

#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define INDEX(i,j,k) (1L*k) + GV.NGRID * ( (1L*j) + GV.NGRID * (1L*i) )



////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES FOR INTERPOLATION IN DAUBECHIES-TYPES SCHEMES //
////////////////////////////////////////////////////////////////////
gsl_interp_accel *acc=NULL;
gsl_spline       *spline=NULL;
int               len_array_D20=0;
double           *x_D20=NULL;
double           *y_D20=NULL;
FILE             *fin_D20=NULL;



/////////////////////////////////////////////
// STRUCTURES OF THE PROGRAM AND NEW TYPES //
/////////////////////////////////////////////

/* Global variales */
struct globalVariables
{
  /* CELL PROPERTIES */
  int           NGRID;          // Number of cell in each axis.
  long int      NGRID3;         // Total number of cells (NGRID3 = NGRID^3)
  char         *SCHEME;         // Scheme used for grid assignation
  double        H;              // Size of the cell
  double        VOL_CELL;       // Volume of each cell

  /* SNAPSHOT INFORMATION */
  char         *SNAP_BASE;      // SNAP_BASE name (including path) of the GADGET snapshot
  int           NSNAPS;         // Numeber of subsnapshots of the simulation
  int           GADGET_VERSION; // GADGET version of the snapshot

  /* SIMULATION PROPERTIES */
  double        L;              // Lenght of the simulation in Mpc
  unsigned long NP_TOT;         // Total number of particles in the simulation
  double        TOTAL_MASS;     // Total mass of all particles in the simulation
  double        RHO_MEAN;       // Mean density of ALL the simulation
  
  
  /* COSMOLOGICAL PARAMETERS OF THE SIMULATION */
  //double        OMEGA_M0;       //Omega matter at present time
  //double        OMEGA_L0;       //Omega Lambda at present time
  //double        ZRS;            //Redshift of the simulation
  //double        HUBBLEPARAM;    //Hubble parameter of the simulation

  /* INFORMATION OF THE OUTPUT DIRECTORY */
  char         *OUTPUT_DIR;     // OUTPUT DIRECTORY of the density contrast
}GV;

/* Gadget header */
struct gadget_head
{
  unsigned int npart[6];
  double       mass[6];
  double       time;
  double       redshift;
  int          flag_sfr;
  int          flag_feedback;
  unsigned int npartTotal[6];
  int          flag_cooling;
  int          num_files;
  double       BoxSize;
  double       Omega0;
  double       OmegaLambda;
  double       HubbleParam;
  int          flag_age;
  int          flag_metals;
  unsigned int nallHW[6];
  char         fill[256-
		    6*4-
		    6*8-
		    2*8-
		    2*4-
		    6*4-
		    2*4-
		    4*8-
		    2*4-
		    6*4]; // Fills to 256 Bytes
}Gheader;
