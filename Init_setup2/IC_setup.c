#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "gsl_randist.h"
//#include "allvars.h"
//#include "proto.h"

size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
int FatalError(int errnum);
void read_parameter_file(char *fname);

struct io_header
{
  int npart[6];//number of particles of each type in current file
  double mass[6];//mass of each particle type
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];//total number of particles of each type in the WHOLE simulation (i.e. in case each snapshot has multiple files)
  int flag_cooling;
  int num_files;//number of files in each snapshot
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char fill[60];	/* fills header structure to 256 Bytes */
} header;

  #define MAXLEN_FILENAME 100
  #define MAXLEN_OUTPUTLIST 500

struct global_data_all_processes
{
  long long TotNumPart;		/*!< total particle numbers (global value) */
  long long TotN_gas;		/*!< total gas particle number (global value) */

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */

  double BoxSize;               /*!< Boxsize in case periodic boundary conditions are used */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;/*!< maximum number of files that may be written simultaneously when
                                     writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density computation */
  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH hydrodynamical force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain decomposition */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				     NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				     the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				     the maximum(!) number of particles.  Note: A typical local tree for N
				     particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  double DesNumNgb;             /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */


  /* some force counters  */

  long long TotNumOfForces;	             /*!< counts total number of force computations  */
  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */


  /* system of units  */

  double G;                        /*!< Gravity-constant in internal units */
  double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h */
  double UnitMass_in_g;            /*!< factor to convert internal mass unit to grams/h */
  double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to cm/sec */
  double UnitLength_in_cm;         /*!< factor to convert internal length unit to cm/h */
  double UnitPressure_in_cgs;      /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
  double UnitDensity_in_cgs;       /*!< factor to convert internal length unit to g/cm^3*h^2 */
  double UnitCoolingRate_in_cgs;   /*!< factor to convert internal cooling rate to cgs units */
  double UnitEnergy_in_cgs;        /*!< factor to convert internal energy to cgs units */
  double UnitTime_in_Megayears;    /*!< factor to convert internal time to megayears/h */
  double GravityConstantInternal;  /*!< If set to zero in the parameterfile, the internal value of the
                                        gravitational constant is set to the Newtonian value based on the system of
                                        units specified. Otherwise the value provided is taken as internal gravity constant G. */


  /* Cosmological parameters */

  double Hubble;       /*!< Hubble-constant in internal units */
  double Omega0;       /*!< matter density in units of the critical density (at z=0)*/
  double OmegaLambda;  /*!< vaccum energy density relative to crictical density (at z=0) */
  double OmegaBaryon;  /*!< baryon density in units of the critical density (at z=0)*/
  double HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute physical values for cooling physics */
  

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified file */


  /* Parameters determining output frequency */

  int SnapshotFileCount;        /*!< number of snapshot that is written next */
  double TimeBetSnapshot;       /*!< simulation time interval between snapshot files */
  double TimeOfFirstSnapshot;   /*!< simulation time of first snapshot files */
  double CpuTimeBetRestartFile; /*!< cpu-time between regularly generated restart files */
  double TimeLastRestartFile;   /*!< cpu-time when last restart-file was written */
  double TimeBetStatistics;     /*!< simulation time interval between computations of energy statistics */
  double TimeLastStatistics;    /*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this point */


  /* Current time of the simulation, global step, and end of simulation */

  double Time;                  /*!< current time of the simulation */
  double TimeBegin;             /*!< time of initial conditions of the simulation */
  double TimeStep;              /*!< difference between current times of previous and current timestep */
  double TimeMax;	        /*!< marks the point of time until the simulation is to be evolved */


  /* variables for organizing discrete timeline */

  double Timebase_interval;     /*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;               /*!< current time on integer timeline */ 
  int Ti_nextoutput;            /*!< next output time on integer timeline */
#ifdef FLEXSTEPS
  int PresentMinStep;           /*!< If FLEXSTEPS is used, particle timesteps are chosen as multiples of the present minimum timestep. */
  int PresentMaxStep;		/*!< If FLEXSTEPS is used, this is the maximum timestep in timeline units, rounded down to the next power 2 division */
#endif
#ifdef PMGRID
  int PM_Ti_endstep;            /*!< begin of present long-range timestep */
  int PM_Ti_begstep;            /*!< end of present long-range timestep */
#endif


  /* Placement of PM grids */

#ifdef PMGRID
  double Asmth[2];              /*!< Gives the scale of the long-range/short-range split (in mesh-cells), both for the coarse and the high-res mesh */
  double Rcut[2];               /*!< Gives the maximum radius for which the short-range force is evaluated with the tree (in mesh-cells), both for the coarse and the high-res mesh */
  double Corner[2][3];          /*!< lower left corner of coarse and high-res PM-mesh */
  double UpperCorner[2][3];     /*!< upper right corner of coarse and high-res PM-mesh */
  double Xmintot[2][3];         /*!< minimum particle coordinates both for coarse and high-res PM-mesh */
  double Xmaxtot[2][3];         /*!< maximum particle coordinates both for coarse and high-res PM-mesh */
  double TotalMeshSize[2];      /*!< total extension of coarse and high-res PM-mesh */
#endif


  /* Variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;          /*!< CPU time limit as defined in parameterfile */
  double CPU_TreeConstruction;  /*!< time spent for constructing the gravitational tree */
  double CPU_TreeWalk;          /*!< actual time spent for pure tree-walks */
  double CPU_Gravity;           /*!< cumulative time used for gravity computation (tree-algorithm only) */
  double CPU_Potential;         /*!< time used for computing gravitational potentials */
  double CPU_Domain;            /*!< cumulative time spent for domain decomposition */
  double CPU_Snapshot;          /*!< time used for writing snapshot files */
  double CPU_Total;             /*!< cumulative time spent for domain decomposition */
  double CPU_CommSum;           /*!< accumulated time used for communication, and for collecting partial results, in tree-gravity */
  double CPU_Imbalance;         /*!< cumulative time lost accross all processors as work-load imbalance in gravitational tree */
  double CPU_HydCompWalk;       /*!< time used for actual SPH computations, including neighbour search */
  double CPU_HydCommSumm;       /*!< cumulative time used for communication in SPH, and for collecting partial results */
  double CPU_HydImbalance;      /*!< cumulative time lost due to work-load imbalance in SPH */
  double CPU_Hydro;             /*!< cumulative time spent for SPH related computations */
  double CPU_EnsureNgb;         /*!< time needed to iterate on correct neighbour numbers */
  double CPU_Predict;           /*!< cumulative time to drift the system forward in time, including dynamic tree updates */
  double CPU_TimeLine;          /*!< time used for determining new timesteps, and for organizing the timestepping, including kicks of active particles */
  double CPU_PM;                /*!< time used for long-range gravitational force */
  double CPU_Peano;             /*!< time required to establish Peano-Hilbert order */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                     timestep is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                     timestep determined by the timestep criteria falls below this limit. */ 
  double MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                     in comoving coordinates.  To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such that the rms displacement
                                     is a fraction of the mean particle separation (determined from the
                                     particle mass and the cosmological parameters). This parameter specifies
                                     this fraction. */

  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */

  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


  /* Gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length).
   * Five groups of particles are supported 0="gas", 1="halo", 2="disk", 3="bulge", 4="stars", 5="bndry"
   */

  double MinGasHsmlFractional;  /*!< minimum allowed SPH smoothing length in units of SPH gravitational softening length */
  double MinGasHsml;            /*!< minimum allowed SPH smoothing length */


  double SofteningGas;          /*!< comoving gravitational softening lengths for type 0 */ 
  double SofteningHalo;         /*!< comoving gravitational softening lengths for type 1 */ 
  double SofteningDisk;         /*!< comoving gravitational softening lengths for type 2 */ 
  double SofteningBulge;        /*!< comoving gravitational softening lengths for type 3 */ 
  double SofteningStars;        /*!< comoving gravitational softening lengths for type 4 */ 
  double SofteningBndry;        /*!< comoving gravitational softening lengths for type 5 */ 

  double SofteningGasMaxPhys;   /*!< maximum physical softening length for type 0 */ 
  double SofteningHaloMaxPhys;  /*!< maximum physical softening length for type 1 */ 
  double SofteningDiskMaxPhys;  /*!< maximum physical softening length for type 2 */ 
  double SofteningBulgeMaxPhys; /*!< maximum physical softening length for type 3 */ 
  double SofteningStarsMaxPhys; /*!< maximum physical softening length for type 4 */ 
  double SofteningBndryMaxPhys; /*!< maximum physical softening length for type 5 */ 

  double SofteningTable[6];     /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];     /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  double MassTable[6];          /*!< Table with particle masses for particle types with equal mass.
                                     If particle masses are all equal for one type, the corresponding entry in MassTable 
                                     is set to this value, allowing the size of the snapshot files to be reduced. */
  


  /* some filenames */
  char InitCondFile[MAXLEN_FILENAME];          /*!< filename of initial conditions */
  char OutputDir[MAXLEN_FILENAME];             /*!< output directory of the code */
  char SnapshotFileBase[MAXLEN_FILENAME];      /*!< basename to construct the names of snapshotf files */
  char EnergyFile[MAXLEN_FILENAME];            /*!< name of file with energy statistics */
  char CpuFile[MAXLEN_FILENAME];               /*!< name of file with cpu-time statistics */
  char InfoFile[MAXLEN_FILENAME];              /*!< name of log-file with a list of the timesteps taken */
  char TimingsFile[MAXLEN_FILENAME];           /*!< name of file with performance metrics of gravitational tree algorithm */
  char RestartFile[MAXLEN_FILENAME];           /*!< basename of restart-files */
  char ResubmitCommand[MAXLEN_FILENAME];       /*!< name of script-file that will be executed for automatic restart */
  char OutputListFilename[MAXLEN_FILENAME];    /*!< name of file with list of desired output times */

  double OutputListTimes[MAXLEN_OUTPUTLIST];   /*!< table with desired output times */
  int OutputListLength;                        /*!< number of output times stored in the table of desired output times */

} All;                                          /*!< a container variable for global variables that are equal on all processors */



/*int NumPart, Ngas;

double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam;
double PrimordialIndex;
double ShapeGamma;
double Box;
int Nmesh, Nsample;
char     GlassFile[500];
char     FileWithInputSpectrum[500];
int      TileFac;
int Seed;
int      SphereMode;
int  NumFilesWrittenInParallel;
char OutputDir[100], FileBase[100];
double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
int      WhichSpectrum;
double InputSpectrum_UnitLength_in_cm;
int    ReNormalizeInputSpectrum;
int      ThisTask, NTask;*/

struct particle_data
{
  float Pos[3];//comoving particle coords (default h^-1 kpc)
  float Vel[3];//particle velocities (default km/s)
  float Mass;//mass of individual particles if not constant
  int Type;//Type of particle e.g. gas

  float Rho, U, Temp, Ne;
  //Rho = comovingdensity of SPH particles [default 10^10 h^-1 Msol / (h^-1 kpc)^3]
  //U = internal energy per unit mass for gas particles 
} *P;

int *Id;//particle ID

double Time, Redshift;

int main(int argc, char *argv[])
{
#define BUFFER 10
  FILE *fd;
  int i, j, k, pc;
  char pfname[200];
  char fname[200];
  //char buf[300];
  float density_ratio, d1, d2, P1, P2, rho1, rho2, gamma, U1, U2;
  int N1, N2;
  size_t bytes;
  float *block;
  int *blockid;
  int blockmaxlen, maxidlen;
  int dummy;

  strcpy(pfname, argv[1]);
  strcpy(fname, argv[2]);
  printf("Parameter file = %s\n", pfname);

  /*printf("Enter parameter filename:\n");
  scanf("%s", &pfname);

  printf("Enter IC filename to write:\n");
  scanf("%s", &fname);*/

 read_parameter_file(pfname);

 printf("Boxsize = %f\n", All.BoxSize);

  for(i=0; i<6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.mass[i] = 0;
      header.npartTotalHighWord[i] = 0;
    }

  header.npart[0] = 200;//number of gas particles
  header.npartTotal[0] = header.npart[0];//assumes no other particle types
  header.mass[0] = 1;//mass of gas particles
  density_ratio = 8;//ratio of density in first half to second half of box
  N2=(header.npart[0]-2)/(density_ratio+1);
  N1=header.npart[0]-N2;
  printf("N_tot=%d N1=%d N2=%d\n", header.npart[0], N1, N2);
  d1=(All.BoxSize/2)/(float)(N1);
  d2=(All.BoxSize/2)/(float)(N2-1);
  printf("d1=%f d2=%f\n", d1, d2);
  rho1 = header.mass[0]/(pow(d1,3));//density in first half of box
  rho2 = header.mass[0]/(pow(d2,3));
  P1 = 30.0 * rho1; //P1 is the pressure in the first half of the box
  P2 = 0.14/0.125 * rho2;
  gamma = 1.4; //the adiabatic index
  U1 = (P1 / rho1) / (gamma - 1);
  U2 = (P2 / rho2) / (gamma - 1);

  header.time = 0;//doesnt actually matter what this is, it's ignored in the IC
  header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;
  header.num_files = All.NumFilesPerSnapshot;//number of files in each snapshot

  /*long seed;
seed = time (NULL) * getpid();
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);*/


  P = malloc(bytes = header.npart[0] * sizeof(struct particle_data));
  Id = malloc(header.npart[0] * sizeof(int));

  j=1;
  for(i=0; i<header.npart[0]; i++)
    {
      Id[i]=i;
      P[i].Mass = header.mass[0];
      printf("Mass %d = %f\t", i+1, P[i].Mass);
      P[i].Type = 0;
      P[i].Vel[0] = 0;
      P[i].Vel[1] = 0;
      P[i].Vel[2] = 0;
      P[i].Pos[1] = 0;
      P[i].Pos[2] = 0;
      if(i<=N1)
	{
	  P[i].Pos[0] = (float)(i) * d1;
	}
      if(i>N1)
	{
	  P[i].Pos[0] = P[N1].Pos[0] + (float)(j) * d2;
	  j++;
	}
      //set internal energies per unit mass
      if(i<N1)
	{
	  P[i].U = U1;
	}
      if(i>=N1)
	{
	  P[i].U = U2;
	}
      printf("Pos%d = %f\t U%d = %f\n", i+1, P[i].Pos[0], i+1, P[i].U);
    }

  /*Start writing to file fd*/
  
  fd = fopen(fname, "w");
  if (!fd)
    {
      printf("Unable to open file.\n");
      return 1;
    }

  /*Write the header*/
  dummy = sizeof(header);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);//this is the same as SKIP in Gadget - it writes a dummy block to the file
  my_fwrite(&header, sizeof(header), 1, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);//another block-size field after the data block

  /*Write particle position coordinates*/
  block = malloc(bytes = BUFFER * 1024 * 1024);
  blockmaxlen = bytes / (3 * sizeof(float));
  dummy = sizeof(float) * 3 * header.npart[0];
  my_fwrite(&dummy, sizeof(dummy), 1, fd);//SKIP
  for(i=0, pc=0; i < header.npart[0]; i++)
    {
      for(k = 0; k< 3; k++)
	{
	  block[3 * pc + k] = P[i].Pos[k];//block contains x,y,z in sequence for all particles
	}
      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
  my_fwrite(block, sizeof(float), 3 * pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  /* Write particle velocities*/
  dummy = sizeof(float) * 3 * header.npart[0];
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i =0, pc =0; i < header.npart[0]; i++)
    {
      for(k = 0; k < 3; k++)
	{

	}
      pc++;
      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    {
      my_fwrite(block, sizeof(float), 3 * pc, fd);
    }
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  /* Write Particle IDs*/
  blockid = (int *) block;
  maxidlen = bytes / (sizeof(int));
  dummy = sizeof(int) * header.npart[0];
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < header.npart[0]; i++)
    {
      blockid[pc] = Id[i];
      pc++;

      if(pc == maxidlen)
	{
	  my_fwrite(blockid, sizeof(int), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    {
      my_fwrite(blockid, sizeof(int), pc, fd);
    }
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  /*Write internal energies*/
  dummy = sizeof(float) * header.npart[0];
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < header.npart[0]; i++)
    {
      block[pc] = P[i].U;
      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    {
      my_fwrite(block, sizeof(float), pc, fd);
    }
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  free(block);

  fclose(fd);

  return 0;
}

/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) has occured.\n");
      fflush(stdout);
      FatalError(777);
    }
  return nwritten;
}

int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  //MPI_Abort(MPI_COMM_WORLD, errnum);
  exit(0);
}

/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;


      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeAllocFactor");
      addr[nt] = &All.TreeAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = DOUBLE;

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	      system(buf3);
	    }
	}
      else
	{
	  printf("\nParameter file %s not found.\n\n", fname);
	  errorFlag = 2;
	}

      if(errorFlag != 2)
	for(i = 0; i < nt; i++)
	  {
	    if(*tag[i])
	      {
		printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
		errorFlag = 1;
	      }
	  }

      /*if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
      All.OutputListLength = 0;*/
    

//MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      //MPI_Finalize();
      exit(0);
    }

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}
