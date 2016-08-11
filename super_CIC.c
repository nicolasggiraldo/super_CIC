#include          <stdio.h>
#include         <stdlib.h>
#include           <math.h>
#include         <string.h>
#include <gsl/gsl_spline.h>
#include            <mpi.h>

#include "constansts_and_structures.h"
#include                 "functions.h"

int main(int argc, char *argv[])
{
  int i,j,k;
  int ii, jj, kk; // Aditional counters for the scheme assignation calculation
  unsigned long ul_i;
  int index_start=0, index_end=0; // Index to start and finish according to the scheme used
  double *denCon;
  double *denCon_recv;
  long int index_cell;
  //long int indexaux;
  unsigned long Npart_snap;
  unsigned int dummy;
  double (*W)(double,
	      double,
	      double,
	      double)=NULL; // Memory addres to the window function
  struct gadget_head header;
  int err;
  int s;
  char buf[300];
  char label[4];
  float pos[3]; //float vel[3];
  double xc, yc, zc;
  unsigned long blksize_pos;
  //unsigned long blksize_vel;
  //unsigned long blksize_ids;

  FILE *fp_pos  =NULL;
  //FILE *fp_vel  =NULL;
  //FILE *fp_ids  =NULL;
  FILE *fp_head = NULL;
  FILE *outfile = NULL;


  
  // MPI variables
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; // Value of the processor id number
  int size; // Number of processors to use
  MPI_Status status;
  int chunk;
  
  
  
  /////////////////////
  //* MPI BEGGINING *//
  /////////////////////
  MPI_Init(&argc, &argv);
  
  // Number of processors
  MPI_Comm_size(comm, &size);
  // Task numbering
  MPI_Comm_rank(comm, &rank);
  
  
  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(argc<2)
    {
      if(rank==0)
	{
	  printf("\n***********************************");
	  printf("***********************************\n");
	  printf("%s: You must specify the name of the parameter file\n",argv[0]);
	  printf("For example: %s path_of_file/parameter_file\n",argv[0]);
	  printf("***********************************");
	  printf("***********************************\n\n");
	}
      exit(0);
    }

  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1],rank) )
    {
    case -1 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad path to the parameter file.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -2 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad settings in the parameter file.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }

  // Reading header
  sprintf(buf,"%s.%d",GV.SNAP_BASE,0);
  if( read_head(buf, rank)==-1 )
    {
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad path to header.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
  
  
  
  ///////////////////////////////////////////
  //* DEFINING THE WINDOW FUNCTION TO USE *//
  ///////////////////////////////////////////
  if(      strcmp(GV.SCHEME, "NGP") == 0 )
    {
      // NGP
      W = W_NGP;
      index_start = 0;
      index_end   = 0;
    }
  else if( strcmp(GV.SCHEME, "CIC") == 0 )
    {
      // CIC
      W = W_CIC;
      index_start = -1;
      index_end   =  1;
    }
  else if( strcmp(GV.SCHEME, "TSC") == 0 )
    {
      // TSC
      W = W_TSC;
      index_start = -2;
      index_end   =  2;
    }
  else if( strcmp(GV.SCHEME, "D20") == 0 )
    {
      // D20
      len_array_D20 = 610;
      x_D20 = (double *) calloc( len_array_D20, sizeof(double));
      y_D20 = (double *) calloc( len_array_D20, sizeof(double));
      
      fin_D20 = fopen("./D20.txt", "r");
      
      for(i=0; i<len_array_D20; i++)
	err=fscanf(fin_D20, "%lf %lf", &x_D20[i], &y_D20[i]);
      
      fclose(fin_D20);
      
      // GSL interpolation allocation
      acc    = gsl_interp_accel_alloc(); // accelerator
      spline = gsl_spline_alloc(gsl_interp_cspline, len_array_D20); // spline
      
      // GSL init
      gsl_spline_init(spline, x_D20, y_D20, len_array_D20);
      
      W = W_D20;
      index_start = 0;
      index_end   = 10; //index_end   = 8;
    }
  
  
  
  //////////////////////////////
  //* CELL STRUCT ALLOCATION *//
  //////////////////////////////
  /* Array of structure Cell, size NGRID^3 */
  //cells = (struct Cell *) calloc( GV.NGRID3, sizeof( struct Cell) );
  denCon = (double *) calloc( GV.NGRID3, sizeof(double) );
  
  
  
  //////////////////////////////
  //* FROM PARTICLES TO GRID *//
  //////////////////////////////
  if(rank==0)
    {
      printf("\nLocating particles into the grid\n");
      printf("-----------------------------------------------\n");
      fflush(stdout);
    }

  //domain decomposition 
  //npart[6]
  for(s=rank; s<GV.NSNAPS; s+=size)
    {
      sprintf(buf,"%s.%d",GV.SNAP_BASE,s);

      // Reading the positions, velocities and ID's of the snapshot
      if( (fp_head=fopen(buf,"r"))==NULL )
	{
	  printf("rank:%d. Cannot open header %s",rank,buf);
	  exit(0);
	}
      fp_pos = fopen(buf,"r");
      //fp_vel = fopen(buf,"r");
      //fp_ids = fopen(buf,"r");
      
      /* reading the header for this sub snap */
      if(GV.GADGET_VERSION==2)
	{
	  err=fread(&dummy, sizeof(dummy), 1, fp_head);
	  err=fread(&label, sizeof(char),  4, fp_head);
	  err=fread(&dummy, sizeof(dummy), 1, fp_head);
	  err=fread(&dummy, sizeof(dummy), 1, fp_head);
	}
      
      err=fread(&dummy,  sizeof(dummy),  1, fp_head);
      err=fread(&header, sizeof(header), 1, fp_head);
      err=fread(&dummy,  sizeof(dummy),  1, fp_head);
      
      fclose(fp_head);

      //Npart_snap = header.npart[GAS]+header.npart[HALO]+
      //header.npart[DISK]+header.npart[BULGE]+header.npart[STARS]+header.npart[BNDRY];

      Npart_snap = header.npart[HALO];
      
      printf("rank:%d   snap:%d\n * N[GAS]  :%u\n * N[HALO] :%u\n * N[DISK] :%u\n * N[BULGE]:%u\n * N[STARS]:%u\n * N[BNDRY]:%u\n\n",
	     rank,s,
	     header.npart[GAS],header.npart[HALO],header.npart[DISK],
	     header.npart[BULGE],header.npart[STARS],header.npart[BNDRY]);
      
      //  Begining with the groups
      if(GV.GADGET_VERSION==2)
	{
	  blksize_pos =  9*sizeof(int) + 2*4*sizeof(char) + sizeof(header);
	  //blksize_vel = 14*sizeof(int) + 3*4*sizeof(char) + sizeof(header) +   3*Npart_snap*sizeof(float);
	  //blksize_ids = 19*sizeof(int) + 4*4*sizeof(char) + sizeof(header) + 2*3*Npart_snap*sizeof(float);
	}
      else
	{
	  blksize_pos = 3*sizeof(int) + sizeof(header);
	  //blksize_vel = 5*sizeof(int) + sizeof(header) +   3*Npart_snap*sizeof(float);
	  //blksize_ids = 7*sizeof(int) + sizeof(header) + 2*3*Npart_snap*sizeof(float);
	}
      /* In case there is GAS jump that extra space :P */
      blksize_pos += header.npart[GAS]*3*sizeof(float);
      
      fseek(fp_pos, blksize_pos, SEEK_SET);
      //fseek(fp_vel, blksize_vel, SEEK_SET);
      //fseek(fp_ids, blksize_ids, SEEK_SET);

      for(ul_i=0; ul_i<Npart_snap ;ul_i++)
	{
	  err=fread(&pos[0], sizeof(float),  3, fp_pos);
	  //err=fread(&vel[0], sizeof(float),  3, fp_vel);
	  //err=fread(&id,     sizeof(idtype), 1, fp_ids);
	  
	  // index in the x, y and z-axis
	  i = (int) floor( (pos[X] / GV.L) * GV.NGRID );
	  j = (int) floor( (pos[Y] / GV.L) * GV.NGRID );
	  k = (int) floor( (pos[Z] / GV.L) * GV.NGRID );
	  
	  // index in C-order
	  //index_cell = INDEX(i,j,k);
	  
	  /* Calculating scheme assignation */
	  for( ii=index_start; ii<=index_end; ii++ )
	    {
	      for( jj=index_start; jj<=index_end; jj++ )
		{
		  for( kk=index_start; kk<=index_end; kk++ )
		    {
		      index_cell = INDEX(mod(i+ii, GV.NGRID),
					 mod(j+jj, GV.NGRID),
					 mod(k+kk, GV.NGRID));
		      xc = GV.H * (0.5 + i+ii);
		      yc = GV.H * (0.5 + j+jj);
		      zc = GV.H * (0.5 + k+kk);
		      denCon[index_cell] += header.mass[HALO] * W(xc-pos[X], yc-pos[Y], zc-pos[Z], GV.H);
		    }
		}
	    }
	}
      fclose(fp_pos);
    }
    
  MPI_Barrier(comm);



  ///////////////////////////////////////////////
  //* SENDING DENCON OF OTHER TASKS TO RANK 0 *//
  ///////////////////////////////////////////////
  
  chunk = (GV.NGRID3)/SIZEDOUBLE;

  if(rank==0)
    {
      denCon_recv = (double *) calloc( GV.NGRID3, sizeof(double) );

      for(i=1;i<size;i++)//FOR IN RANKS
	{
	  for(j=0; j<SIZEDOUBLE; j++)//FOR IN CHUNKS
	    {
	      MPI_Recv(denCon_recv+(j*chunk), chunk, MPI_DOUBLE, i, j, comm, &status);
	    }
	  
	  for(index_cell=0; index_cell<GV.NGRID3; index_cell++)
	    {
	      denCon[index_cell] += denCon_recv[index_cell];
	    }
	  
	  printf("Recived density constrast from rank:%d to rank:%d\n",i,rank);
	  fflush(stdout);
	  
	}
      free(denCon_recv);
    }
  else
    {
      for(j=0; j<SIZEDOUBLE; j++)//FOR IN CHUNKS
	    {
	      MPI_Send(denCon+(j*chunk), chunk, MPI_DOUBLE, 0, j, comm);
	    }
    }
  
  /////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////
  //* TERMINATING THE CALCULATION AND SAVING IN FILE *//
  //////////////////////////////////////////////////////
  if(rank==0)
    { 
      char *last = strrchr(GV.SNAP_BASE,'/');
      
      if(last==NULL)
	last=GV.SNAP_BASE;
      else
	last=last+1;

      if( strcmp(GV.SCHEME, "NGP") == 0 )
	{
	  // NGP
	  sprintf(buf, "%s/%s_NGRID_%d_NGP.dens", GV.OUTPUT_DIR, last, GV.NGRID);
	}
      else if( strcmp(GV.SCHEME, "CIC") == 0 )
	{
	  // CIC
	  sprintf (buf, "%s/%s_NGRID_%d_CIC.dens", GV.OUTPUT_DIR, last, GV.NGRID);
	}
      else if( strcmp(GV.SCHEME, "TSC") == 0 )
	{
	  // TSC
	  sprintf (buf, "%s/%s_NGRID_%d_TSC.dens", GV.OUTPUT_DIR, last, GV.NGRID);
	}
      else if( strcmp(GV.SCHEME, "D20") == 0 )
	{
	  // D20
	  sprintf (buf, "%s/%s_NGRID_%d_D20.dens", GV.OUTPUT_DIR, last, GV.NGRID);
	}
  
      // Opening file for output of the cell
      outfile = fopen(buf, "wb");

      printf("Saving data in %s\n", buf);
      printf("-----------------------------------------------\n");
      fflush(stdout);
      
      /* Saving cosmological parameters of the simulation */
      fwrite( &(Gheader.Omega0),      sizeof(double), 1, outfile);
      fwrite( &(Gheader.OmegaLambda), sizeof(double), 1, outfile);
      fwrite( &(Gheader.redshift),    sizeof(double), 1, outfile);
      fwrite( &(Gheader.HubbleParam), sizeof(double), 1, outfile);
      
      /* Saving simulation parameters */
      fwrite(          &(GV.NGRID),           sizeof(int), 1, outfile);
      fwrite( &(GV.GADGET_VERSION),           sizeof(int), 1, outfile);
      fwrite(              &(GV.L),        sizeof(double), 1, outfile);
      fwrite(         &(GV.NP_TOT), sizeof(unsigned long), 1, outfile);
      fwrite(     &(GV.TOTAL_MASS),        sizeof(double), 1, outfile);
      fwrite(       &(GV.RHO_MEAN),        sizeof(double), 1, outfile);
      fwrite(       &(GV.VOL_CELL),        sizeof(double), 1, outfile);
      fwrite(              &(GV.H),        sizeof(double), 1, outfile);
      fwrite(      &(GV.SCHEME[0]),          sizeof(char), 1, outfile);
      fwrite(      &(GV.SCHEME[1]),          sizeof(char), 1, outfile);
      fwrite(      &(GV.SCHEME[2]),          sizeof(char), 1, outfile);
      
      // Setting total mass to zero
      double totalMass = 0.0;

      /* Saving cell data */
      for(i=0; i<GV.NGRID; i++)
	{
	  for(j=0; j<GV.NGRID; j++)
	    {
	      for(k=0; k<GV.NGRID; k++)
		{	  
		  index_cell = INDEX(i,j,k);
		  
		  /* Calculating the final density in the cell. This 
		     is made as a verification of the scheme, the 
		     mass must conservate with any scheme used. */
		  totalMass += denCon[index_cell];
		  
		  /* Calculating the density in the cell */
		  denCon[index_cell] /= GV.VOL_CELL;
		  
		  /* Calculating the final density contrast in the cell */
		  denCon[index_cell] = (denCon[index_cell] / GV.RHO_MEAN) - 1.0;
		  
		  // Writing density contrast
		  fwrite( &(denCon[index_cell]), sizeof(double), 1, outfile );
		}
	    }
	}
      fclose(outfile);
      
      /* Validation of the scheme */
      if( strcmp(GV.SCHEME, "NGP") == 0 )
	{
	  // NGP
	  printf("Mass NGP        = %lf\n", totalMass);
	}
      else if( strcmp(GV.SCHEME, "CIC") == 0 )
	{
	  // CIC
	  printf("Mass CIC        = %lf\n", totalMass);
	}
      else if( strcmp(GV.SCHEME, "TSC") == 0 )
	{
	  // TSC
	  printf("Mass TSC        = %lf\n", totalMass);
	}
      else if( strcmp(GV.SCHEME, "D20") == 0 )
	{
	  // D20
	  printf("Mass D20        = %lf\n", totalMass);
	}
      
      printf("Mass Simulation = %lf\n", GV.TOTAL_MASS);
      printf("%% Difference    = %.4e%%\n",
	     (100.0 * (totalMass-GV.TOTAL_MASS)) / GV.TOTAL_MASS);
    }
  

  
  ///////////////////
  //* FREE MEMORY *//
  ///////////////////
  if( strcmp(GV.SCHEME, "D20") == 0 )
    {
      // D20 
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      free(x_D20);
      free(y_D20);
    }
  free(denCon);

  if(err){}

  MPI_Finalize();
  return 0;
}
