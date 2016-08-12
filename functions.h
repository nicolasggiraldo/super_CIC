/*
 * Function:  read_parameters
 * --------------------
 * Reads the parameter file in which are the main parameters 
 * necessary to run the code.
 *
 * The information loaded in the same order are: 
 * SNAP_BASE:      Snapbase name of the snapshot.
 * NSNAPS          Number of files related with the snapshot.
 * GADGET_VERSION  Version of gadget used.
 * NPARTICLES      Number of particles to use. The number of 
 *                 particles is writed in Npar^3.
 * NGRID           Number of grids in each side.
 * SCHEME          Mass assignment to use, valid params are
 *                 NGP, CIC, TSC, D20.
 * OUTPUT_DIR      Output directory to put the density contrast 
 *                 output.
 * 
 *  param_file_name: String with the name of the parameter file.
 *  rank:            Rank of the process.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the parameter file.
 *           -2 --> There is an error whith the settings of the 
 *                  parameter file.
 */
int read_parameters(char param_file_name[], int rank)
{
  FILE *cfg  = NULL; // Stream to the parameter (config) file
  int   len  = LENCHAR;  // Len of the read parameter
  char *buf  = NULL; // buffer variables used to read strings
  char *buf1 = NULL;
  char *buf2 = NULL; 
  char *dumb = NULL;
  
  char ngp[]    = "NGP";
  char cic[]    = "CIC";
  char tsc[]    = "TSC";
  char daub20[] = "D20";
  
  
  if( (cfg=fopen(param_file_name,"r"))==NULL )
    {
      printf("%s not found.\n", param_file_name);
      // Value -1 means there is an error loading the param file
      return -1;
    }
  
  buf  = (char *) malloc( len*sizeof(char) );
  buf1 = (char *) malloc( len*sizeof(char) );
  buf2 = (char *) malloc( len*sizeof(char) );
  
  /* Reading SNAP_BASE parameter */  
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'SNAP_BASE' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.SNAP_BASE = strdup(buf2);
      if(rank==0)
	printf("Snapshot base name: %s\n", GV.SNAP_BASE);
    }
  
  /* Reading NSNAPS parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'NSNAPS' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.NSNAPS = atoi(buf2);
      if(GV.NSNAPS>0)
	{
	  if(rank==0)
	    printf("Number of snapshots to read: %d\n", GV.NSNAPS);
	}
      else
	{
	  if(rank==0)
	    printf("Invalid 'NSNAPS' setting in configuration file.\n");
	  return -2;
	}
    }

  /* GADGET_VERSION parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'GADGET_VERSION' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.GADGET_VERSION=atoi(buf2);
      if(GV.GADGET_VERSION == 1 || GV.GADGET_VERSION == 2)
	{
	  if(rank==0)
	    printf("GADGET snapshot version: %d\n", GV.GADGET_VERSION);
	}
      else
	{
	  if(rank==0)
	    printf("Invalid 'GADGET_VERSION' setting in configuration file.\n");
	  return -2;
	}
    }
  
  /* NPARTICLES parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'NPARTICLES' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.NP_TOT=atoi(buf2);
      if(GV.NP_TOT > 0)
	{
	  GV.NP_TOT=POW3(GV.NP_TOT);
	  if(rank==0)
	    printf("Total number of particles: %lu\n", GV.NP_TOT);
	}
      else
	{
	  if(rank==0)
	    printf("Invalid 'NPARTICLES' setting in configuration file.\n");
	  return -2;
	}
    }
  
  /* NGRID parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'NGRID' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.NGRID=atoi(buf2);
      if(GV.NGRID > 0)
	{
	  GV.NGRID3=POW3(GV.NGRID);
	  if(rank==0)
	    printf("Number of grids in each axis: %d\n", GV.NGRID);
	}
      else
	{
	  if(rank==0)
	    printf("Invalid 'NGRID' setting in configuration file.\n");
	  return -2;
	}
    }
  
  /* Reading SCHEME parameter */  
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'SCHEME' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.SCHEME = strdup(buf2);
      if(strcmp(GV.SCHEME,ngp)    == 0 || strcmp(GV.SCHEME,cic)    == 0 ||
	 strcmp(GV.SCHEME,tsc)    == 0 || strcmp(GV.SCHEME,daub20) == 0 )
	{
	  if(rank==0)
	    printf("Scheme used: %s\n", GV.SCHEME);
	}
      else
	{
	  if(rank==0)
	    printf("Invalid 'SCHEME' setting in configuration file.\n");
	  return -2;
	}
    }

  /* Reading OUTPUT_DIR parameter */  
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      if(rank==0)
	printf("No 'OUTPUT_DIR' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.OUTPUT_DIR = strdup(buf2);
      if(rank==0)
	printf("Output directory: %s\n", GV.OUTPUT_DIR);
    }

  fclose(cfg);
  free(buf);
  free(buf1);
  free(buf2);

  return 0;
}



/*
 * Function:  read_head
 * --------------------
 * Reads the header of a gadget snapshot.
 *
 *  infile: String with the name of the snapshot.
 *  rank:   Rank of the process.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the header.
 *           -2 --> HALO masses are individual. As a condition of the 
 *                  code, HALO masses must be global.   
 */
int read_head(char *infile, int rank)
{  
  int i;
  unsigned int dummy;
  char label[4];
  FILE *fp_inp=NULL;
  int err;
  
  if( (fp_inp=fopen(infile,"r"))==NULL )
    {
      if(rank==0)
	printf("read_head: cannot open %s\n",infile);
      return -1;
    }
  
  
  if(GV.GADGET_VERSION==2)
    { 
      err=fread(&dummy, sizeof(dummy), 1, fp_inp);
      err=fread(&label, sizeof(char),  4, fp_inp);
      err=fread(&dummy, sizeof(dummy), 1, fp_inp);
      err=fread(&dummy, sizeof(dummy), 1, fp_inp);
    }

  err=fread(&dummy,   sizeof(dummy),   1, fp_inp);
  err=fread(&Gheader, sizeof(Gheader), 1, fp_inp);
  err=fread(&dummy,   sizeof(dummy),   1, fp_inp);
  
  fclose(fp_inp);

  /***************************************************************/
  /* Storing parameters of the simulation in the gloval structure */
  /***************************************************************/
  
  GV.L = Gheader.BoxSize;                       // Simulation lenght in Mpc
  GV.H = GV.L / (1.0*GV.NGRID);                 // Size of the cell
  GV.VOL_CELL = POW3(GV.H);                     // Volume of each cell
  GV.TOTAL_MASS = GV.NP_TOT*Gheader.mass[HALO]; // Total simulation mass 
  GV.RHO_MEAN = GV.TOTAL_MASS / POW3(GV.L);     // Mean density of ALL the simulation
  
  /* Cosmological parameters  */
  //GV.OMEGA_M0 = Header.Omega0;           //Omega matter at present time
  //GV.OMEGA_L0 = Header.OmegaLambda;      //Omega Lambda at present time
  //GV.ZRS = Header.redshift;              //Redshift of the simulation
  //GV.HUBBLEPARAM = Header.HubbleParam;   //Hubble parameter of the simulation
  
  //if(rank==0){printf("\n");}

  for(i=0; i<6; i++)
    { 
      if((Gheader.npart[i] != 0) && (Gheader.mass[i] != 0.0))
	{
	  if(rank==0)
	    {
	      printf(" * The mass of each particle is %d es %g\n",i,Gheader.mass[i]);
	    }
	}
      
      if((Gheader.npart[i] != 0) && (Gheader.mass[i] == 0.0))
	{
	  if(rank==0)
	    {
	      printf(" * There are individual mases for this particle set %d\n",i);
	    }
	  
	  if(i==HALO)
	    {
	      if(rank==0)
		{
		  printf("read_head Error: HALO masses must be global not individual\n");
		}
	      return -2;
	    }
	}
    } 
  
  if(rank==0)
    {
      printf("\n");
          
      printf(" * Frame's Time... %g\n", Gheader.time); 
      printf(" * Redshift... %g\n",     Gheader.redshift);
      printf(" * Flagsfr... %d\n",      Gheader.flag_sfr);
      printf(" * Flagfed... %d\n",      Gheader.flag_feedback);
      
      //printf("\n");
      
      //GV.NPART_TOTAL = 0;
      //for(i=0; i<6; i++)
      //{
      //printf(" * Header nall[%d] is: %d\n",i,Gheader.npartTotal[i]);
      //GV.NPART_TOTAL += Gheader.npartTotal[i];
      //}
      
      //GV.NPART_SELECTED = (unsigned long) llround( (GV.PERCENTAGE_TAKEN/100.0)*GV.NPART_TOTAL );
      
      printf("\n");
      
      printf(" * Flagcool...       %d\n",  Gheader.flag_cooling);
      printf(" * numfiles...       %d\n",  Gheader.num_files);
      printf(" * Boxsize...        %g\n",  Gheader.BoxSize);
      printf(" * Omega0...         %g\n",  Gheader.Omega0);
      printf(" * OmageLa...        %g\n",  Gheader.OmegaLambda);
      printf(" * Hubbleparam...    %g\n",  Gheader.HubbleParam);
      printf(" * NPART_TOTAL...    %lu\n", GV.NP_TOT);
      
    }
  
  if(err){}
  
  return 0;
}



/*
 * Function:  W_NGP
 * --------------------
 * Computes the window function for the Nearest GRID POINT (NGP) scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: window function value according to the distribution scheme.
 */
double W_NGP(double x, double y, double z, double H)
{
  //double Wx, Wy, Wz;
  
  /* Nearest Grid Point */
  // One dimensional window function in the X-axis
  /*
    if( fabs(x) < H*0.5 ){
    Wx = 1.0;
    }else if( fabs(x) == H*0.5  ){
    Wx = 0.5;
    }else{
    Wx = 0.0;
    }
  */
  
  // One dimensional window function in the Y-axis
  /*
    if( fabs(y) < H*0.5 ){
    Wy = 1.0;
    }else if( fabs(y) == H*0.5  ){
    Wy = 0.5;
    }else{
    Wy = 0.0;
    }
  */
  
  // One dimensional window function in the Z-axis
  /*
    if( fabs(z) < H*0.5 ){
    Wz = 1.0;
    }else if( fabs(z) == H*0.5  ){
    Wz = 0.5;
    }else{
    Wz = 0.0;
    }
  */

  
  //return Wx * Wy * Wz;
  /* As we use a regular cubic grid, the 
     three dimensional window function is 
     given as the multiplication of three 
     one dimensional window functions.
     WindowFunction = W(x) * W(y) * W(z).
  */
  return 1.0;
}



/*
 * Function:  W_CIC
 * --------------------
 * Computes the window function for the Cloud In Cell (CIC) scheme.
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: window function value according to the distribution scheme.
 */
double W_CIC(double x, double y, double z, double H)
{
  double Wx, Wy, Wz;
  
  /* Cloud In Cell */
  // One dimensional window function in the X-axis
  if( fabs(x) < H )
    {
      Wx = 1.0 - fabs(x)/H;
    }
  else
    {
      Wx = 0.0;
    }
  
  // One dimensional window function in the Y-axis
  if( fabs(y) < H )
    {
      Wy = 1.0 - fabs(y)/H;
    }
  else
    {
      Wy = 0.0;
    }
  
  // One dimensional window function in the Z-axis
  if( fabs(z) < H )
    {
      Wz = 1.0 - fabs(z)/H;
    }
  else
    {
      Wz = 0.0;
    }
  
  return Wx * Wy * Wz; /* As we use a regular cubic grid, the 
			  three dimensional window function is 
			  given as the multiplication of three 
			  one dimensional window functions.
			  WindowFunction = W(x) * W(y) * W(z).
		       */
}



/*
 * Function:  W_TSC
 * --------------------
 * Computes the window function for the Triangular Shaped Cloud (TSC) 
 * scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: window function value according to the distribution scheme.
 */
double W_TSC(double x, double y, double z, double H)
{
  double Wx, Wy, Wz;
  
  /* Triangular Shaped Cloud */
  // One dimensional window function in the X-axis
  if( fabs(x) <= H*0.5 )
    {
      Wx = 0.75 - ( (x*x)/(H*H) );
    }
  else if( H*0.5 <= fabs(x) && fabs(x) <= 1.5*H  )
    {
      Wx = 0.5*(1.5 - fabs(x)/H)*(1.5 - fabs(x)/H);
    }
  else
    {
      Wx = 0.0;
    }
  
  // One dimensional window function in the Y-axis
  if( fabs(y) <= H*0.5 )
    {
      Wy = 0.75 - ( (y*y)/(H*H) );
    }
  else if( H*0.5 <= fabs(y) && fabs(y) <= 1.5*H  )
    {
      Wy = 0.5*(1.5 - fabs(y)/H)*(1.5 - fabs(y)/H);
    }
  else
    {
      Wy = 0.0;
    }
  
  // One dimensional window function in the Z-axis
  if( fabs(z) <= H*0.5 )
    {
      Wz = 0.75 - ( (z*z)/(H*H) );
    }
  else if( H*0.5 <= fabs(z) && fabs(z) <= 1.5*H  )
    {
      Wz = 0.5*(1.5 - fabs(z)/H)*(1.5 - fabs(z)/H);
    }
  else
    {
      Wz = 0.0;
    }
  
  return Wx * Wy * Wz; /* As we use the regular cubic grid, the 
			  three dimensional window function is 
			  given as the multiplication of three 
			  one dimensional window functions.
			  WindowFunction = W(x) * W(y) * W(z).
		       */
}



/*
 * Function:  W_D20
 * --------------------
 * Computes the window function for the Daubechies D20 scaling 
 * function scheme. 
 *
 *  x, y, z: position in x, y and z.
 *
 *  H: separation of grids ( H = L / N_grid ).
 *
 *  returns: value of the window function according to the distribution 
 *           scheme.
 */
double W_D20(double x, double y, double z, double H)
{
  double Wx, Wy, Wz;
  
  // One dimensional window function in the X-axis
  Wx = gsl_spline_eval(spline, fabs(x)/GV.H, acc);
  
  // One dimensional window function in the Y-axis
  Wy = gsl_spline_eval(spline, fabs(y)/GV.H, acc);
  
  // One dimensional window function in the Z-axis
  Wz = gsl_spline_eval(spline, fabs(z)/GV.H, acc);
  
  return Wx * Wy * Wz; /* As we use the regular cubic grid, the 
			  three dimensional window function is 
			  given as the multiplication of three 
			  one dimensional window functions.
			  WindowFunction = W(x) * W(y) * W(z).
		       */
}



/*
 * Function:  mod
 * --------------------
 * Calculate the modulo operation for two numbers a and b (a%b) 
 * includding negative numbers.
 *
 *  a: Numerator of the division.
 *  b: Denominator of the division.
 *
 *  returns: The modulo a%b includding the option for negative numbers.
 */
int mod(int a, int b)
{
  int mod = a%b;
  while(mod<0)
    {
      mod += b;
    }
  return mod;
}
