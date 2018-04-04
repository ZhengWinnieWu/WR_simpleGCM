/*
  mppnccombine - joins together netCDF (eventually also plain IEEE) data
                 files representing a decomposed domain into a unified
                 netCDF file.  It was originally designed to be used as a
                 postprocessor for the parallel I/O programming interface
                 "mpp_io_mod" (http://www.gfdl.gov/~vb/mpp_io.html) by
                 V. Balaji.

  V1.1.1: Added a fix for dimensions that are not also variables.
  V1.1: Changed loop order for increased I/O efficiency; records are now the
        innermost loop then the variables loop.
  V1.0: Original release.

  Written by Hans Vahlenkamp (hnv@gfdl.gov)
  Geophysical Fluid Dynamics Laboratory/NOAA
  Princeton University Forrestal Campus
  Last Updated: 4/25/01
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <netcdf.h>

/* Information structure for a input file */
struct fileinfo
  {
   int ncfid;  /* ID of the input netCDF file */
   int ndims;  /* Number of dimensions */
   int nvars;  /* Number of variables */
   int ngatts;  /* Number of global attributes */
   int recdim;  /* IDs of the record dimensions */
   char varname[MAX_NC_VARS][MAX_NC_NAME];  /* Names of the variables */
   nc_type datatype[MAX_NC_VARS]; /* Data types of the variables */
   int varndims[MAX_NC_VARS];  /* Number of dimensions for each variable */
   int vardim[MAX_NC_VARS][MAX_NC_DIMS];  /* Dimensions for each variable */
   int natts[MAX_NC_VARS];  /* Number of attributes for each variable */
   char dimname[MAX_NC_DIMS][MAX_NC_NAME];  /* Names of the dimensions */
   long dimsize[MAX_NC_DIMS];  /* Sizes of the dimensions (decomposed) */
   long dimfullsize[MAX_NC_DIMS];  /* Full sizes of the dimensions */
   long dimstart[MAX_NC_DIMS];  /* Start positions within full dimensions */
   long dimend[MAX_NC_DIMS];  /* End positions within full dimensions */
  };

/* Auxiliary function prototypes */
void usage();
int process_ncinfile(char *, unsigned char, int, char *, int *,
                     unsigned char);
int copy_nc_data(struct fileinfo *, int, unsigned char, unsigned char);
#if DEBUG==1
void print_debug(struct fileinfo *, unsigned char);
char *nc_type_to_str(nc_type);
#endif


int main(int argc, char *argv[])
  {
   unsigned char verbose=0;  /* Print some progress information? */
   unsigned char appendnc=0;  /* Append to an existing netCDF file? */
   unsigned char removein=0;  /* Remove the ".####" decomposed input files? */
   int outputarg=(-1);  /* Argument # of the output netCDF file */
   int inputarg=(-1);  /* Argument # of first input netCDF file */
   struct stat statbuf;  /* Dummy structure for file-testing "stat" call */
   int outncfid;  /* ID of the output netCDF file */
   char outfilename[2048], *strptr;  /* Name of the output netCDF file */
   int outlen;  /* Length of the output filename */
   char infilename[2048];  /* Name of an input file */
   unsigned char infileerror=0;  /* Were there any file errors? */
   int nfiles=(-1), nfiles2;  /* Number of files in the decomposed domain */
   int a;  /* Loop variable */

   /* Check the command-line arguments */
   if (argc < 2)
     {
      usage(); return(1);
     }
   for (a=1; a < argc; a++)
     {
      if (!strcmp(argv[a],"-v")) verbose=1;
      else if (!strcmp(argv[a],"-a")) appendnc=1;
      else if (!strcmp(argv[a],"-r")) removein=1;
      else
        {
         outputarg=a; break;
        }
     }
   if (outputarg==(-1))
     {
      usage(); return(1);
     }
   if (argc-1 > outputarg) inputarg=outputarg+1;
   sprintf(outfilename,argv[outputarg]); outlen=strlen(outfilename);
   if (outlen > 4)
     {
      strptr=outfilename+outlen-5;
      if (!strcmp(strptr,".0000")) outfilename[outlen-5]='\0';
     }

   /* Disable fatal returns from netCDF library functions */
   ncopts=0;

   /* Create a new netCDF output file */
   if (!appendnc)
     {
      if (stat(outfilename,&statbuf)==0)
        {
         fprintf(stderr,"Error: output file seems to exist already!\n");
         return(1);
        }
      if ((outncfid=nccreate(outfilename,NC_NOCLOBBER))==(-1))
        {
         fprintf(stderr,"Error: cannot create the output netCDF file!\n");
         return(1);
        }
      ncsetfill(outncfid,NC_NOFILL);
     }
   /* Open an existing netCDF file for appending */
   else
     {
      if ((outncfid=ncopen(outfilename,NC_WRITE))==(-1))
        {
         fprintf(stderr,"Error: cannot open the output netCDF file for appending!\n");
         return(1);
        }
     }

   /* No input files are specified on the command-line */
   if (inputarg==(-1))
     {
      nfiles=1;
      for (a=0; a < nfiles; a++)
        {
         sprintf(infilename,"%s.%04d",outfilename,a);
         if (verbose) printf("Processing... \"%s\"",infilename);
#if DEBUG==1
         else if (!verbose) printf("\nfile=%s\n",infilename);
#endif
         if (stat(infilename,&statbuf)!=0)
           {
            if (verbose) printf("\n");
            fprintf(stderr,"Error: cannot read the input file \"%s\"!\n",
                    infilename);
            if (a==0) infileerror=3;
            else infileerror=1;
            break;
           }
         infileerror=process_ncinfile(infilename,appendnc,outncfid,
                                      outfilename,&nfiles2,verbose);
         if (infileerror==2)
           {
            printf("IEEE input files are not currently supported!\n");
            infileerror=3; break;
           }
         else if (infileerror==1) break;
         if (a==0 && nfiles2 > 0) nfiles=nfiles2;
         appendnc=1;
         if (verbose) printf("\n (%d files to go)\n",nfiles-a-1);
        }
     }
   /* Loop over all the specified input files */
   else
     for (a=inputarg; a < argc; a++)
       {
        if (verbose) printf("Processing... \"%s\"",argv[a]);
#if DEBUG==1
        else if (!verbose) printf("\nfile=%s\n",argv[a]);
#endif
        if (stat(argv[a],&statbuf)!=0)
          {
           if (verbose) printf("\n");
           fprintf(stderr,"Error: cannot read the input file \"%s\"!\n",
                   argv[a]);
           if (a==inputarg) infileerror=3;
           else infileerror=1;
           break;
          }
        infileerror=process_ncinfile(argv[a],appendnc,outncfid,outfilename,
                                     &nfiles2,verbose);
        if (infileerror==2)
          {
           printf("IEEE input files are not currently supported!\n");
           infileerror=3; break;
          }
        else if (infileerror==1) return(1);
        appendnc=1;
        if (verbose) printf("\n (%d files to go)\n",argc-a-1);
       }

   /* Clean up... return 1 on error, otherwise 0 */
   ncclose(outncfid);
   if (!infileerror)
     {
      if (removein)
        {
         /* No input files are specified on the command-line */
         if (inputarg==(-1))
           {
            for (a=0; a < nfiles; a++)
              {
               sprintf(infilename,"%s.%04d",outfilename,a);
               if (verbose)
                 printf("Removing... \"%s\" (%d files to go)\n",infilename,
                        nfiles-a-1);
               unlink(infilename);
              }
           }
         /* Loop over all the specified input files */
         else
           for (a=inputarg; a < argc; a++)
             {
              printf("Removing... \"%s\" (%d files to go)\n",argv[a],
                     argc-a-1);
              unlink(argv[a]);
             }
        }
     }
   else if (infileerror==3)
     {
      unlink(outfilename); infileerror=1;
     }
   return(infileerror);
  }


/* Print the usage message for mppnccombine */
void usage()
  {
   printf("Usage:  mppnccombine [-v] [-a] [-r] output.nc [input ...]\n\n");
   printf("   -v   Print some progress information.\n");
   printf("   -a   Append to an existing netCDF file.\n");
   printf("   -r   Remove the \".####\" decomposed files after a successful run.\n\n");
   printf("mppnccombine joins together an arbitrary number of data files (netCDF only for\n");
   printf("now, but IEEE binary may be supported in the future), containing chunks of a\n");
   printf("decomposed domain, into a unified netCDF file.  An output file must be\n");
   printf("specified and it is assumed to be the first filename argument.  If the output\n");
   printf("file already exists, then it will not be modified unless the option is chosen\n");
   printf("to append to it.  If no input files are specified then their names will be\n");
   printf("based on the name of the output file plus the extensions \".0000\", \".0001\", etc.\n");
   printf("If input files are specified then they are assumed to be absolute filenames.\n\n");
   printf("A value of 0 is returned if execution completed successfully; a value of 1\n");
   printf("otherwise.\n");
  }


/* Open an input netCDF file and get some information about it */
int process_ncinfile(char *ncname, unsigned char appendnc, int outncfid,
                     char *outncname, int *nfiles, unsigned char verbose)
  {
   struct fileinfo ncinfile;  /* Information about an input netCDF file */
   int nfiles2;  /* Number of files in the decomposed domain */
   int d, v, n;  /* Loop variables */
   int dimid;  /* ID of a dimension */
   int decomp[4];  /* "domain_decomposition" information */
   char attname[MAX_NC_NAME];  /* Name of a global or variable attribute */
   unsigned char ncinfileerror=0;  /* Were there any file errors? */

   /* Open an input netCDF file; return if not openable - possibly IEEE */
   if ((ncinfile.ncfid=ncopen(ncname,NC_NOWRITE))==(-1)) return(2);

   /* Determine the number of files in the decomposed domain */
   if (ncattget(ncinfile.ncfid,NC_GLOBAL,"NumFilesInSet",
                (void *)&nfiles2)==(-1))
     {
      if (*nfiles==1)
        {
         fprintf(stderr,"Error: missing the \"NumFilesInSet\" global attribute!\n");
         return(1);
        }
      else if (*nfiles==(-1))
        {
         fprintf(stderr,"Warning: missing the \"NumFilesInSet\" global attribute.\n");
        }
     }
   *nfiles=nfiles2;

   /* Get some general information about the input netCDF file */
   if (ncinquire(ncinfile.ncfid,&(ncinfile.ndims),&(ncinfile.nvars),
                 &(ncinfile.ngatts),&(ncinfile.recdim))==(-1))
     {
      fprintf(stderr,"Error: cannot read the file's metadata!\n");
      ncclose(ncinfile.ncfid); return(1);
     }

   /* Get some information about the dimensions */
   for (d=0; d < ncinfile.ndims; d++)
     {
      if ((ncdiminq(ncinfile.ncfid,d,ncinfile.dimname[d],
                    &(ncinfile.dimsize[d])))==(-1))
        {
         fprintf(stderr,"Error: cannot read dimension #%d's metadata!\n",d);
         ncclose(ncinfile.ncfid); return(1);
        }
      ncinfile.dimfullsize[d]=ncinfile.dimsize[d];
      ncinfile.dimstart[d]=1; ncinfile.dimend[d]=(-1);
     }

   /* Get some information about the variables */
   for (v=0; v < ncinfile.nvars; v++)
     {
      if ((ncvarinq(ncinfile.ncfid,v,ncinfile.varname[v],
                    &(ncinfile.datatype[v]),&(ncinfile.varndims[v]),
                    ncinfile.vardim[v],&(ncinfile.natts[v])))==(-1))
        {
         fprintf(stderr,"Error: cannot read variable #%d's metadata!\n",v);
         ncclose(ncinfile.ncfid); return(1);
        }

      /* If the variable is also a dimension then get decomposition info */
      if ((dimid=ncdimid(ncinfile.ncfid,ncinfile.varname[v]))!=(-1))
        {
         if (ncattget(ncinfile.ncfid,v,"domain_decomposition",
             (void *)decomp)!=(-1))
           {
            ncinfile.dimfullsize[dimid]=decomp[1]-decomp[0]+1;
            ncinfile.dimstart[dimid]=decomp[2]-(decomp[0]-1);
            ncinfile.dimend[dimid]=decomp[3]-(decomp[0]-1);
           }
         else
           {
            ncinfile.dimfullsize[dimid]=ncinfile.dimsize[dimid];
            ncinfile.dimstart[dimid]=1; ncinfile.dimend[dimid]=(-1);
           }
        }
     }

#if DEBUG==1
   print_debug(&ncinfile,verbose);
#endif

   /* If the output netCDF file was just created then define its structure */
   if (!appendnc)
     {
#if DEBUG==1
      printf("Creating output netCDF file... \"%s\"\n",outncname);
#endif
      /* Define the dimensions */
      for (d=0; d < ncinfile.ndims; d++)
        {
         if (d==ncinfile.recdim)
           ncdimdef(outncfid,ncinfile.dimname[d],NC_UNLIMITED);
         else ncdimdef(outncfid,ncinfile.dimname[d],ncinfile.dimfullsize[d]);
        }

      /* Define the variables and copy their attributes */
      for (v=0; v < ncinfile.nvars; v++)
        {
         ncvardef(outncfid,ncinfile.varname[v],ncinfile.datatype[v],
                  ncinfile.varndims[v],ncinfile.vardim[v]);
         for (n=0; n < ncinfile.natts[v]; n++)
           {
            ncattname(ncinfile.ncfid,v,n,attname);
            if (!strcmp(attname,"domain_decomposition")) continue;
            else
              {
               if (ncattcopy(ncinfile.ncfid,v,attname,outncfid,v)==(-1))
                 {
                  fprintf(stderr,"Error: cannot copy variable \"%s\"'s attributes!\n",
                          ncinfile.varname[v]);
                  return(1);
                 }
              }
           }
        }

      /* Copy the global attributes */
      for (n=0; n < ncinfile.ngatts; n++)
        {
         ncattname(ncinfile.ncfid,NC_GLOBAL,n,attname);
         if (!strcmp(attname,"NumFilesInSet")) continue;
         else if (!strcmp(attname,"filename"))
           ncattput(outncfid,NC_GLOBAL,attname,NC_CHAR,strlen(outncname),
                    (void *)outncname);
         else
           {
            if (ncattcopy(ncinfile.ncfid,NC_GLOBAL,attname,outncfid,
                          NC_GLOBAL)==(-1))
              {
               fprintf(stderr,"Error: cannot copy the file's global attributes!\n");
               return(1);
              }
           }
        }

      /* Definitions done */
      ncendef(outncfid);
     }

   /* Copy all the data values of the dimensions and variables */
   ncinfileerror=copy_nc_data(&ncinfile,outncfid,appendnc,verbose);

   /* Done */
   ncclose(ncinfile.ncfid); return(ncinfileerror);
  }


/* Copy all the data values in an input netCDF file */
int copy_nc_data(struct fileinfo *ncinfile, int outncfid,
                 unsigned char appendnc, unsigned char verbose)
  {
   int v, d, r;  /* Loop variables */
   int dimid;  /* ID of a dimension */
   void *values;  /* Data values to copy */
   long instart[MAX_NC_DIMS], outstart[MAX_NC_DIMS];  /* Copy array sizes */
   long count[MAX_NC_DIMS];                           /*        "         */
   long nrecs;  /* Number of records */
   long recsize;  /* Number of values in a record */
   int varrecdim;  /* Position of a variable's record dimension */

   /* Loop over all the records */
   nrecs=ncinfile->dimsize[ncinfile->recdim];
   for (r=0; r < nrecs; r++)
     {
#if DEBUG==0
      if (verbose) printf("\n  record=%d",r+1);
#endif

      /* Loop over all the variables */
      for (v=0; v < ncinfile->nvars; v++)
        {
#if DEBUG==0
         if (verbose) printf("\n    variable=%s",ncinfile->varname[v]);
#endif

         /* Avoid multiple reads/writes of non-decomposed dimensions */
         if ((dimid=ncdimid(ncinfile->ncfid,ncinfile->varname[v]))!=(-1))
           if (appendnc && ncinfile->dimend[dimid]==(-1)) continue;

         /* Get read/write dimension sizes for the variable */
         recsize=1; varrecdim=(-1);
         for (d=0; d < ncinfile->varndims[v]; d++)
           {
            if (ncinfile->vardim[v][d]==ncinfile->recdim)
              {
               count[d]=1; varrecdim=d;
              }
            else
              {
               count[d]=ncinfile->dimsize[ncinfile->vardim[v][d]];
               recsize*=count[d]; instart[d]=0;
               outstart[d]=ncinfile->dimstart[ncinfile->vardim[v][d]]-1;
              }
#if DEBUG==1
            printf("%d:  instart=%ld  ouststart=%ld  count=%ld\n",d,
                   instart[d],outstart[d],count[d]);
#endif
           }

         /* Avoid multiple reads/writes of non-record variables */
         if (varrecdim==(-1) && r > 0) continue;

#if DEBUG==0
         if (verbose) printf("  (read/write)");
#endif

         /* Allocate a buffer for the variable's record */
         if ((values=malloc(nctypelen(ncinfile->datatype[v])*recsize))==NULL)
           {
            fprintf(stderr,"Error: cannot allocate memory for variable \"%s\"'s values!\n",
                    ncinfile->varname[v]);
            return(1);
           }

         /* Copy the record */
         if (varrecdim!=(-1)) instart[varrecdim]=outstart[varrecdim]=r;
         if (ncvarget(ncinfile->ncfid,v,instart,count,values)==(-1))
           {
            fprintf(stderr,"Error: cannot read variable \"%s\"'s values!\n",
                    ncinfile->varname[v]);
            return(1);
           }
         if (ncvarput(outncfid,v,outstart,count,values)==(-1))
           {
            fprintf(stderr,"Error: cannot write variable \"%s\"'s values!\n",
                    ncinfile->varname[v]);
            return(1);
           }

         /* Deallocate the record buffer */
         free(values);
        }
     }
   return(0);
  }


#if DEBUG==1

/* Print some debugging information */
void print_debug(struct fileinfo *infile, unsigned char verbose)
  {
   int v, d;  /* Loop variables */

   if (verbose) printf("\n");
   printf("ncfid=%d  ndims=%d  nvars=%d  ngatts=%d  recdim=%d\n",
          infile->ncfid,infile->ndims,infile->nvars,infile->ngatts,
          infile->recdim);
   printf("variables:\n");
   for (v=0; v < infile->nvars; v++)
     {
      printf("  varname=%s  datatype=%s  varndims=%d  natts=%d\n",
             infile->varname[v],nc_type_to_str(infile->datatype[v]),
             infile->varndims[v],infile->natts[v]);
      printf("    vardim =");
      for (d=0; d < infile->varndims[v]; d++)
        printf(" %d",infile->vardim[v][d]);
      printf("\n");
     }
   printf("dimensions:\n");
   for (d=0; d < infile->ndims; d++)
     {
      printf("  dimname=%s  dimsize=%ld",infile->dimname[d],
             infile->dimsize[d]);
      printf("  dimfullsize=%ld  dimstart=%ld  dimend=%ld\n",
             infile->dimfullsize[d],infile->dimstart[d],infile->dimend[d]);
     }
  }


/* Convert a netCDF datatype number to a string */
char *nc_type_to_str(nc_type datatype)
  {
   switch (datatype)
     {
      case NC_BYTE:
        return("NC_BYTE"); break;
      case NC_CHAR:
        return("NC_CHAR"); break;
      case NC_SHORT:
        return("NC_SHORT"); break;
      case NC_LONG:
        return("NC_LONG"); break;
      case NC_FLOAT:
        return("NC_FLOAT"); break;
      case NC_DOUBLE:
        return("NC_DOUBLE"); break;
      default:
        return("Private!"); break;
     }
  }

#endif
