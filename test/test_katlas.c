/* 

   test_katlas: 

   This test matches .katlas files with .pdstor files to make sure that the diagrams
   in the katlas files are represented in the .pdstor files.

*/




#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_SYS_STAT_H
  #include<sys/stat.h>
#endif

#ifdef HAVE_SYS_TYPES_H
  #include<sys/types.h>
#endif

#ifdef HAVE_ASSERT_H
  #include<assert.h>
#endif

#ifdef HAVE_STDBOOL_H
  #include<stdbool.h>
#endif

#ifdef HAVE_STDLIB_H
  #include<stdlib.h>
#endif

#ifdef HAVE_STDIO_H
  #include<stdio.h>
#endif

#ifdef HAVE_TIME_H
  #include<time.h>
#endif

#ifdef HAVE_STDINT_H
  #include<stdint.h>
#endif

#ifdef HAVE_ARGTABLE2_H
  #include<argtable2.h>
#endif

#define MAXVERTS 12
#define MAXEDGES 25

int VERBOSE;
int PD_VERBOSE;

#include<ordie.h>
#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>
#include<pd_storage.h>

#define DEBUG 1            /* Turn on asserts */
int PD_VERBOSE={0};       /* Turn on debugging info */

struct arg_lit  *verbose;
struct arg_lit  *help;
struct arg_end  *end;
struct arg_end  *helpend;
  
int main(int argc,char *argv[])
{
  int            nerrors,ncross;
   
  void *argtable[] = 
    {
     verbose = arg_lit0("v","verbose","print debugging information"),
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};

  printf("test_katlas (%s)\n",PACKAGE_STRING);

 /* First, we parse the command-line arguments using argtable. */

  if (arg_nullcheck(argtable) != 0)
    printf("test_katlas: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"test_katlas");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("test_katlas matches .katlas files with .pdstor files and checks to \n"
	     "see if all the pd codes in the .katlas files are found somewhere in \n"
	     "the corresponding .pdstor file. This is a sanity check to make sure \n"
	     "that our algorithm is at least generating all the diagrams we already\n"
	     "know about.\n"
	     "usage: \n\n");

      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (verbose->count > 0) {

    VERBOSE = 50;

  } else { 

    VERBOSE = 0;

  }

  clock_t big_start,big_end;
  double  big_cpu_time_used;

  big_start = clock();

  for(ncross=3;ncross < 10;ncross++) {

    FILE         *katlas;
    FILE         *pdstor;
    char         katlas_name[256],pdstor_name[256];
    clock_t      start,end;
    double       cpu_time_used;

    sprintf(katlas_name,"%d.katlas",ncross);
    sprintf(pdstor_name,"%d.pdstor",ncross);

    katlas = fopen(katlas_name,"r");
    pdstor = fopen(pdstor_name,"r");

    if (katlas != NULL && pdstor != NULL) { 

      printf("Matched files %s and %s. Starting scan.\n",katlas_name,pdstor_name);
      
      PD_VERBOSE = 60;
      pd_stor_t *existing_pdstor = pd_read_pdstor(pdstor);
      PD_VERBOSE = VERBOSE;

      unsigned int nhashes, nelts;
      pd_stor_stats(existing_pdstor,&nhashes,&nelts);
      printf("\tRead %d hash, %d element pdstor from %s.\n",nhashes,nelts,pdstor_name);

      int katlas_ncross,katlas_ncodes;
      if (fscanf(katlas," katlas ncrossings %d ncodes %d ",&katlas_ncross, &katlas_ncodes) != 2) {

	printf("=============================\n");
	printf("Failed to read %s. Test FAIL.\n",katlas_name);
	exit(1);

      }
      if (katlas_ncross != ncross) { 

	printf("===========================\n");
	printf("Corrupt %s file. Test FAIL.\n",katlas_name);
	exit(1);

      }
      
      int katlas_code;
      for(katlas_code = 0;katlas_code < katlas_ncodes;katlas_code++) { 

	pd_code_t katlas_pd;
	int cr;

	/* Read the crossing data from this pdcode in the katlas file. */

	katlas_pd.ncross = katlas_ncross;

	for(cr=0;cr<katlas_ncross;cr++) { 

	  int epos;
	  for(epos = 0;epos < 4;epos++) {

	    if (fscanf(katlas,"%hi",&(katlas_pd.cross[cr].edge[epos])) != 1) { 

	      printf("=========================\n");
	      printf("Corrupt %s file. Test FAIL\n",katlas_name);
	      exit(1);

	    }
	  }
	}
	
	/* Now read the knot name as a string */

	char katlas_code_name[256];
	
	if (fscanf(katlas,"%200s",katlas_code_name) != 1) { 

	   printf("=========================\n");
	   printf("Corrupt %s file. Test FAIL\n",katlas_name);
	   exit(1);
	   
	}

	/* Now regenerate the pd_code_t katlas_pd from this crossing info. */

	pd_regenerate(&katlas_pd);

	/* We now search the pdstor for an isomorphic copy of this pd. */

	pd_iso_t **isos;
	unsigned int nisos;
	pd_code_t *matching_pd;

	matching_pd = pd_search_pdstor_by_isomorphism(existing_pdstor,&katlas_pd,
						      &isos,&nisos);

	if (matching_pd == NULL) { 

	  printf("=====================================================\n");
	  printf("Couldn't match KnotAtlas diagram for %s given by pd\n",katlas_code_name);
	  pd_write(stdout,&katlas_pd);
	  printf("\nTest FAIL.\n");
	  exit(1);

	}
	  
	printf("\t%s -> Present.\n",katlas_code_name);
	pd_free_isos(&nisos,&isos);
	  
      }	
      
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      
      printf("Found all elements of %s in %s. Elapsed time %4.2g seconds.\n",
	     katlas_name,pdstor_name,cpu_time_used);
      
      fclose(katlas); 
      fclose(pdstor);
      
    }

  }
	
  big_end = clock();
  big_cpu_time_used = ((double) (big_end - big_start)) / CLOCKS_PER_SEC;
    
  printf("PASS. Elapsed time %4.2g seconds.\n",
	 big_cpu_time_used);
 
  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
  free(helpend);

  exit(0);

}
