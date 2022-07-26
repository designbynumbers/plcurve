/* tangle interior test Z assembled c code */

pd_code_t *pd_create_tangle_slide_input_testZ_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(5);
assert(pd != NULL);
pd->ncross = 5;
pd->nedges = 10;
pd->ncomps = 3;
pd->nfaces = 7;

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 3;
pd->cross[0].edge[2] = 3;
pd->cross[0].edge[3] = 2;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 6;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 5;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 6;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 7;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 4;
pd->cross[3].edge[1] = 8;
pd->cross[3].edge[2] = 5;
pd->cross[3].edge[3] = 9;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 4;
pd->cross[4].edge[1] = 9;
pd->cross[4].edge[2] = 7;
pd->cross[4].edge[3] = 8;
pd->cross[4].sign = 0;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 0;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 0;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 1;

pd->edge[4].head = 3;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 0;

pd->edge[5].head = 1;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 3;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 2;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 4;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 2;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 3;
pd->edge[8].headpos = 1;
pd->edge[8].tail = 4;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 4;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 3;
pd->edge[9].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 4;
pd->comp[1].edge[1] = 5;
pd->comp[1].edge[2] = 6;
pd->comp[1].edge[3] = 7;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 8;
pd->comp[2].edge[1] = 9;


/* Face data */

pd->face[0].nedges = 6;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].orient = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].orient != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].orient[0] = 1;

pd->face[0].edge[1] = 5;
pd->face[0].orient[1] = 0;

pd->face[0].edge[2] = 8;
pd->face[0].orient[2] = 0;

pd->face[0].edge[3] = 7;
pd->face[0].orient[3] = 0;

pd->face[0].edge[4] = 2;
pd->face[0].orient[4] = 1;

pd->face[0].edge[5] = 3;
pd->face[0].orient[5] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].orient = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].orient != NULL);

pd->face[1].edge[0] = 1;
pd->face[1].orient[0] = 1;

pd->face[1].edge[1] = 7;
pd->face[1].orient[1] = 1;

pd->face[1].edge[2] = 9;
pd->face[1].orient[2] = 0;

pd->face[1].edge[3] = 5;
pd->face[1].orient[3] = 1;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].orient = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].orient != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].orient[0] = 0;

pd->face[2].edge[1] = 2;
pd->face[2].orient[1] = 0;

pd->face[2].edge[2] = 6;
pd->face[2].orient[2] = 0;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].orient = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].orient != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].orient[0] = 0;

pd->face[3].edge[1] = 6;
pd->face[3].orient[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].orient = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].orient != NULL);

pd->face[4].edge[0] = 4;
pd->face[4].orient[0] = 0;

pd->face[4].edge[1] = 8;
pd->face[4].orient[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].orient = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].orient != NULL);

pd->face[5].edge[0] = 4;
pd->face[5].orient[0] = 1;

pd->face[5].edge[1] = 9;
pd->face[5].orient[1] = 1;

pd->face[6].nedges = 1;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].orient = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].orient != NULL);

pd->face[6].edge[0] = 3;
pd->face[6].orient[0] = 1;


/* End of data. */

pd_regenerate_hash(pd);
assert(pd_ok(pd));
return pd;

}


bool pdint_check_tslide_data_ok_and_find_te_testZ() {

  /*                                                         */
/*                                                         */
/*                   /--3--\                               */
/*                   A     ^                               */
/*                   | (6) |c0                             */
/*                   \-3>--|-->-\          (0)             */
/*                         2    0                          */
/*                      c2 | (2)|c1                        */
/*                /----7-----<6----5----\                  */
/*                |        |    |       |                  */
/*                |  T     ^ (3)v       |                  */
/*                B  +---------------+  ^                  */
/*                |  |     \--1-/    |  |                  */
/*                |  |       (1)     |  |                  */
/*                7  |               |  |                  */
/*                |  |  /-C--9----\  |  |                  */
/*                v  +---------------+  |                  */
/*                |     v   (5)   ^     |                  */
/*                |     |         |     |                  */
/*                \7---------4-->-----5-/                  */
/*                   c4 |         |c3                      */
/*                      C   (4)   |                        */
/*                      \----8-->-/                        */
		       
pd_idx_t nedges = 4 ;            
pd_idx_t tangle_faces[4] = {3,1,5,1};
pd_idx_t tangle_edges[4] = {1,9,9,1};

pd_idx_t noverstrand_edges = 2;
pd_idx_t overstrand_edges[2] = {7,4};
pd_idx_t border_faces[2] = {1,5};

bool valid_ts_input = true ; 

pd_idx_t tangle_slide_edges[1] = {9};
pd_idx_t complementary_edges[3] = {1,1,9};
pd_boundary_or_t complementary_or[3] = {out,in,in};
bool overstrand_goes_OVER = true;
pd_or_t overstrand_orientation = PD_POS_ORIENTATION;


  printf("--------------------------------------------------\n"
    	 "pdint_check_tslide_data_ok_and_find_te test Z\n"
	     "--------------------------------------------------\n");

  printf("creating pd...");
  pd_code_t *pd = pd_create_tangle_slide_input_testZ_0();
  if (!pd_ok(pd)) {

    printf("fail (doesn't pass pd_ok)\n");
    return false;

  }
  
  printf("pass (passes pd_ok)\n");

  printf("creating basic tangle information...");
  pd_tangle_t *t = pd_tangle_new(nedges);

  pd_idx_t i;
  for(i=0;i<t->nedges;i++) {
  
    t->edge[i] = tangle_edges[i];
    t->face[i] = tangle_faces[i];

  }

  printf("done (tangle has %d edges)\n",nedges);

  printf("running pd_regenerate_tangle...");
  pd_regenerate_tangle(pd,t);

  if (!pd_tangle_ok(pd,t)) {

    printf("fail (doesn't pass pd_tangle_ok)\n");
    return false;
    
  }

  printf("pass (didn't crash, passes pd_tangle_ok)\n");
  printf("running pdint_check_tslide_data_ok_and_find_te ");

  pd_idx_t completely_fake_memory_address;
  pd_idx_t *computed_tangle_slide_edges = &completely_fake_memory_address;
  /* We want to initialize this to something that's not NULL, since we need it
     to be NULL if the input data is invalid. */

  bool computed_overstrand_goes_OVER;
  pd_or_t computed_overstrand_orientation;

  pd_idx_t *computed_complementary_edges = &completely_fake_memory_address;
  pd_boundary_or_t *computed_complementary_or = &completely_fake_memory_address;
  
  bool computed_ts_input;

  if (valid_ts_input) {

    printf("(expect true)...");
    
  } else {

    printf("(expect false)...");

  }



  computed_ts_input
     = pdint_check_tslide_data_ok_and_find_te(pd,t,noverstrand_edges,
                                              overstrand_edges,
                                              border_faces,
                                              &computed_tangle_slide_edges,
                                              &computed_complementary_edges,
                                              &computed_complementary_or,
                                              &computed_overstrand_goes_OVER,
                                              &computed_overstrand_orientation);
                                                             
  if (computed_ts_input != valid_ts_input) {

     printf("FAIL");

     if (computed_ts_input) {

        printf(" (got true)\n");

     } else {

        printf(" (got false)\n");

     }

     return false;

  }

  printf("pass");

  if (computed_ts_input) {

        printf(" (got true)\n");

  } else {

        printf(" (got false)\n");

  }

  if (computed_ts_input) { /* If there ARE tangle edges, try to match them... */

    /*************** Tangle edges *******************/

    printf("checking address of computed_tangle_slide_edges...");
    if (computed_tangle_slide_edges == &completely_fake_memory_address) {

      printf("FAIL (not updated in call)\n");
      return false;

    }

    if (computed_tangle_slide_edges == NULL) {

      printf("FAIL (set to NULL, even though input is valid)\n");
      return false;

    }

    printf("pass (set to a new memory address != NULL)\n");

    printf("checking computed tangle_slide_edges against expected...");

    pd_idx_t i;

    for(i=0;i<noverstrand_edges-1;i++) {

     if (computed_tangle_slide_edges[i] != tangle_slide_edges[i]) {

        printf("FAIL.\nExpected and computed tangle edges don't match at pos %d\n",i);

        pd_idx_t j;

        printf("Expected tangle edges: ");
        for(j=0;j<noverstrand_edges;j++) { printf("%4d ",tangle_slide_edges[j]); }
        printf("\nComputed tangle edges: ");
        for(j=0;j<noverstrand_edges;j++) { printf("%4d ",computed_tangle_slide_edges[j]); }
        printf("\n                       ");
        for(j=0;j<i;j++) { printf("     "); }
        printf("-----\n");

        return false;

      }

    }

    printf("pass (lists of edges match)\n");

    /*************** Complementary Edges ******************/

    printf("checking address of computed_complementary_edges...");
    if (computed_complementary_edges == &completely_fake_memory_address) {

      printf("FAIL (not updated in call)\n");
      return false;

    }

    if (computed_complementary_edges == NULL) {

      printf("FAIL (set to NULL, even though input is valid)\n");
      return false;

    }

    printf("pass (set to a new memory address != NULL)\n");

    printf("checking computed complementary_edges against expected...");

    pd_idx_t Ncomplementary = t->nedges-(noverstrand_edges-1);

    for(i=0;i<Ncomplementary;i++) {

      if (computed_complementary_edges[i] != complementary_edges[i]) {

        printf("FAIL.\nExpected and computed complementary edges don't match at pos %d\n",i);

        pd_idx_t j;

        printf("Expected complementary edges: ");
        for(j=0;j<Ncomplementary;j++) { printf("%4d ",complementary_edges[j]); }
        printf("\nComputed complementary edges: ");
        for(j=0;j<Ncomplementary;j++) { printf("%4d ",computed_complementary_edges[j]); }
        printf("\n                       ");
        for(j=0;j<i;j++) { printf("     "); }
        printf("-----\n");

        return false;

      }

    }

    printf("pass (lists of edges match)\n");


    /*************** Complementary Orientations ******************/

    printf("checking address of computed_complementary_or...");
    if (computed_complementary_or == &completely_fake_memory_address) {

      printf("FAIL (not updated in call)\n");
      return false;

    }

    if (computed_complementary_or == NULL) {

      printf("FAIL (set to NULL, even though input is valid)\n");
      return false;

    }

    printf("pass (set to a new memory address != NULL)\n");

    printf("checking computed complementary_or against expected...");

    for(i=0;i<Ncomplementary;i++) {

      if (computed_complementary_edges[i] != complementary_edges[i]) {

        printf("FAIL.\nExpected and computed complementary orientations don't match at pos %d\n",i);

        pd_idx_t j;

        printf("Expected complementary or: ");
        for(j=0;j<Ncomplementary;j++) {
           pd_printf("%BDY_OR ",NULL,complementary_or[j]); }
        printf("\nComputed complementary or: ");
        for(j=0;j<Ncomplementary;j++) {
           pd_printf("%BDY_OR ",NULL,computed_complementary_or[j]); }
        printf("\n                       ");
        for(j=0;j<i;j++) { printf("     "); }
        printf("-----\n");

        return false;

      }

    }

    printf("pass (lists match)\n");

    /************** overstrand_goes_OVER ***************/

    printf("checking computed overstrand_goes_OVER...");

    if (computed_overstrand_goes_OVER == overstrand_goes_OVER) {

       printf("pass (both are %s)\n",(overstrand_goes_OVER ? "true" : "false"));

    } else {

       printf("fail (expected %s, got %s)\n",
       (overstrand_goes_OVER ? "true" : "false"),
       (computed_overstrand_goes_OVER ? "true" : "false"));

       return false;

     }

     /************* overstrand_orientation **************/

    printf("checking computed overstrand_orientation...");

    if (computed_overstrand_orientation == overstrand_orientation) {

       pd_printf("pass (both are %OR)\n",NULL,&overstrand_orientation);

    } else {

       pd_printf("fail (expected %OR, got %OR)\n",NULL,
       &overstrand_orientation,&computed_overstrand_orientation);

       return false;

    }

  } else {

    printf("checking computed_tangle_slide_edges == NULL ...");

    if (computed_tangle_slide_edges != NULL) {

       printf("FAIL.\n");
       return false;

    }

    printf("pass\n");

    printf("checking computed_complementary_edges == NULL ...");

    if (computed_complementary_edges != NULL) {

       printf("FAIL.\n");
       return false;

    }

    printf("pass\n");

    printf("checking computed_complementary_or == NULL ...");

    if (computed_complementary_or != NULL) {

       printf("FAIL.\n");
       return false;

    }

    printf("pass\n");

  }
  
  printf("housecleaning...");
  
  pd_code_free(&pd);
  pd_tangle_free(&t);
  free(computed_tangle_slide_edges);
  free(computed_complementary_edges);
  free(computed_complementary_or);

  printf("done\n");

  printf("--------------------------------------------------------\n"
         "pdint_check_tslide_data_ok_and_find_te test Z : PASS\n"
	     "--------------------------------------------------------\n");
  
  return true;

}