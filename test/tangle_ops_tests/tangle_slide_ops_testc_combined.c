/* tangle slide operations c assembled c code */

pd_code_t *pd_create_tangle_slide_operation_testc_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 3;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGBAQCAgICAwQCAgAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 4;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 5;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 6;
pd->cross[1].edge[2] = 3;
pd->cross[1].edge[3] = 7;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 6;
pd->cross[3].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 1;
pd->edge[3].headpos = 2;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 1;

pd->edge[5].head = 0;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 3;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 1;
pd->edge[7].headpos = 3;
pd->edge[7].tail = 3;
pd->edge[7].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 4;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 4;
pd->comp[1].edge[1] = 5;

pd->comp[2].nedges = 2;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 6;
pd->comp[2].edge[1] = 7;


/* Face data */

pd->face[0].nedges = 4;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 5;
pd->face[0].or[1] = 0;

pd->face[0].edge[2] = 2;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 0;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;

pd->face[1].edge[1] = 7;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 2;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 4;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 4;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 5;
pd->face[3].or[1] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 3;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 6;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_tangle_slide_operation_testc_after_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 2;
pd->nedges = 4;
pd->ncomps = 2;
pd->nfaces = 4;
sprintf(pd->hash,"%s","AgQEAgICAgICAgAAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 2;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 3;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 2;
pd->cross[1].sign = 0;


/* Edge data */

pd->edge[0].head = 0;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 1;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 1;
pd->edge[1].headpos = 2;
pd->edge[1].tail = 0;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 1;
pd->edge[2].headpos = 3;
pd->edge[2].tail = 0;
pd->edge[2].tailpos = 1;

pd->edge[3].head = 0;
pd->edge[3].headpos = 3;
pd->edge[3].tail = 1;
pd->edge[3].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 2;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;

pd->comp[1].nedges = 2;
pd->comp[1].tag = 'C';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 2;
pd->comp[1].edge[1] = 3;


/* Face data */

pd->face[0].nedges = 2;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 2;
pd->face[0].or[1] = 0;

pd->face[1].nedges = 2;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 3;
pd->face[1].or[1] = 0;

pd->face[2].nedges = 2;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 2;
pd->face[2].or[1] = 1;

pd->face[3].nedges = 2;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 1;

pd->face[3].edge[1] = 3;
pd->face[3].or[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

pd_code_t *pd_create_tangle_slide_operation_testc_after_1() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(2);
assert(pd != NULL);
pd->ncross = 0;
pd->nedges = 1;
pd->ncomps = 1;
pd->nfaces = 2;
sprintf(pd->hash,"%s","AAECAQEBAQAAAAAAAAAAAAAAAAAAAAA");

/* Crossing data. */


/* Edge data */

pd->edge[0].head = -1;
pd->edge[0].headpos = 4;
pd->edge[0].tail = -1;
pd->edge[0].tailpos = 4;


/* Component Data */

pd->comp[0].nedges = 1;
pd->comp[0].tag = 'B';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;


/* Face data */

pd->face[0].nedges = 1;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[1].nedges = 1;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}


bool pdint_check_tangle_slide_ops_testc() {

  /*                                                          */
/*              /----C--7---<-----\          (1)            */
/*           c1 |       (5)       | c3                      */
/*      /-----------<---3---<----------2-\                  */
/*      |       |                 |      |                  */
/*      |   T   6                 6      |                  */
/*      |       v                 ^      |                  */
/*      v   +------------------------+   ^                  */
/*      |   |   |       (4)       |  |   0                  */
/*      |   |   \-->--C--6--->----/  |   A                  */
/*     e0   |     (0)          (0)   |   0                  */
/*      |   |                        |   2                  */
/*      |   |   /----<-B-5---<----\  |   0                  */
/*      A   |   |                 |  |   |                  */
/*      |   |   5       (3)       5  |   ^                  */
/*      v   |   v                 ^  |   0                  */
/*      0   +------------------------+   0                  */
/*      0    c0 |                 | c2   |                  */
/*      \--0--000000->--0-1--00->-0000-2-/                  */
/*              |                 |                         */
/*              4       (2)       4         (1)             */
/*              \-->-B---4---->---/                         */
/*                                                          */
 						
pd_idx_t nedges = 4;	      			
pd_idx_t tangle_faces[4] = {4,0,3,0};	
pd_idx_t tangle_edges[4] = {6,5,5,6};	

pd_idx_t noverstrand_edges = 3;
pd_idx_t overstrand_edges[3] = {0,1,2};
pd_idx_t border_faces[3] = {0,3,0};

/*                                                          */
/*      /-----<--A-----0----\                               */
/*      |                   0                               */
/*      |                   |                               */
/*      |                   |c1-                             */
/*      0         /----3---------<-----\                    */
/*      |         3         |          |                    */
/*      |         v         1          |                    */
/*      |       c0|+        |          |                    */
/*      \---0->---|----1->--/          2                    */
/*                |                    |                    */
/*                2                    ^       /-----B--\   */
/*                C                    |       |        |   */
/*                |                    |       |        B   */
/*                \---->-----2---------/       B        |   */
/*                                             |        |   */
/*                                             \---B----/   */
/*                                                          */


  printf("--------------------------------------------------\n"
    	 "pdint check tangle slide ops test c\n"
         "--------------------------------------------------\n");

  printf("creating 'before' pd...");
  pd_code_t *pd = pd_create_tangle_slide_operation_testc_before_0();
  if (!pd_ok(pd)) {

    printf("fail (doesn't pass pd_ok)\n");
    return false;

  }
  
  printf("pass (passes pd_ok)\n");

  printf("parsing tangle information...");
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
  printf("parsing overstrand_edges and border_faces...");
  printf("%d overstrand edges\n",noverstrand_edges);

  printf("running pd_tangle_slide...");

  pd_code_t **children;
  pd_idx_t nchildren;

  pd_tangle_slide(pd,t,
                  noverstrand_edges,overstrand_edges,border_faces,
                  &nchildren,&children);

  printf("pass (didn't crash, returned %d children)\n",nchildren);
  
  printf("checking children for pd_ok...\n");
  
  for(i=0;i<nchildren;i++) {

     if (!pd_ok(children[i])) {

        pd_printf("child %d fails pd_ok. \n %PD",children[i],i);
        return false;

     } else {

        printf("\t%d crossing child pd #%d is pd_ok...pass\n",
               children[i]->ncross,i);

     }

   }
  
   pd_idx_t nxchildren = 2;
   pd_code_t *xchildren[2];
   xchildren[0] = pd_create_tangle_slide_operation_testc_after_0();
   xchildren[1] = pd_create_tangle_slide_operation_testc_after_1();


   printf("checking %d expected children for pd_ok...\n",nxchildren);

   for(i=0;i<nxchildren;i++) {

    if (!pd_ok(xchildren[i])) {

       pd_printf("expected child %d failed pd_ok:\n %PD",xchildren[i],i);
       return false;

    } else {

       printf("\t%d crossing expected child pd #%d is pd_ok...pass\n",
              xchildren[i]->ncross,i);

    }

   }

   printf("comparing list of xchildren against actual children...\n");

   if (!compare_list_of_pds(nchildren,children,nxchildren,xchildren)) {

       printf("FAIL. Couldn't match expected and actual children\n"
              "after tangle_slide\n");
       return false;

   } else {

       printf("comparing list of xchildren against actual children...pass\n");

    }

  printf("housecleaning...\n");

  printf("\tfreeing parent pd and tangle...");
  pd_code_free(&pd);
  pd_tangle_free(&t);
  printf("done\n");

  printf("\tfreeing list of children and children buffer...");
  for(i=0;i<nchildren;i++) {
     pd_code_free(&(children[i]));
  }
  free(children);
  printf("done\n");

  printf("\tfreeing list of xchildren...");
  for(i=0;i<nxchildren;i++) {
     pd_code_free(&(xchildren[i]));
  }
  printf("done\n");

  printf("--------------------------------------------------------\n"
         "pdint check tangle slide ops testc : PASS\n"
	     "--------------------------------------------------------\n");
  
  return true;

}