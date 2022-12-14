/* tangle interior test a assembled c code */

pd_code_t *pd_create_tangle_testa_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(11);
assert(pd != NULL);
pd->ncross = 11;
pd->nedges = 22;
pd->ncomps = 4;
pd->nfaces = 13;
sprintf(pd->hash,"%s","0000000000pd0unset0hash00000000");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 5;
pd->cross[0].edge[2] = 9;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 21;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 18;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 16;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 15;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 16;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 17;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 3;
pd->cross[4].edge[1] = 21;
pd->cross[4].edge[2] = 4;
pd->cross[4].edge[3] = 20;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 5;
pd->cross[5].edge[1] = 18;
pd->cross[5].edge[2] = 6;
pd->cross[5].edge[3] = 19;
pd->cross[5].sign = 1;

pd->cross[6].edge[0] = 6;
pd->cross[6].edge[1] = 11;
pd->cross[6].edge[2] = 7;
pd->cross[6].edge[3] = 12;
pd->cross[6].sign = 1;

pd->cross[7].edge[0] = 7;
pd->cross[7].edge[1] = 13;
pd->cross[7].edge[2] = 8;
pd->cross[7].edge[3] = 12;
pd->cross[7].sign = 1;

pd->cross[8].edge[0] = 8;
pd->cross[8].edge[1] = 20;
pd->cross[8].edge[2] = 9;
pd->cross[8].edge[3] = 19;
pd->cross[8].sign = 1;

pd->cross[9].edge[0] = 10;
pd->cross[9].edge[1] = 14;
pd->cross[9].edge[2] = 11;
pd->cross[9].edge[3] = 15;
pd->cross[9].sign = 0;

pd->cross[10].edge[0] = 10;
pd->cross[10].edge[1] = 17;
pd->cross[10].edge[2] = 13;
pd->cross[10].edge[3] = 14;
pd->cross[10].sign = 1;


/* Edge data */

pd->edge[0].head = 1;
pd->edge[0].headpos = 0;
pd->edge[0].tail = 0;
pd->edge[0].tailpos = 0;

pd->edge[1].head = 2;
pd->edge[1].headpos = 0;
pd->edge[1].tail = 1;
pd->edge[1].tailpos = 2;

pd->edge[2].head = 3;
pd->edge[2].headpos = 0;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 4;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 0;
pd->edge[4].headpos = 3;
pd->edge[4].tail = 4;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 5;
pd->edge[5].headpos = 0;
pd->edge[5].tail = 0;
pd->edge[5].tailpos = 1;

pd->edge[6].head = 6;
pd->edge[6].headpos = 0;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 2;

pd->edge[7].head = 7;
pd->edge[7].headpos = 0;
pd->edge[7].tail = 6;
pd->edge[7].tailpos = 2;

pd->edge[8].head = 8;
pd->edge[8].headpos = 0;
pd->edge[8].tail = 7;
pd->edge[8].tailpos = 2;

pd->edge[9].head = 0;
pd->edge[9].headpos = 2;
pd->edge[9].tail = 8;
pd->edge[9].tailpos = 2;

pd->edge[10].head = 9;
pd->edge[10].headpos = 0;
pd->edge[10].tail = 10;
pd->edge[10].tailpos = 0;

pd->edge[11].head = 6;
pd->edge[11].headpos = 1;
pd->edge[11].tail = 9;
pd->edge[11].tailpos = 2;

pd->edge[12].head = 7;
pd->edge[12].headpos = 3;
pd->edge[12].tail = 6;
pd->edge[12].tailpos = 3;

pd->edge[13].head = 10;
pd->edge[13].headpos = 2;
pd->edge[13].tail = 7;
pd->edge[13].tailpos = 1;

pd->edge[14].head = 9;
pd->edge[14].headpos = 1;
pd->edge[14].tail = 10;
pd->edge[14].tailpos = 3;

pd->edge[15].head = 2;
pd->edge[15].headpos = 3;
pd->edge[15].tail = 9;
pd->edge[15].tailpos = 3;

pd->edge[16].head = 3;
pd->edge[16].headpos = 1;
pd->edge[16].tail = 2;
pd->edge[16].tailpos = 1;

pd->edge[17].head = 10;
pd->edge[17].headpos = 1;
pd->edge[17].tail = 3;
pd->edge[17].tailpos = 3;

pd->edge[18].head = 5;
pd->edge[18].headpos = 1;
pd->edge[18].tail = 1;
pd->edge[18].tailpos = 3;

pd->edge[19].head = 8;
pd->edge[19].headpos = 3;
pd->edge[19].tail = 5;
pd->edge[19].tailpos = 3;

pd->edge[20].head = 4;
pd->edge[20].headpos = 3;
pd->edge[20].tail = 8;
pd->edge[20].tailpos = 1;

pd->edge[21].head = 1;
pd->edge[21].headpos = 1;
pd->edge[21].tail = 4;
pd->edge[21].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 10;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;
pd->comp[0].edge[6] = 6;
pd->comp[0].edge[7] = 7;
pd->comp[0].edge[8] = 8;
pd->comp[0].edge[9] = 9;

pd->comp[1].nedges = 4;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 10;
pd->comp[1].edge[1] = 11;
pd->comp[1].edge[2] = 12;
pd->comp[1].edge[3] = 13;

pd->comp[2].nedges = 4;
pd->comp[2].tag = 'C';

pd->comp[2].edge = calloc(pd->comp[2].nedges,sizeof(pd_idx_t));
assert(pd->comp[2].edge != NULL);

pd->comp[2].edge[0] = 14;
pd->comp[2].edge[1] = 15;
pd->comp[2].edge[2] = 16;
pd->comp[2].edge[3] = 17;

pd->comp[3].nedges = 4;
pd->comp[3].tag = 'D';

pd->comp[3].edge = calloc(pd->comp[3].nedges,sizeof(pd_idx_t));
assert(pd->comp[3].edge != NULL);

pd->comp[3].edge[0] = 18;
pd->comp[3].edge[1] = 19;
pd->comp[3].edge[2] = 20;
pd->comp[3].edge[3] = 21;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].orient = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].orient != NULL);

pd->face[0].edge[0] = 1;
pd->face[0].orient[0] = 1;

pd->face[0].edge[1] = 15;
pd->face[0].orient[1] = 0;

pd->face[0].edge[2] = 11;
pd->face[0].orient[2] = 1;

pd->face[0].edge[3] = 6;
pd->face[0].orient[3] = 0;

pd->face[0].edge[4] = 18;
pd->face[0].orient[4] = 0;

pd->face[1].nedges = 5;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].orient = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].orient != NULL);

pd->face[1].edge[0] = 3;
pd->face[1].orient[0] = 1;

pd->face[1].edge[1] = 20;
pd->face[1].orient[1] = 0;

pd->face[1].edge[2] = 8;
pd->face[1].orient[2] = 0;

pd->face[1].edge[3] = 13;
pd->face[1].orient[3] = 1;

pd->face[1].edge[4] = 17;
pd->face[1].orient[4] = 0;

pd->face[2].nedges = 4;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].orient = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].orient != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].orient[0] = 0;

pd->face[2].edge[1] = 21;
pd->face[2].orient[1] = 0;

pd->face[2].edge[2] = 3;
pd->face[2].orient[2] = 0;

pd->face[2].edge[3] = 16;
pd->face[2].orient[3] = 0;

pd->face[3].nedges = 4;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].orient = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].orient != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].orient[0] = 1;

pd->face[3].edge[1] = 17;
pd->face[3].orient[1] = 1;

pd->face[3].edge[2] = 10;
pd->face[3].orient[2] = 1;

pd->face[3].edge[3] = 15;
pd->face[3].orient[3] = 1;

pd->face[4].nedges = 4;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].orient = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].orient != NULL);

pd->face[4].edge[0] = 6;
pd->face[4].orient[0] = 1;

pd->face[4].edge[1] = 12;
pd->face[4].orient[1] = 1;

pd->face[4].edge[2] = 8;
pd->face[4].orient[2] = 1;

pd->face[4].edge[3] = 19;
pd->face[4].orient[3] = 0;

pd->face[5].nedges = 4;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].orient = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].orient != NULL);

pd->face[5].edge[0] = 7;
pd->face[5].orient[0] = 0;

pd->face[5].edge[1] = 11;
pd->face[5].orient[1] = 0;

pd->face[5].edge[2] = 14;
pd->face[5].orient[2] = 0;

pd->face[5].edge[3] = 13;
pd->face[5].orient[3] = 0;

pd->face[6].nedges = 3;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].orient = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].orient != NULL);

pd->face[6].edge[0] = 0;
pd->face[6].orient[0] = 0;

pd->face[6].edge[1] = 4;
pd->face[6].orient[1] = 0;

pd->face[6].edge[2] = 21;
pd->face[6].orient[2] = 1;

pd->face[7].nedges = 3;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].orient = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].orient != NULL);

pd->face[7].edge[0] = 0;
pd->face[7].orient[0] = 1;

pd->face[7].edge[1] = 18;
pd->face[7].orient[1] = 1;

pd->face[7].edge[2] = 5;
pd->face[7].orient[2] = 0;

pd->face[8].nedges = 3;
pd->face[8].edge = calloc(pd->face[8].nedges,sizeof(pd_idx_t));
pd->face[8].orient = calloc(pd->face[8].nedges,sizeof(pd_or_t));
assert(pd->face[8].edge != NULL);
assert(pd->face[8].orient != NULL);

pd->face[8].edge[0] = 4;
pd->face[8].orient[0] = 1;

pd->face[8].edge[1] = 9;
pd->face[8].orient[1] = 0;

pd->face[8].edge[2] = 20;
pd->face[8].orient[2] = 1;

pd->face[9].nedges = 3;
pd->face[9].edge = calloc(pd->face[9].nedges,sizeof(pd_idx_t));
pd->face[9].orient = calloc(pd->face[9].nedges,sizeof(pd_or_t));
assert(pd->face[9].edge != NULL);
assert(pd->face[9].orient != NULL);

pd->face[9].edge[0] = 5;
pd->face[9].orient[0] = 1;

pd->face[9].edge[1] = 19;
pd->face[9].orient[1] = 1;

pd->face[9].edge[2] = 9;
pd->face[9].orient[2] = 1;

pd->face[10].nedges = 2;
pd->face[10].edge = calloc(pd->face[10].nedges,sizeof(pd_idx_t));
pd->face[10].orient = calloc(pd->face[10].nedges,sizeof(pd_or_t));
assert(pd->face[10].edge != NULL);
assert(pd->face[10].orient != NULL);

pd->face[10].edge[0] = 2;
pd->face[10].orient[0] = 0;

pd->face[10].edge[1] = 16;
pd->face[10].orient[1] = 1;

pd->face[11].nedges = 2;
pd->face[11].edge = calloc(pd->face[11].nedges,sizeof(pd_idx_t));
pd->face[11].orient = calloc(pd->face[11].nedges,sizeof(pd_or_t));
assert(pd->face[11].edge != NULL);
assert(pd->face[11].orient != NULL);

pd->face[11].edge[0] = 7;
pd->face[11].orient[0] = 1;

pd->face[11].edge[1] = 12;
pd->face[11].orient[1] = 0;

pd->face[12].nedges = 2;
pd->face[12].edge = calloc(pd->face[12].nedges,sizeof(pd_idx_t));
pd->face[12].orient = calloc(pd->face[12].nedges,sizeof(pd_or_t));
assert(pd->face[12].edge != NULL);
assert(pd->face[12].orient != NULL);

pd->face[12].edge[0] = 10;
pd->face[12].orient[0] = 0;

pd->face[12].edge[1] = 14;
pd->face[12].orient[1] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}


bool tangle_testa() {

  /*                                                                            */
/*     ____13______^_______________    ________________________^_17_          */
/*    /                            \ /                              \         */
/*   /                              /10                              \        */
/*  /              (5)             / \               (3)              \       */
/* /                              14  10                                \      */
/* \                 +-------------+---+----<-----------+-+            /      */
/*  \   _____^____7__|_            v(12)v        ____2___v_|__         /       */
/*   \ /        (11) | \            \ /        / (10)     |  \       /        */
/*    \              |  \            \__15_^_________16_^_|___\ ____/         */
/*   /7\____^____12__|__ \  __^_11__/9       /2           |   3\              */
/*  /                |   6\                 ^             |     \             */
/*  \                v     ^      (0)      1              |     /             */
/*   \        (4)    |      6             /           (2) |    /              */
/*    \             _|______ \ ___<__18_______<__         ^   /               */
/*     \          19 |       5\         /1       \        |  3     (1)        */
/*      8         /  |         ^  (7)  ^          \       | /                 */
/*       \       /   |          \     0            \      |/                  */
/*        \     /    |           5   /             21     |                   */
/*         ^   v     | (9)        \ /                \    |                   */
/*          \ /      |             \        (6)       \ / |                   */
/*           /       |            /0\                  /  |<--tangle          */
/*          /8\      |           /   \                /4\ |                   */
/*         /   \     |          /     ^              /   \|                   */
/*        /     \____|___^_9___/      \_________4___/     |                   */
/*       /           |                                    |\                  */
/*      /            |          (8)                       | \                 */
/*      \            +-------------------->---------------+ /                 */
/*       \                                                 /                  */
/*        \_____________________20_______>________________/                   */

pd_idx_t nedges = 10;
pd_idx_t tangle_faces[10] = {8,1,2,10,3,12,5,11,4,9};
pd_idx_t tangle_edges[10]        = {20,3 ,16 ,2  ,10,14,7  , 12,19 ,9};
pd_boundary_or_t edge_bdy_or[10] = {in,in,out,out,in,in,out,out,out,in};
/*                                   0  1  2   3   4  5  6   7   8  9  */
pd_idx_t ninterior_cross = 7;
pd_idx_t interior_cross[7] = {0,1,2,4,5,6,9};

pd_idx_t ninterior_edges = 9;
pd_idx_t interior_edge[9] = {0,1,4,5,6,18,11,15,21};

pd_idx_t nstrands = 5;
pd_tangle_strand_t strand_data[5] = {{0,8,4,3},
				     {1,6,5,0},
				     {4,7,3,1},
				     {5,2,3,2},
				     {9,3,4,0}};


  printf("---------------------------------------\n"
	 "tangle_regenerate test a\n"
	 "---------------------------------------\n");

  printf("creating pd...");
  pd_code_t *pd = pd_create_tangle_testa_0();
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

  printf("testing edge boundary orientations...");

  for(i=0;i<t->nedges;i++) {

    if (t->edge_bdy_or[i] != edge_bdy_or[i]) {

      pd_idx_t j;

      printf("fail (list of boundary orientations for tangles doesn't match expected at pos %d)\n",i);

      printf("found:   ");

      for(j=0;j<t->nedges;j++) {

        if (t->edge_bdy_or[j] == in) { printf("  in "); }
        else { printf(" out "); }

      }

      printf("\n");

      printf("expected:");

      for(j=0;j<t->nedges;j++) {

        if (edge_bdy_or[j] == in) { printf("  in "); }
        else { printf(" out "); }

      }

      printf("\n");

      return false;

    }

  }

  printf("pass (edge boundary orientations match expected)\n");

  printf("checking interior crossings...");
  qsort(interior_cross,(size_t)(ninterior_cross),sizeof(pd_idx_t),pd_idx_cmp);

  if (t->ninterior_cross != ninterior_cross) {

    printf("fail (# interior crossings found (%d) != # expected (%d))\n",
    t->ninterior_cross,ninterior_cross);
    return false;

  }

  for(i=0;i<t->ninterior_cross;i++) {

     if (t->interior_cross[i] != interior_cross[i]) {

        pd_idx_t j;

        printf("fail (list of interior crossings doesn't match expected at pos %d)\n",i);

        printf("found:   ");

        for(j=0;j<t->ninterior_cross;j++) {

          printf(" %4d ",t->interior_cross[j]);

        }

        printf("\n");

        printf("expected:");

        for(j=0;j<t->ninterior_cross;j++) {

          printf(" %4d ",interior_cross[j]);

        }

        printf("\n");

        return false;

    }

  }

  printf("pass (interior crossings match)\n");


  printf("checking interior edges...");
  qsort(interior_edge,(size_t)(ninterior_edges),sizeof(pd_idx_t),pd_idx_cmp);

  if (t->ninterior_edges != ninterior_edges) {

    printf("fail (# interior edges found (%d) != # expected (%d))\n",
    t->ninterior_edges,ninterior_edges);
    return false;

  }

  for(i=0;i<t->ninterior_edges;i++) {

     if (t->interior_edge[i] != interior_edge[i]) {

        pd_idx_t j;

        printf("fail (list of interior edges doesn't match expected at pos %d)\n",i);

        printf("found:    ");

        for(j=0;j<t->ninterior_edges;j++) {

          printf(" %4d ",t->interior_edge[j]);

        }

        printf("\n");

        printf("expected:");

        for(j=0;j<t->ninterior_edges;j++) {

          printf(" %4d ",interior_edge[j]);

        }

        printf("\n");

        return false;

    }

  }

  printf("pass (interior edges match)\n");

  printf("testing strand data (%d) strands...\n",nstrands);

  for(i=0;i<nstrands;i++) {

     printf("\tstrand %d...",i);

     if(strand_data[i].start_edge != t->strand[i].start_edge) {

        printf("FAIL. (start_edge %d != expected start_edge %d)\n",
               strand_data[i].start_edge,t->strand[i].start_edge);
        return false;

     }

     if(strand_data[i].end_edge != t->strand[i].end_edge) {

        printf("FAIL. (end_edge %d != expected end_edge %d)\n",
               strand_data[i].end_edge,t->strand[i].end_edge);
        return false;

     }

     if(strand_data[i].nedges != t->strand[i].nedges) {

        printf("FAIL. (nedges %d != expected nedges %d)\n",
               strand_data[i].nedges,t->strand[i].nedges);
        return false;

     }

      if(strand_data[i].comp != t->strand[i].comp) {

        printf("FAIL. (comp %d != expected comp %d)\n",
               strand_data[i].comp,t->strand[i].comp);
        return false;

     }

     printf("pass\n");

  }

  printf("testing strand data (%d) strands...pass\n",nstrands);
  
  printf("housecleaning...");
  
  pd_code_free(&pd);
  pd_tangle_free(&t);

  printf("done\n");

  printf("---------------------------------------\n"
         "tangle test a : PASS\n"
	     "---------------------------------------\n");
  
  return true;

}