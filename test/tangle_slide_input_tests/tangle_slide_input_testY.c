pd_code_t *pd_create_tangle_slide_input_testY_0() { 

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

