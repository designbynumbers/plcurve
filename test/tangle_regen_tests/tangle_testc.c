pd_code_t *pd_create_tangle_testc_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(4);
assert(pd != NULL);
pd->ncross = 4;
pd->nedges = 8;
pd->ncomps = 1;
pd->nfaces = 6;
sprintf(pd->hash,"%s","BAgGAwMDAwICAQgAAAAAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 2;
pd->cross[0].edge[2] = 7;
pd->cross[0].edge[3] = 3;
pd->cross[0].sign = 1;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 6;
pd->cross[1].edge[2] = 1;
pd->cross[1].edge[3] = 5;
pd->cross[1].sign = 0;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 4;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 5;
pd->cross[2].sign = 1;

pd->cross[3].edge[0] = 3;
pd->cross[3].edge[1] = 7;
pd->cross[3].edge[2] = 4;
pd->cross[3].edge[3] = 6;
pd->cross[3].sign = 1;


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
pd->edge[2].headpos = 1;
pd->edge[2].tail = 2;
pd->edge[2].tailpos = 2;

pd->edge[3].head = 3;
pd->edge[3].headpos = 0;
pd->edge[3].tail = 0;
pd->edge[3].tailpos = 3;

pd->edge[4].head = 2;
pd->edge[4].headpos = 1;
pd->edge[4].tail = 3;
pd->edge[4].tailpos = 2;

pd->edge[5].head = 1;
pd->edge[5].headpos = 3;
pd->edge[5].tail = 2;
pd->edge[5].tailpos = 3;

pd->edge[6].head = 3;
pd->edge[6].headpos = 3;
pd->edge[6].tail = 1;
pd->edge[6].tailpos = 1;

pd->edge[7].head = 0;
pd->edge[7].headpos = 2;
pd->edge[7].tail = 3;
pd->edge[7].tailpos = 1;


/* Component Data */

pd->comp[0].nedges = 8;
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


/* Face data */

pd->face[0].nedges = 3;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 0;

pd->face[0].edge[1] = 3;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 6;
pd->face[0].or[2] = 0;

pd->face[1].nedges = 3;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 0;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 5;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 2;
pd->face[1].or[2] = 1;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 1;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 6;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 4;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 2;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 4;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 7;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 2;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 5;
pd->face[4].or[1] = 1;

pd->face[5].nedges = 2;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 3;
pd->face[5].or[0] = 0;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 0;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

