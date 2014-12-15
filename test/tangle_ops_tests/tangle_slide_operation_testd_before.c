pd_code_t *pd_create_tangle_slide_operation_testd_before_0() { 

/* This procedure is machine generated by pd_write_c */
/* and probably shouldn't be hand-edited. */

pd_code_t *pd;
pd = pd_code_new(6);
assert(pd != NULL);
pd->ncross = 6;
pd->nedges = 12;
pd->ncomps = 2;
pd->nfaces = 8;
sprintf(pd->hash,"%s","BgwIBQQDAwMDAgECBgYAAAAAAAAAAAA");

/* Crossing data. */

pd->cross[0].edge[0] = 0;
pd->cross[0].edge[1] = 3;
pd->cross[0].edge[2] = 1;
pd->cross[0].edge[3] = 4;
pd->cross[0].sign = 0;

pd->cross[1].edge[0] = 0;
pd->cross[1].edge[1] = 7;
pd->cross[1].edge[2] = 5;
pd->cross[1].edge[3] = 8;
pd->cross[1].sign = 1;

pd->cross[2].edge[0] = 1;
pd->cross[2].edge[1] = 9;
pd->cross[2].edge[2] = 2;
pd->cross[2].edge[3] = 10;
pd->cross[2].sign = 0;

pd->cross[3].edge[0] = 2;
pd->cross[3].edge[1] = 9;
pd->cross[3].edge[2] = 3;
pd->cross[3].edge[3] = 8;
pd->cross[3].sign = 1;

pd->cross[4].edge[0] = 4;
pd->cross[4].edge[1] = 10;
pd->cross[4].edge[2] = 5;
pd->cross[4].edge[3] = 11;
pd->cross[4].sign = 0;

pd->cross[5].edge[0] = 6;
pd->cross[5].edge[1] = 6;
pd->cross[5].edge[2] = 11;
pd->cross[5].edge[3] = 7;
pd->cross[5].sign = 1;


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

pd->edge[3].head = 0;
pd->edge[3].headpos = 1;
pd->edge[3].tail = 3;
pd->edge[3].tailpos = 2;

pd->edge[4].head = 4;
pd->edge[4].headpos = 0;
pd->edge[4].tail = 0;
pd->edge[4].tailpos = 3;

pd->edge[5].head = 1;
pd->edge[5].headpos = 2;
pd->edge[5].tail = 4;
pd->edge[5].tailpos = 2;

pd->edge[6].head = 5;
pd->edge[6].headpos = 1;
pd->edge[6].tail = 5;
pd->edge[6].tailpos = 0;

pd->edge[7].head = 1;
pd->edge[7].headpos = 1;
pd->edge[7].tail = 5;
pd->edge[7].tailpos = 3;

pd->edge[8].head = 3;
pd->edge[8].headpos = 3;
pd->edge[8].tail = 1;
pd->edge[8].tailpos = 3;

pd->edge[9].head = 2;
pd->edge[9].headpos = 1;
pd->edge[9].tail = 3;
pd->edge[9].tailpos = 1;

pd->edge[10].head = 4;
pd->edge[10].headpos = 1;
pd->edge[10].tail = 2;
pd->edge[10].tailpos = 3;

pd->edge[11].head = 5;
pd->edge[11].headpos = 2;
pd->edge[11].tail = 4;
pd->edge[11].tailpos = 3;


/* Component Data */

pd->comp[0].nedges = 6;
pd->comp[0].tag = 'A';

pd->comp[0].edge = calloc(pd->comp[0].nedges,sizeof(pd_idx_t));
assert(pd->comp[0].edge != NULL);

pd->comp[0].edge[0] = 0;
pd->comp[0].edge[1] = 1;
pd->comp[0].edge[2] = 2;
pd->comp[0].edge[3] = 3;
pd->comp[0].edge[4] = 4;
pd->comp[0].edge[5] = 5;

pd->comp[1].nedges = 6;
pd->comp[1].tag = 'B';

pd->comp[1].edge = calloc(pd->comp[1].nedges,sizeof(pd_idx_t));
assert(pd->comp[1].edge != NULL);

pd->comp[1].edge[0] = 6;
pd->comp[1].edge[1] = 7;
pd->comp[1].edge[2] = 8;
pd->comp[1].edge[3] = 9;
pd->comp[1].edge[4] = 10;
pd->comp[1].edge[5] = 11;


/* Face data */

pd->face[0].nedges = 5;
pd->face[0].edge = calloc(pd->face[0].nedges,sizeof(pd_idx_t));
pd->face[0].or = calloc(pd->face[0].nedges,sizeof(pd_or_t));
assert(pd->face[0].edge != NULL);
assert(pd->face[0].or != NULL);

pd->face[0].edge[0] = 0;
pd->face[0].or[0] = 1;

pd->face[0].edge[1] = 4;
pd->face[0].or[1] = 1;

pd->face[0].edge[2] = 11;
pd->face[0].or[2] = 1;

pd->face[0].edge[3] = 6;
pd->face[0].or[3] = 0;

pd->face[0].edge[4] = 7;
pd->face[0].or[4] = 1;

pd->face[1].nedges = 4;
pd->face[1].edge = calloc(pd->face[1].nedges,sizeof(pd_idx_t));
pd->face[1].or = calloc(pd->face[1].nedges,sizeof(pd_or_t));
assert(pd->face[1].edge != NULL);
assert(pd->face[1].or != NULL);

pd->face[1].edge[0] = 2;
pd->face[1].or[0] = 1;

pd->face[1].edge[1] = 8;
pd->face[1].or[1] = 0;

pd->face[1].edge[2] = 5;
pd->face[1].or[2] = 0;

pd->face[1].edge[3] = 10;
pd->face[1].or[3] = 0;

pd->face[2].nedges = 3;
pd->face[2].edge = calloc(pd->face[2].nedges,sizeof(pd_idx_t));
pd->face[2].or = calloc(pd->face[2].nedges,sizeof(pd_or_t));
assert(pd->face[2].edge != NULL);
assert(pd->face[2].or != NULL);

pd->face[2].edge[0] = 0;
pd->face[2].or[0] = 0;

pd->face[2].edge[1] = 8;
pd->face[2].or[1] = 1;

pd->face[2].edge[2] = 3;
pd->face[2].or[2] = 1;

pd->face[3].nedges = 3;
pd->face[3].edge = calloc(pd->face[3].nedges,sizeof(pd_idx_t));
pd->face[3].or = calloc(pd->face[3].nedges,sizeof(pd_or_t));
assert(pd->face[3].edge != NULL);
assert(pd->face[3].or != NULL);

pd->face[3].edge[0] = 1;
pd->face[3].or[0] = 0;

pd->face[3].edge[1] = 3;
pd->face[3].or[1] = 0;

pd->face[3].edge[2] = 9;
pd->face[3].or[2] = 1;

pd->face[4].nedges = 3;
pd->face[4].edge = calloc(pd->face[4].nedges,sizeof(pd_idx_t));
pd->face[4].or = calloc(pd->face[4].nedges,sizeof(pd_or_t));
assert(pd->face[4].edge != NULL);
assert(pd->face[4].or != NULL);

pd->face[4].edge[0] = 1;
pd->face[4].or[0] = 1;

pd->face[4].edge[1] = 10;
pd->face[4].or[1] = 1;

pd->face[4].edge[2] = 4;
pd->face[4].or[2] = 0;

pd->face[5].nedges = 3;
pd->face[5].edge = calloc(pd->face[5].nedges,sizeof(pd_idx_t));
pd->face[5].or = calloc(pd->face[5].nedges,sizeof(pd_or_t));
assert(pd->face[5].edge != NULL);
assert(pd->face[5].or != NULL);

pd->face[5].edge[0] = 5;
pd->face[5].or[0] = 1;

pd->face[5].edge[1] = 7;
pd->face[5].or[1] = 0;

pd->face[5].edge[2] = 11;
pd->face[5].or[2] = 0;

pd->face[6].nedges = 2;
pd->face[6].edge = calloc(pd->face[6].nedges,sizeof(pd_idx_t));
pd->face[6].or = calloc(pd->face[6].nedges,sizeof(pd_or_t));
assert(pd->face[6].edge != NULL);
assert(pd->face[6].or != NULL);

pd->face[6].edge[0] = 2;
pd->face[6].or[0] = 0;

pd->face[6].edge[1] = 9;
pd->face[6].or[1] = 0;

pd->face[7].nedges = 1;
pd->face[7].edge = calloc(pd->face[7].nedges,sizeof(pd_idx_t));
pd->face[7].or = calloc(pd->face[7].nedges,sizeof(pd_or_t));
assert(pd->face[7].edge != NULL);
assert(pd->face[7].or != NULL);

pd->face[7].edge[0] = 6;
pd->face[7].or[0] = 1;


/* End of data. */

assert(pd_ok(pd));
return pd;

}

