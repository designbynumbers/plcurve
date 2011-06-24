/* 

   plc_nearest_neighbor: This is a set of functions which include a reasonable primitive 
   nearest-neighbor query for plCurves. Since the data in a plCurve is so spatially correlated,
   it makes the most sense to use a skipping strategy.

*/

int GLOB_SEARCH_DIMENSION = 0;
plc_vector *GLOB_NN_BUFFER;

int exhaustive_nearest_neighbor(plc_vector query,int n,plc_vector *buffer,int stride)

/* This is a simple nearest neighbor loop, designed to be run on a subset of the buffer. */
/* This is a utility function which is NOT exposed in the interface, so it does not get */
/* the plc_ name convention. */

{
  int pos;
  double bestyet = DBL_MAX,this;
  int bestpos;

  for(pos=0;pos<n;pos+=stride) {

    this = plc_M_sq_dist(query,buffer[pos]);
    if (this < bestyet) { bestyet = this; bestpos = pos; }

  }

  return bestpos;

}

int max_dimension(int n, plc_vector *buffer) 

/* Finds the largest dimension of bbox of the buffer. Not exposed in the interface. */

{
  plc_vector bbox[2];
  bbox[0] = bbox[1] = buffer[0];
  int i,j;

  for(i=0;i<n;i++) {

    for(j=0;j<3;j++) {

      bbox[0].c[j] = (bbox[0].c[j] < buffer[i].c[j]) ? bbox[0].c[j] : buffer[i].c[j];
      bbox[1].c[j] = (bbox[1].c[j] > buffer[i].c[j]) ? bbox[1].c[j] : buffer[i].c[j];
      
    }

  }

  int bbdims[3];
  for(i=0;i<3;i++) { bbdims[i] = bbox[1].c[i] - bbox[0].c[i]; }

  if (bbdims[0] > bbdims[1]) {

    return (bbdims[0] > bbdims[2]) ? 0 : 2;

  } else {

    return (bbdims[1] > bbdims[2]) ? 1 : 2;

  }

}
   
int  nn_cmp (const int *i1, const int *i2)
{
  return (GLOB_NN_BUFFER[i1].c[GLOB_SEARCH_DIMENSION] > GLOB_NN_BUFFER[i2].c[GLOB_SEARCH_DIMENSION]) ? 1 : -1;
}     


int plc_nearest_neighbor(const plc_vector pt, const int n, const plc_vector *buffer,void **pc_data,int *plc_error)

/* This is an (at the moment) crude nearest neighbor code. Our model is that the
   code, whatever it ends up being, will give the nearest neighbor of the given 
   point among the points in the search_buffer.
   
   This involves some precomputation, so the pp_data buffer is expected to point
   to a void pointer set to NULL if we've never computed with this search_buffer 
   before and to a void pointer set non-NULL if we are using precomputed data.

   If we generate new precomputation data (and pc_data != NULL), we'll store the 
   precomputation data in *pp_data. It is the caller's responsibility to free this pointer 
   later when discarding the data. 

   If you don't want to preserve precomputation data, then pass a NULL to pp_data
   and this routine will free its own memory afterwards. 

   ---

   Model: We assume that the buffer is spatially coherent data which probably 
   comes from a polyline. So the method is to first do a brute-force search 
   among a small number of evenly spaced points in the buffer to get 
   some candidates. 

   We then sort the buffer by the maximum dimension of its bounding box.
   (Finding this dimension and sorting is the preprocessing step, and it 
   is this data which is saved). 

   Then we find the corresponding plane to the coordinate of pt using binary
   search and search linearly forward and backward from this point in the 
   sorted array, eliminating points as we find closer candidates.

*/
   
{
  plc_vector *check_pointer = NULL;
  int search_dimension;
  int *sorted_buffer = NULL;

  /* The first thing we do is check for precomputed data. */

  if (pc_data != NULL) {  /* We are using the pc data interface. */

    if (*pc_data != NULL) { /* We have passed in some pc data. */

      check_pointer = (plc_vector *)(**pc_data); /* The first thing in pc_data is supposed to be a pointer. */

      if (check_pointer != buffer) {

	(*plc_error) = PLC_E_STALE_DATA; /* The check pointer doesn't match. We will recompute the correct pcdata. */
	free(*pc_data);
	*pc_data = NULL;

      } else {

	/* The check pointer passed, so we assume that *pc_data is a pointer to valid precomputed data. */
	/* This data is supposed to look like a plc_vector check pointer followed by an int giving the  */
	/* sort dimension for the buffer followed by a buffer of ints giving the sort. To get at the    */
	/* later stuff requires a little pointer arithmetic. */

	int *raw_pcdata;

	raw_pcdata = (int *)((plc_vector *)(**pc_data) + 1); /* ?! We have to type the pointer to make the 
								arithmetic step "the size of a plc_vector *
								past **pcdata, then REtype the resulting 
								address so that now we step by ints. */

	search_dimension = raw_pcdata[0];
	sorted_buffer = raw_pcdata + 1;

      }
	
    }

  }

  /* If we had valid pc_data, and are using the interface, we loaded it. Let's check. */

  if (sorted_buffer == NULL) {  /* If we don't have a sort, we need to precompute. */

    search_dimension = max_dimension(n,buffer);
    GLOB_SEARCH_DIMENSION = search_dimension;
    GLOB_NN_BUFFER = buffer;

    sorted_buffer = malloc((n+1)*sizeof(int));
    for(i=1;i<=n;i++) { sorted_buffer[i] = i-1; }

    qsort(sorted_buffer+1,n,sizeof (int),nn_cmp);

  }

  /* We now have a sorted buffer which contains n in the 0th entry and a sorted set of 
     buffer indices in positions 1 to n. Our goal now is to find an entry i in the list 
     whose search_dimension coordinate is <= pt.c[search_dimension] while entry i+1 has
     search_dimension coordinate >= pt.c[search_dimension]. 

     This is basically bisection search on the list of coordinates. */

  int start; 

  if (pt.c[search_dimension] < buffer[sorted_buffer[1]].c[search_dimension]) {

    start = 0;

  } else if (pt.c[search_dimension] > buffer[sorted_buffer[n]].c[search_dimension]) {

    start = n-1;

  } else {

    int low = 1, high = n;
    int mid = (high-low)/2;  /* This will be an integer division */
    
    for(;high - low > 1;mid = (high-low)/2) {

      if (pt.c[search_dimension] > buffer[mid].c[search_dimension]) {

	low = mid;

      } else {

	high = mid;

      }

    }
    
    start = low;

  }

  /* We now have a point which has close to the same search_dimension coordinate as
     our search point. We will start our search for a nearest neighbor there. We know
     that as soon as the difference in search_dimension coordinates is larger than 
     the distance to our current winner, we can stop searching. */

  double best_dist = plc_distance(pt,buffer[sorted_buffer[start]]);
  int    best_index = sorted_buffer[start];
  int    up=start,down=start;

  for(;buffer[sorted_buffer[up]].c[search_dimension] - pt.c[search_dimension] < best_dist && up <= n;up++) {
    
    double this_dist = plc_distance(pt,buffer[sorted_buffer[up]]);

    if (this_dist < best_dist) {

      best_dist = this_dist;
      best_index = sorted_buffer[up];

    } 

  }

  for(;pt.c[search_dimension] - buffer[sorted_buffer[down]].c[search_dimension] < best_dist && down >= 1;down--) {
    
    double this_dist = plc_distance(pt,buffer[sorted_buffer[down]]);

    if (this_dist < best_dist) {

      best_dist = this_dist;
      best_index = sorted_buffer[down];

    } 

  }
  
  /* We now have a nearest neighbor in the buffer! */

  if (pc_data != NULL) {  /* We want to preserve precomputed data */

    if (*pc_data == NULL) { /* We weren't passed precomputed data (or we killed it and recomputed) */

      *pc_data = sorted_buffer;

    } 

  } else { /* We must do our own housekeeping */

    free(sorted_buffer);

  }

  return(best_index);

}

struct nearest_vertex_data {

  plCurve *check_curve;
  void **component_data;

}

void plc_nearest_vertex_pc_data_free(void **pc_data) 

/* This is long, because we have to be wary of double-frees or null frees 
   if this data has already been cleared at some other point in the program. */

{
  int cp;

  if (pc_data != NULL) {

    if (*pc_data != NULL) {

      if ((*pc_data)->check_curve != NULL) {

	if ((*pc_data)->component_data != NULL) {

	  for(cp=0;cp<(*pc_data)->check_curve->nc;cp++) {
	    
	    if ((*pc_data)->component_data[cp] != NULL) {
	      
	      free((*pc_data)->component_data[cp]);
	      (*pc_data)->component_data[cp] = NULL;
	      
	    }

	  }

	  free((*pc_data)->component_data);
	  (*pc_data)->component_data = NULL;

	}

      }

      free(*pc_data);
      *pc_data = NULL;

    }

  }

}
   
plc_vector plc_nearest_vertex(const plc_vector pt,const plCurve * const L,int *cp, int *vt, void **pc_data)

/* Finds a nearest vertex on L to pt using calls to plc_nearest_neighbor. */

{
  struct nearest_vertex_data *cpdata = NULL;

  if (pc_data != NULL) {  /* We are using the pc_data interface. */

    if (*pc_data != NULL) { /* We are passed pcdata. */

      cpdata = (struct nearest_vertex_data *)(*pc_data);

      if (cpdata->check_curve != L) { /* Check the curve pointer to provide some level of 
					 confidence that the data was precomputed for THIS curve. */

	plc_nearest_vertex_pc_data_free(pc_data);
	*pc_data = NULL;
	cpdata = NULL;

      }

    }

  }

  /* We loaded the data if we could. If we couldn't, allocate space for it now. */

  if (cpdata == NULL) {
    
    cpdata = calloc(1,sizeof(struct nearest_vertex_data));
    cpdata->check_curve = L;
    cpdata->component_data = calloc(L->nc,sizeof(void *));

  }

  /* Now we call nearest_neighbor on the various components, using the pc_data interface each time. */

  int cp;
  int *best_vt;
  plc_vector *best_vector;

  best_vt = calloc(L->nc,sizeof(int));
  best_vector = calloc(L->nc,sizeof(plc_vector));

  for(cp=0;cp<L->nc;cp++) {

    best_vt[cp] = plc_nearest_neighbor(pt,L->cp[cp].nv,L->cp[cp].vt,&(cpdata->component_data[cp]));
    best_vector[cp] = L->cp[cp].vt[best_vt[cp]];

  }

  /* We now need only figure out the closest of these points over the various components of L. */

  int best_cp; 
  best_cp = exhaustive_nearest_neighbor(pt,L->nc,best_vector,1);

  (*cp) = best_cp;
  (*vt) = best_vt[best_cp];

  /* Now we do some housecleaning. */

  free(best_vt);
  free(best_vector);

  if (pc_data != NULL) { /* We are using the precomputed data interface. */

    if (*pc_data == NULL) { /* We generated precomputation data. */

      *pc_data = (void *)(cpdata);

    }

  } else { /* We are responsible for our own cleanup, so kill the precomputation data. */

    plc_nearest_vertex_pc_data_free(&cpdata);

  }

  /* Now we return the final position. */

  return L->cp[*cp].vt[*vt];

}
    

  
  
    

    

      

	
  
    
