/* 

   Code for deleting edges and crossings. This isn't part of 
   the exposed API for plCurve because it results in diagrams
   that don't pass pd_ok. 

   Jason Cantarella, November 2014.

*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif 

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_STDARG_H
   #include<stdarg.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#ifdef HAVE_CTYPE_H
   #include<ctype.h>
#endif

#ifdef HAVE_MATH_H
  #include<math.h>
#endif

#include"plcTopology.h"
//#include<libcassie/cassie.h>
//#include"/usr/local/include/thrift/Thrift.h"
//#include<python2.7/Python.h>
#include"pd_multidx.h"
#include"pd_dihedral.h"
#include"pd_perm.h"
#include"pd_orientation.h"
#include"pd_isomorphisms.h"
#include"pd_sortedbuf.h"
#include"pd_deletions.h"

bool pd_edge_deleted(pd_code_t *pd, pd_idx_t edge)
/* Returns true if edge has been deleted (by calling edge_delete) */
{

  return (pd->edge[edge].head == PD_UNSET_IDX &&
	  pd->edge[edge].tail == PD_UNSET_IDX &&
	  pd->edge[edge].headpos == PD_UNSET_POS &&
	  pd->edge[edge].tailpos == PD_UNSET_POS);
}


void pd_edge_delete(pd_code_t *pd, pd_idx_t edge)
/* 'Delete' edge from pd. 

   1) Set each crossing which refers to the edge to 
   refer to PD_UNSET_IDX, and 
  
   2) setting each reference to edge in the faces to 
   PD_UNSET_IDX, and

   3) set the head and tail of the edge record to 
   PD_UNSET_IDX, while setting the headpos and tailpos
   to PD_UNSET_POS.

   This does not compact the buffer of edges, and probably
   results in a pd_code_t which doesn't pass pd_ok.

*/
{
  pd_check_edge(SRCLOC,pd,edge);

  pd->cross[pd->edge[edge].head].edge[pd->edge[edge].headpos] = PD_UNSET_IDX;
  pd->cross[pd->edge[edge].tail].edge[pd->edge[edge].tailpos] = PD_UNSET_IDX;
  pd->cross[pd->edge[edge].head].sign = PD_UNSET_ORIENTATION;
  pd->cross[pd->edge[edge].tail].sign = PD_UNSET_ORIENTATION;
  
  pd_idx_t i,j;

  for(i=0;i<pd->nfaces;i++) { 

    for(j=0;j<pd->face[i].nedges;j++) { 

      if (pd->face[i].edge[j] == edge) { 
	
	pd->face[i].edge[j] = PD_UNSET_IDX;
	pd->face[i].orient[j] = PD_UNSET_ORIENTATION;

      }

    }

  }

  for(i=0;i<pd->ncomps;i++) { 

    for(j=0;j<pd->comp[i].nedges;j++) { 

      if (pd->comp[i].edge[j] == edge) { 

	pd->comp[i].edge[j] = PD_UNSET_IDX;

      }

    }

  }

  pd->edge[edge].head = PD_UNSET_IDX;
  pd->edge[edge].headpos = PD_UNSET_POS;
  pd->edge[edge].tail = PD_UNSET_IDX;
  pd->edge[edge].tailpos = PD_UNSET_POS;

}
  
void pd_crossing_delete(pd_code_t *pd,pd_idx_t cross)

/* "Deletes" the crossing "cross" from the pd code. 
    In practice, this means that the following actions
    are taken:

    1) Look at the edges coming into the crossing. If 
       a pair of opposite positions 0-2 or 1-3 has both 
       edges set, splice the edges together. For example,
       given

                0
                ^
             1--+--3
                | 
                2

       and alter the tail of the edge in position 0 to be equal to
       the tail of the edge in position 2 and we'll delete the edge
       in position 2. 

       Set instances of the deleted edges in the component records
       to PD_UNSET_IDX to record that they've been deleted. 

       If a pair of opposite positions has only one edge
       set, don't do any edge deletions, but change the head 
       or tail of the corresponding edge to PD_UNSET_IDX.

    2) Set the data associated to the crossing (edges, sign) to 
       PD_UNSET_IDX and PD_UNSET_OR, as appropriate. 

*/

{
  pd_check_cr(SRCLOC,pd,cross);

  pd_idx_t posbase;

  for(posbase=0;posbase<2;posbase++) { 

    if (pd->cross[cross].edge[posbase] != PD_UNSET_IDX &&
	pd->cross[cross].edge[posbase+2] != PD_UNSET_IDX) {

      /* Both edges in a pair are SET; we'll delete the one coming IN. */

      pd_idx_t in_edge, out_edge;

      if (pd->edge[pd->cross[cross].edge[posbase]].head == cross &&
	  pd->edge[pd->cross[cross].edge[posbase+2]].tail == cross) { 

	in_edge = pd->cross[cross].edge[posbase];
	out_edge = pd->cross[cross].edge[posbase+2];
	  
      } else if (pd->edge[pd->cross[cross].edge[posbase+2]].head == cross &&
		 pd->edge[pd->cross[cross].edge[posbase]].tail == cross) { 

	in_edge = pd->cross[cross].edge[posbase+2];
	out_edge = pd->cross[cross].edge[posbase];

      } else {  /* Add a safeguard in case things aren't right. */

	pd_error(SRCLOC,"trying to delete %CROSS, and found edges set at positions\n"
		 "%d and %d, but those edge records %EDGE and %EDGE don't pass through\n"
		 "this crossing correctly\n",pd,
		 cross,posbase,posbase+2,pd->cross[cross].edge[posbase],
		 pd->cross[cross].edge[posbase+2]);
	exit(1);

      }
	
      /* Now we're ready to splice edges together. */

      pd->edge[out_edge].tail = pd->edge[in_edge].tail;
      pd->edge[out_edge].tailpos = pd->edge[in_edge].tailpos;

      /* Purge references to the deleted edge from crossing records. */

      pd->cross[pd->edge[out_edge].tail].edge[pd->edge[out_edge].tailpos] = 
	out_edge;

      /* Erase data from the in_edge record in memory */

      pd->edge[in_edge].head = PD_UNSET_IDX;
      pd->edge[in_edge].headpos = PD_UNSET_POS;
      pd->edge[in_edge].tail = PD_UNSET_IDX;
      pd->edge[in_edge].tailpos = PD_UNSET_POS;

      /* Purge references to the deleted edge from the component records. */
      
      pd_idx_t i,j;

      for(i=0;i<pd->ncomps;i++) { 

	for(j=0;j<pd->comp[i].nedges;j++) { 

	  if (pd->comp[i].edge[j] == in_edge) { 

	    pd->comp[i].edge[j] = PD_UNSET_IDX;

	  }

	}

      }
  
    } else { /* At least one edge is unset. */

      pd_idx_t set_edge = PD_UNSET_IDX;

      if (pd->cross[cross].edge[posbase] != PD_UNSET_IDX) { 
	
	set_edge = pd->cross[cross].edge[posbase];

      } else if (pd->cross[cross].edge[posbase+2] != PD_UNSET_IDX) {

	set_edge = pd->cross[cross].edge[posbase+2];

      }

      if (set_edge != PD_UNSET_IDX) { /* There's exactly one set edge */

	if (pd->edge[set_edge].head == cross) { 

	  pd->edge[set_edge].head = PD_UNSET_IDX;
	  pd->edge[set_edge].headpos = PD_UNSET_POS;

	} else if (pd->edge[set_edge].tail == cross) { 

	  pd->edge[set_edge].tail = PD_UNSET_IDX;
	  pd->edge[set_edge].tailpos = PD_UNSET_POS;

	} 
	
      }

    }

  }

  /* Now we set the actual crossing data to UNSET */

  pd_idx_t i;

  for(i=0;i<4;i++) { 

    pd->cross[cross].edge[i] = PD_UNSET_IDX;
    pd->cross[cross].sign = PD_UNSET_ORIENTATION;

  }

}
