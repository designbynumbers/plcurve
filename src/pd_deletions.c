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

#include<plcTopology.h>
//#include<libcassie/cassie.h>
//#include"/usr/local/include/thrift/Thrift.h"
//#include<python2.7/Python.h>
#include<pd_multidx.h>
#include<pd_dihedral.h>
#include<pd_perm.h>
#include<pd_orientation.h>
#include<pd_isomorphisms.h>
#include<pd_sortedbuf.h>
#include<pd_deletions.h>

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

  pd->cross[pd->edges[edge].head].edge[pd->edges[edge].headpos] = PD_UNSET_IDX;
  pd->cross[pd->edges[edge].tail].edge[pd->edges[edge].tailpos] = PD_UNSET_IDX;
  pd->cross[pd->edges[edge].head].sign = PD_UNSET_ORIENTATION;
  pd->cross[pd->edges[edge].tail].sign = PD_UNSET_ORIENTATION;
  
  pd_idx_t i,j;

  for(i=0;i<pd->nfaces;i++) { 

    for(j=0;j<pd->face[i].nedges;j++) { 

      if (pd->face[i].edge[j] == edge) { 
	
	pd->face[i].edge[j] = PD_UNSET_IDX;
	pd->face[i].or[j] = PD_UNSET_ORIENTATION;

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
  
