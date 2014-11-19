/* 

   Code for deleting edges and crossings. This isn't part of 
   the exposed API for plCurve because it results in diagrams
   that don't pass pd_ok. 

   Jason Cantarella, November 2014.

*/

#ifndef __PD_DELETIONS_H__
#define __PD_DELETIONS_H__ 1

bool pd_edge_deleted(pd_code_t *pd, pd_idx_t edge); 
/* Returns true if edge has been deleted (by calling edge_delete) */

void pd_edge_delete(pd_code_t *pd, pd_idx_t edge);
/* 'Delete' edge from pd. 

   1) Set each crossing which refers to the edge to 
   refer to PD_UNSET_IDX, and 
  
   2) setting each reference to edge in the faces to 
   PD_UNSET_IDX, and

   3) set each reference to edge in the component 
   records to PD_UNSET_IDX, and

   4) set the head and tail of the edge record to 
   PD_UNSET_IDX, while setting the headpos and tailpos
   to PD_UNSET_POS.

   This does not compact the buffer of edges, and probably
   results in a pd_code_t which doesn't pass pd_ok.

*/

#endif
