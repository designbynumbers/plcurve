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

void pd_crossing_delete(pd_code_t *pd,pd_idx_t cross);

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

#endif
