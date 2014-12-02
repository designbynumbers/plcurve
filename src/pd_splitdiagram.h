/* 

   Code for splitting diagrams into connected components.
   This isn't part of the exposed API because it operates 
   on diagrams that wouldn't pass pd_ok (on the other hand,
   it does _produce_ diagrams which should pass pd_ok).

   Jason Cantarella, November 2014.

*/

#ifndef __PD_SPLITDIAGRAM_H__
#define __PD_SPLITDIAGRAM_H__ 1

pd_idx_t pd_split_diagram(pd_code_t *pd,pd_code_t ***pd_children);

/* Splits a disconnected diagram with valid (though maybe disordered)
   components, edges, and crossings into a collection of child
   diagrams.  The child diagrams WILL have valid faces (and
   canonicalized ordering for everything), and should pass pd_ok. Some
   of the children may be 0-crossing unknots. Component tags are
   transferred by this operation. The case where there is only one
   child is not special-- we just return 1 and the buffer of pd_children
   contains a single pointer.

   The array pd_children will be allocated in this procedure (and is the 
   caller's responsibility to dispose of), and the number of children 
   is returned. 
*/

#endif
