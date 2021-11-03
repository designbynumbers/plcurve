/* 

   Code for embedding pd_codes as planar graphs. This code first generates
   an augmented graph with only a single (spherical) embedding, then 
   uses Boyer's planarity library in order to compute a visibility representation
   of the augmented graph in the plane, then modifies the vrep to extract the 
   diagram embedding information. 

   We modify Boyer's ASCII art generator to produce diagrams on the terminal, 
   and also give several different output formats (an embedded graph in Mathematica
   format, for instance). 

   Jason Cantarella, October 2021.

*/

#ifndef __PD_EMBED_H__
#define __PD_EMBED_H__ 1

plCurve *plCurve_from_pd_code(pd_code_t *pdC);
/* Returns a plCurve representation of a pd_code_t which is guaranteed to reproduce
   the same pd_code when viewed from the z-axis. The vertices in the plCurve will 
   not be related to the crossings in pdC. */

char *graph_from_pd_code_Mathematica(pd_code_t *pdC);
/* Generates an embedded graph in Mathematica's style which has the same planar
   embedding type as pdC. The vertices in pdC will be the first vertices in the 
   graph, but additional vertices of degree 2 will be added to maintain the 
   graph drawing type. */

char *ascii_art_from_pd_code(pd_code_t *pdC);
/* Generates an ASCII art picture of a pdCode which is labeled with vertex and 
   edge numbers appropriately. */

#endif
