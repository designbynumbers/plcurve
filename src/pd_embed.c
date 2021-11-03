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

#include"plCurve.h"
#include"plcTopology.h"

/* In this version of the world, Judy.h is an installed header
   under libpdcode/, or a local header in plcurve/judy/src. */
#include<graph.h>

#include"pd_multidx.h"
#include"pd_dihedral.h"
#include"pd_perm.h"
#include"pd_orientation.h"
#include"pd_isomorphisms.h"

/* There are some unexposed primitives here which are going to
   be useful in all of the drawing: */

int augment_crossvert(pd_code_t *pdC,pd_idx_t vt)
/* The crossing vertices are the first pdC->ncross vertices in the augmented graph*/
/* However, the indexing in the graph library is 1-based, so we add one. */
{
  return vt+1;
}

  int augment_edgevert(pd_code_t *pdC,pd_idx_t edge,bool headQ)
/* The edge vertices in tail-head order by edge number are the next 2*pdC->nedges
   vertices in the augmented graph */
{
  return (pdC->ncross+1)+2*edge+(headQ ? 1:0)
}

int augment_facevert(pd_code_t *pdC,pd_idx_t face)
/* The face vertices in order by face number are the last pdC->nfaces vertices
   in the augmented graph */
{
  return (pdC->ncross+1)+(2*pdC->nedges)+(face)
}

int build_augmented_graph(pd_code_t *pdC,FILE *outfile)
/*
   Constructs a an adjacency list for the augmented graph from a pd_code.

   Let's call the edges of the pd_code "arcs" and stick with the name
   "crossings" for where the arcs intersect. We'll use "vertices" and
   "edges" for the augmented graph. 

   The algorithm is to replace each crossing of the pd_code with a 
   3x3 lattice subgraph of 9 vertices and 12 edges as below
                                          
             |                           
             b                            X 
             |                         v2--v1
    ---c-----v----a->                 X | X | X 
             |                         v3--v0
             d                            X  
             |                          
                                       
     We then replace each arc of the pd_code with a subgraph with 
     two vertices (joined by an edge) each with four "dangling edges"

                                         \  /
                                        --a2--
         v----a-----w       ->            |
                                        --a1--
                                         /  \

     The middle two dangling edges will be joined to corners of the
     lattice graph replacing each crossing, while the outer two will
     be joined to sides of the lattice graph replaceing the crossing,
     as below. 

             |                          |   |
             |                          b2--b1
             b                          | X |
             |                   ---c1--v2-v1--a2---
    ---c-----v-----a->              | X | X | X |
             |                   ---c2--v3-v0--a1---
             d                          | X |
             |                          d1-d2
             |                          |   |

    Proposition. The resulting graph is 3-connected. 

    Proof. By inspection, the only potential 2-vertex disconnections
    of this subgraph delete a pair of vertices replacing an arc of the
    pd_code (e.g. a1, a2) or delete a pair of side vertices on the
    lattice graph (e.g. v0, v1). 

    In either case, we either cut a single edge of the pd_code or
    disconnect an edge from a vertex. Neither operation can
    disconnect the pd_code as a whole. Q.E.D.

    The function will output the graph in adjacency list format (from
    the documentation in graphIO.c).

    The file format is

     On the first line    : N= number of vertices
     On N subsequent lines: #: a b c ... -1
     where # is a vertex number and a, b, c, ... are its neighbors.

     NOTE:  The vertex number is for file documentation only.  It is an
     error if the vertices are not in sorted order in the file.

     NOTE:  If a loop edge is found, it is ignored without error.
*/
{
  typedef struct adjlist_struct {

    int nadj;
    int *adj;

  } adjlist;

  adjlist *adj;

  adj = calloc(pdC->ncross + 2*pdC->nedges + pdC->nfaces);
  if (adj == NULL) {
    fprintf(stderr,
	    "build_augmented_graph: Couldn't allocate space for %d\n"
	    "adjacency lists for augmented graph from pdCode with \n"
	    "%d crossings, %d edges, %d faces\n",
	    pdC->ncross + 2*pdC->nedges + pdC->nfaces,
	    pdC->ncross,pdC->nedges,pdC->nfaces);
    exit(1);
  }

  /* We now figure out the degree of each vertex. The

  /* We now go ahead and start filling in edges with addEdge. We want to
     be somewhat careful about this, because we'll be keeping track of
     which edges are "real" edges when we delete them later. */





}
