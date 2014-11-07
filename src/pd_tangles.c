/* 

   Code for separating n-string tangles from the middle of 
   pd_codes and doing various operations on them. These include
   the "flype" and the "tangle_pass", which is a generalization
   of the RIII move (among others).

   Jason Cantarella, June 2014.

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

bool pd_xor(bool a, bool b) {

  return ((a && !b) || (b && !a));

}

pd_edge_t pd_unset_edge() /* Generates a special edge with everything
			      equal to the corresponding "UNSET" for use 
			      in edge deletion. */
{
  pd_edge_t unset_edge;

  unset_edge.head    = PD_UNSET_IDX;
  unset_edge.headpos = PD_UNSET_POS;
  unset_edge.tail    = PD_UNSET_IDX;
  unset_edge.tailpos = PD_UNSET_POS;

  return unset_edge;
}

pd_tangle_t *pd_tangle_new(pd_idx_t nedges) 
/* Creates space for a new tangle, including strands and orientations. */
{

  if (nedges == 0) { 

    pd_error(SRCLOC,"can't allocate tangle with 0 edges\n",NULL);
    exit(1);

  }

  if (nedges % 2 == 1) { 

    pd_error(SRCLOC,"can't allocate a tangle with an odd number "
	     "of boundary edges %d\n",NULL,nedges);
    exit(1);

  }

  pd_tangle_t *t = calloc(1,sizeof(pd_tangle_t));
  assert(t != NULL);

  t->nedges = nedges;
  t->nstrands = nedges/2;

  t->edge = calloc(t->nedges,sizeof(pd_idx_t));
  assert(t->edge != NULL);

  t->face = calloc(t->nedges,sizeof(pd_idx_t));
  assert(t->face != NULL);

  t->edge_bdy_or = calloc(t->nedges,sizeof(pd_boundary_or_t));
  assert(t->edge_bdy_or != NULL);

  t->strand = calloc(t->nstrands,sizeof(pd_tangle_strand_t));
  assert(t->strand != NULL);

  t->ninterior_cross = 0;
  t->interior_cross = NULL;
  
  t->ninterior_edges = 0;
  t->interior_edge = NULL;

  return t;

}

void pd_tangle_free(pd_tangle_t **T)
/* Carefully frees all memory associate with a tangle */
{
  if ((*T) == NULL) { return; } 
  pd_tangle_t *t = (*T);
  
  t->nedges = 0;

  if (t->edge != NULL) { 
  
    free(t->edge);
    t->edge = NULL;

  }

  if (t->face != NULL) { 

    free(t->face);
    t->face = NULL;

  }
  
  if (t->edge_bdy_or != NULL) { 

    free(t->edge_bdy_or);
    t->edge_bdy_or = NULL;

  }

  t->nstrands = 0;

  if (t->strand != NULL) { 

    free(t->strand);
    t->strand = NULL;

  }

  t->ninterior_cross = 0;

  if (t->interior_cross != NULL) {

    free(t->interior_cross);
    t->interior_cross = NULL;

  }

  t->ninterior_edges = 0;

  if (t->interior_edge != NULL) {

    free(t->interior_edge);
    t->interior_edge = NULL;

  }

  free(t);

}
  
bool pd_tangle_ok(pd_code_t *pd,pd_tangle_t *t) 
/* 
   Checks to see that the edge and face data in t
   defines a tangle within the pd code as in

               f[0]

          +-------------+     
          |             |     
    e[0]--+             +---e[n-1]
          |             |
     f[1] |             | f[n-1]
          |             |     
    e[1]--+             +---e[n-2]
          |             |     
    f[2]  +-------------+ f[n-2]    
             | ..... |
           e[2]     e[n-3]   
              
  and then consistency-checks remaining 
  tangle information.
*/
{
  pd_check_notnull(SRCLOC,"input pd code",pd);
  pd_check_notnull(SRCLOC,"input tangle",t);
  
  if (t->nedges == PD_UNSET_IDX) { 

    pd_error(SRCLOC,"tangle t->nedges is set to PD_UNSET_IDX\n",pd);
    return false;

  } 

  if (t->nstrands == PD_UNSET_IDX) { 

    pd_error(SRCLOC,"tangle t->nstrands is set to PD_UNSET_IDX (maybe tangle wasn't regenerated?)\n",pd);
    return false;

  }

  if (t->nstrands != t->nedges/2) { 

    pd_error(SRCLOC,"t->nstrands = %d is not equal to t->nedges = %d / 2 == %d.\n",
	     pd,t->nstrands,t->nedges,t->nedges/2);
    return false;

  }

  pd_check_notnull(SRCLOC,"t->edge",t->edge);
  pd_check_notnull(SRCLOC,"t->face",t->face);

  pd_idx_t *pos_face,pfp,*neg_face,nfp;
  pd_idx_t i,j;

  pos_face = calloc(t->nedges,sizeof(pd_idx_t));
  assert(pos_face != NULL);
  for(i=0;i<t->nedges;i++) { pos_face[i] = PD_UNSET_IDX; }

  neg_face = calloc(t->nedges,sizeof(pd_idx_t));
  assert(neg_face != NULL);
  for(i=0;i<t->nedges;i++) { neg_face[i] = PD_UNSET_IDX; }


  for(i=0;i<t->nedges;i++) { 

    pd_face_and_pos(pd,t->edge[i],&(pos_face[i]),&pfp,&(neg_face[i]),&nfp); 
    /* We keep overwriting the positions because we don't care about them*/ 

  }
  
  for(i=0;i<t->nedges;i++) { 

    if (!pd_xor(pos_face[i] == t->face[i],neg_face[i] == t->face[i])) { 
    
      pd_printf("edge t->edge[%d] = %EDGE is not on face t->face[%d] %FACE exactly once: tangle test fails in %PD\n",pd,i,t->edge[i],i,t->face[i]);
      return false; 

    }

    if (!pd_xor(pos_face[i] == t->face[(i+1)%t->nedges],neg_face[i] == t->face[(i+1)%t->nedges])) { 

      pd_printf("edge t->e[%d] = %EDGE is not on face t->face[%d] %FACE exactly once: tangle test fails in %PD\n",pd,i,t->edge[i],(i+1)%(t->nedges),t->face[(i+1)%t->nedges]);
      return false; 

    }
 
  }

  free(pos_face); pos_face = NULL;
  free(neg_face); neg_face = NULL;

  /* We now check edge_bdy_or */

  pd_check_notnull(SRCLOC,"t->edge_bdy_or",t->edge_bdy_or);
  pd_idx_t in_count = 0, out_count = 0;

  for(i=0;i<t->nedges;i++) { 
   
    if (t->edge_bdy_or[i] == in) { 
      
      if (in_count >= t->nstrands) { 

	pd_error(SRCLOC,"tangle has %d strands, but more than %d flow into the tangle\n",
		 pd,t->nstrands,in_count);
	return false;
      }

      in_count++; 

    } else if (t->edge_bdy_or[i] == out) {

      if (out_count >= t->nstrands) { 

	pd_error(SRCLOC,"tangle has %d strands, but more than %d flow out of the tangle\n",
		 pd,t->nstrands,out_count);
	return false;
      }
       
      out_count++; 

    } else { 

      pd_error(SRCLOC,"edge_bdy_or[%d] is set to %d, which is not in (== %d) or out (== %d)\n",
	       pd,i,in,out);

      return false;

    }

  }

  /* Ok, we've now checked that we have the right number of innies and outies. 
     Let's check the strand data itself for each strand. Remember that the start_edge
     and end_edge are now POSITIONS ON THE TANGLE.*/

  for(i=0;i<t->nstrands;i++) { 

    if (t->strand[i].start_edge >= t->nedges) { 
      
      pd_error(SRCLOC,"%TANGLE strand %d start edge has illegal index %d\n",pd,
	       t,i,t->strand[i].start_edge);
      exit(1);

    }

    if (t->strand[i].end_edge >= t->nedges) { 
      
      pd_error(SRCLOC,"%TANGLE strand %d end edge has illegal index %d\n",pd,
	       t,i,t->strand[i].end_edge);
      exit(1);

    }

    pd_check_edge(SRCLOC,pd,t->edge[t->strand[i].start_edge]);
    pd_check_edge(SRCLOC,pd,t->edge[t->strand[i].end_edge]);
    pd_check_cmp(SRCLOC,pd,t->strand[i].comp);

    if (t->strand[i].nedges == PD_UNSET_IDX) { 

      pd_error(SRCLOC,"strand[%d] of tangle has number of edges == PD_UNSET_IDX (maybe tangle not regenerated?)\n",
	       pd,i);
      return false;

    }

    if (t->strand[i].nedges > pd->comp[t->strand[i].comp].nedges) {

      pd_error(SRCLOC,"strand[%d] of tangle claims to have %d edges on %COMP"
	       "but comp only has %d edges\n",
	       pd,i,t->strand[i].nedges,t->strand[i].comp,pd->comp[t->strand[i].comp].nedges);
      return false;

    }

  }

  /* We're now going to try to test-walk each strand to verify that the
     number of edges and the start and end data is ok. Remember that the 
     start and end edge DO count in the total count of edges in a strand,
     so we should start with the number of found edges (j) equal to 1.*/

  for(i=0;i<t->nstrands;i++) { 

    pd_idx_t edge;
    for(edge = t->edge[t->strand[i].start_edge], j=1; 
	edge != t->edge[t->strand[i].end_edge] && j != t->strand[i].nedges;
	j++,edge = pd_next_edge(pd,edge)) { }

    if (edge != t->edge[t->strand[i].end_edge]) { 

      pd_error(SRCLOC,
	       "strand[%d] of tangle claims to have %d edges\n"
	       "but counting %d edges from start edge %EDGE\n"
	       "leaves us at %EDGE != the end_edge of %EDGE\n",
	       pd,i,t->strand[i].nedges,j,t->edge[t->strand[i].start_edge],
	       edge,t->edge[t->strand[i].end_edge]);

      return false;

    }

    if (j != t->strand[i].nedges) { 

      pd_error(SRCLOC,
	       "strand[%d] of tangle claims to have %d edges\n"
	       "but counting edges from start_edge %EDGE\n"
	       "finds the end_edge %EDGE after %d edges\n",
	       pd,i,t->strand[i].nedges,t->edge[t->strand[i].start_edge],
	       t->edge[t->strand[i].end_edge],j);

      return false;

    }

  }

  /* Now we check that comp makes sense. */

  for(i=0;i<t->nstrands;i++) { 

    pd_idx_t comp,pos;
    pd_component_and_pos(pd,t->edge[t->strand[i].start_edge],&comp,&pos);
    
    if (t->strand[i].comp != comp) { 

      pd_error(SRCLOC,
	       "tangle strand %d claims to be on component %COMP, but this\n"
	       "component does not include start_edge %EDGE\n",pd,i,
	       t->strand[i].comp,t->edge[t->strand[i].start_edge]);
      return false;

    }

    pd_component_and_pos(pd,t->edge[t->strand[i].end_edge],&comp,&pos);
    
    if (t->strand[i].comp != comp) { 

      pd_error(SRCLOC,
	       "tangle strand %d claims to be on component %COMP, but this\n"
	       "component does not include end_edge %EDGE\n",pd,i,
	       t->strand[i].comp,t->edge[t->strand[i].end_edge]);
      return false;

    }

  }
 
  /* There's no really practical way to check that the interior crossings
     are, in fact, interior (same with the edges) without actually regenerating
     the data here as part of the check, which seems like too much redundant 
     code. 

     We can check that they are allocated, however, and in fact that they refer 
     to valid crossing and edge numbers, and that they are ordered and unique, 
     and we do this in the spirit of checking for data corruption wherever it 
     might crop up. 

     We need to be careful here-- if there are no interior edges (or crossings)
     it is correct for the buffer to be NULL. */
     
  if (t->ninterior_cross != 0) {
    pd_check_notnull(SRCLOC,"t->interior_cross",t->interior_cross);
  } else {
    if (t->interior_cross != NULL) { 
      pd_error(SRCLOC,"t->ninterior_cross is 0, but the buffer t->interior_cross is allocated\n",pd);
    }
  }

  if (t->ninterior_edges != 0) { 
    pd_check_notnull(SRCLOC,"t->interior_edge",t->interior_edge);
  } else {
    if (t->interior_edge != NULL) { 
      pd_error(SRCLOC,"t->ninterior_edges is 0, but the buffer t->interior_edge is allocated\n",pd);
    }
  }

  if (t->ninterior_cross == PD_UNSET_IDX) { 

    pd_error(SRCLOC,"t->ninterior_cross is set to PD_UNSET_IDX (maybe not regenerated?)\n",pd);
    return false;

  } 

  if (t->ninterior_cross > pd->ncross) { 

    pd_error(SRCLOC,"t->ninterior_cross = %d > # of crossings in entire pd code (%d)\n",pd,
	     t->ninterior_cross,pd->ncross);
    return false;

  }

  if (t->ninterior_edges == PD_UNSET_IDX) { 

    pd_error(SRCLOC,"t->ninterior_edges is set to PD_UNSET_IDX (maybe not regenerated?)\n",pd);
    return false;

  } 

  if (t->ninterior_edges > pd->nedges) { 

    pd_error(SRCLOC,"t->ninterior_edges = %d > # of edges in entire pd code (%d)\n",pd,
	     t->ninterior_edges,pd->nedges);
    return false;

  }
  
  /* Now we run through the buffers and check them: */

  for(i=0;i<t->ninterior_cross;i++) { 

    pd_check_cr(SRCLOC,pd,t->interior_cross[i]);
    
  }

  for(i=0;i<t->ninterior_edges;i++) { 

    pd_check_edge(SRCLOC,pd,t->interior_edge[i]);

  }

  /* Now we check ordering. */

  for(i=1;i<t->ninterior_cross;i++) { 

    if (t->interior_cross[i-1] > t->interior_cross[i]) { 

      pd_error(SRCLOC,"list of interior crossings in %TANGLE "
	       "is not in order",pd,t);

    } 

  }

  for(i=1;i<t->ninterior_cross;i++) { 

    if (t->interior_cross[i-1] == t->interior_cross[i]) { 

      pd_error(SRCLOC,"list of interior crossings in %TANGLE "
	       "contains duplicate element at positions %d-%d",pd,t,i-1,i);

    } 

  }

  for(i=1;i<t->ninterior_edges;i++) { 

    if (t->interior_edge[i-1] > t->interior_edge[i]) { 

      pd_error(SRCLOC,"list of interior crossings in %TANGLE "
	       "is not in order",pd,t);

    }

  }

  for(i=1;i<t->ninterior_edges;i++) { 

    if (t->interior_edge[i-1] == t->interior_edge[i]) { 

      pd_error(SRCLOC,"list of interior crossings in %TANGLE "
	       "contains duplicate element at positions %d-%d",pd,t,i-1,i);

    }

  }
    
  /* We've checked everything that we can! */

  return true;
   
}

void tangle_contents_worker(pd_code_t *pd,pd_tangle_t *t,pd_idx_t edge,pd_idx_t depth);
/* The actual code is deferred to the end of the function below for 
   clarity. */

	       
void pd_regenerate_tangle(pd_code_t *pd, pd_tangle_t *t) 

/* The main purpose of this code is to fill in the "interior" data of
   a tangle from the boundary data. (If the boundary data isn't
   specified correctly, there's nothing we can do, because we won't
   know what part of the pd code is supposed to be the tangle in the
   first place). */


{
  pd_idx_t i,j;

  /* First, we check that the list of edges and faces are ok. */

  if (t->nedges == 0 || t->nedges > pd->nedges) { 

    pd_error(SRCLOC,"value for # of tangle edges (%d) doesn't "
	     "make sense for %d edge PD code",pd,t->nedges);
    exit(1);

  }

  for(i=0;i<t->nedges;i++) { 

    pd_check_edge(SRCLOC,pd,t->edge[i]);

  }

  for(i=0;i<t->nedges;i++) { 

    pd_check_face(SRCLOC,pd,t->face[i]);

  }

  /* We now check that the edges are, in fact, on the corresponding
     faces. */

  for(i=0;i<t->nedges;i++) { 

    bool found_edge = false;

    for(j=0;j<pd->face[t->face[i]].nedges;j++) { 

      if (pd->face[t->face[i]].edge[j] == t->edge[i]) { found_edge = true; }

    }

    if (!found_edge) { 

      pd_error(SRCLOC,"%EDGE and %FACE are paired on tangle, but this edge is not on this face in %PD",pd,t->edge[i],t->face[i]);
      exit(1);

    }

  }

  /* Ok, if we've survived this far, we can hope to reconstruct orientations and 
     interior crosing/edge data. We'll start with the innie/outie data. There are
     basically two cases here:

                           +                           
			   |                           
                           |                           
     +--------<---------------------<-----------------+
     |                     |                          |
     |                     |  f[i]                    ^
     V                     |                          |
     |           +---------+                          |
     |           |                                    |
     |           |                                    |
     +----*---edge e[i]--->*-----face boundary-->-----+
                 |                                     
		 |                                     
		 +                                     
		 
       (orientation of e[i] is positive along f[i]: OUT)

     and the other case, 
 
                           +                           
			   |                           
			   |                           
     +--------<---------------------<-----------------+
     |                     |                          |
     |                     |  f[i]                    ^
     V                     |                          |
     |           +---------+                          |
     |           |                                    |
     |           |                                    |
     +----*<--edge e[i]---*-----face boundary-->-----+
                 |                                     
		 |                                     
		 +                    
		 
       (orientation of e[i] is negative along f[i]: IN)
       
  */

  for(i=0;i<t->nedges;i++) { 

    pd_idx_t pos_face,pfp,neg_face,nfp;
    pd_face_and_pos(pd,t->edge[i],&pos_face,&pfp,&neg_face,&nfp);

    if (pos_face == t->face[i]) { 
      t->edge_bdy_or[i] = in;
    } else if (neg_face == t->face[i]) { 
      t->edge_bdy_or[i] = out;
    } else {
      pd_error(SRCLOC,"edge e[%d] = %EDGE in tangle is "
	       "not on face f[%d] = %FACE in pd\n",
	       pd,i,t->edge[i],i,t->face[i]);
      exit(1);
    }

  }

  /* We now have the inward and outward orientations set. We now go ahead 
     and construct the strands of the tangle. */

  assert(t->strand != NULL);  // We're assuming that the buffer of strands
                              // was allocated on creation of the tangle.

  pd_idx_t s;

  for(i=0,s=0;i<t->nedges;i++) { 

    if (t->edge_bdy_or[i] == in) { /* This is the start of a strand! */

      pd_idx_t pos;
      t->strand[s].nedges = 1;
      t->strand[s].start_edge = i;  // Start and end are positions on TANGLE
      pd_component_and_pos(pd,t->edge[i],&(t->strand[s].comp),&pos);

      pd_idx_t j;
      bool found_end;

      for(j=pd_next_edge(pd,t->edge[t->strand[s].start_edge]),found_end = false;
	  !found_end && j != t->edge[t->strand[s].start_edge]; // failsafe
	  j=pd_next_edge(pd,j),
	    t->strand[s].nedges++) {

	/* Now we have to search the list of edges of the tangle 
	   to try to find this particular edge j. Keep in mind that j may 
	   appear twice (in which case, we want the one with the OUT 
	   orientation). */

	pd_idx_t si,nfound=0,tpos[2];

	for(si=0;si<t->nedges;si++) { 
	  if (t->edge[si] == j) { 
	    if (nfound == 2) { 
	      pd_error(SRCLOC,"found edge %d more than 2 times in edges of %TANGLE",pd,j,t);
	    } else {
	      tpos[nfound] = si;
	      nfound++;
	    }
	  }
	}

	if (nfound!=0) { 
	  
	  /* First, we make sure that if nfound == 2, we have only one OUT */
	  if (nfound == 2) { 
	    
	    if (t->edge_bdy_or[tpos[0]] == t->edge_bdy_or[tpos[1]]) { 
	    
	      pd_error(SRCLOC,"found %EDGE twice on boundary of %TANGLE, "
		       "but with the same orientation %EDGE_BDY_OR both times",
		       pd,t->edge[tpos[0]],t,t->edge_bdy_or[tpos[0]]);
	      exit(1);

	    } 

	    if (t->edge_bdy_or[tpos[1]] == out) { /* Make sure the OUT instance is first */

	      tpos[0] = tpos[1];

	    }

	  } 

	  pd_idx_t bdypos = tpos[0];

	  if (!(t->edge_bdy_or[bdypos] == out)) { /* Just to be safe */
	    
	    pd_error(SRCLOC,"tangle strand which started (in) with %EDGE "
		     "has encountered tangle edge %d (%EDGE) with boundary "
		     "orientation not 'out' \n",pd,
		     t->strand[s].start_edge,bdypos,t->edge[bdypos]);
	    exit(1);

	  }

	  found_end = true;
	  t->strand[s].end_edge = bdypos;

	}

      }

      if (!found_end) { 

	pd_error(SRCLOC,
		 "component %COMP of pd appears "
		 "to enter tangle at least once (at tangle edge %EDGE), "
		 " but never leave it",
		 pd,t->strand[s].comp,t->strand[s].start_edge);

	exit(1); 

      }

      /* If we've survived to here, we completed the strand. */ 

      s++;

    }

  }

  /* 
     We've now built all the strands. It's going to be relatively easy
     to detect interior crossings and edges by a variant of the
     "contagion via crossings" algorithm that we used for the
     Reidemeister II code.

     The basic idea is that we'll start at a single (inward) strand, and 
     march along. At each (interior) crossing, we'll recurse, exploring 
     the forward direction along each edge leaving the crossing. 

     The recursion stops if we encounter a crossing we've already
     seen.  We will seed from each inward edge. The danger of the
     thing is that it's slow: we're going to have to maintain the list
     of interior crossings and edges in sorted order and do in-order
     inserts into the lists in order to keep the searching from
     consuming too much time.
  */

  t->ninterior_cross = 0;
  t->ninterior_edges = 0;
  
  if (t->interior_cross != NULL) { 

    free(t->interior_cross);
    t->interior_cross = NULL;

  }

  if (t->interior_edge != NULL) { 

    free(t->interior_edge);
    t->interior_edge = NULL;

  }

  for(i=0;i<t->nedges;i++) { 

    if (t->edge_bdy_or[i] == in) { 

      tangle_contents_worker(pd,t,t->edge[i],0);

    }

  }

  /* If we got here, we've set the interior stuff (hopefully correctly!). */

  assert(pd_tangle_ok(pd,t));

}

void tangle_contents_worker(pd_code_t *pd,pd_tangle_t *t,pd_idx_t edge,pd_idx_t depth)
/* 
   We are given an edge which is either inside the tangle or incident
   to the boundary of the tangle. 

   If the edge is leading out of the tangle, we simply return. 

   Otherwise, we proceed to the next crossing and try to add it to the 
   list of interior crossings. 

          If it's already listed, we return;

	  Otherwise, we add it and recurse along outward edges.

*/
{
  pd_idx_t tangle_bdy_pos;

  /* First, make sure the failsafe is working. */

  if (depth > pd->nedges) { 

    pd_error(SRCLOC,
	     "attempting to compute contents of tangle recursively\n"
	     "reached %d levels of recursion in a %d edge pd code\n"
	     "something is wrong (bug in tangle_contents_worker?)\n\n"
	     "pd is %PD\n",
	     pd,depth,pd->nedges);
    exit(1);

  }

  /* Now see if we're exiting the tangle here. */

  if (pd_buf_contains(t->edge,t->nedges,edge,&tangle_bdy_pos)) { 

    if (t->edge_bdy_or[tangle_bdy_pos] == out) { 

      return;

    }

  } else { /* We're _in_ the tangle, so add this edge to the interior edges. */

    t->interior_edge = pd_sortedbuf_insert(t->interior_edge,&(t->ninterior_edges),edge);
    
  }
    

  /* Ok, we're entering or in the tangle. We first set the crossing, 
     and see if we've already seen it. */

  pd_idx_t cr = pd->edge[edge].head;

  if (pd_sortedbuf_contains(t->interior_cross,t->ninterior_cross,cr)) { 

    return;

  }

  /* If we haven't seen it, we need to add the crossing to the list of interior crossings. */

  t->interior_cross = pd_sortedbuf_insert(t->interior_cross,&(t->ninterior_cross),cr);

  /* Now we'll need to figure out the outgoing edges and recurse on them. 
     The first recursion target is easy-- it's just the next edge. */

  tangle_contents_worker(pd,t,pd_next_edge(pd,edge),depth+1);

  /* The second recursion target is more challenging. First, we obtain
     the edge numbers of the "transverse" edges to our guy at cr. */

  pd_idx_t transverse_edge[2];

  transverse_edge[0] = pd->cross[cr].edge[(pd->edge[edge].headpos+1)%4];
  transverse_edge[1] = pd->cross[cr].edge[(pd->edge[edge].headpos+3)%4];
    
  if (pd->edge[transverse_edge[0]].tail == cr) { 

    tangle_contents_worker(pd,t,transverse_edge[0],depth+1);

  } else if (pd->edge[transverse_edge[1]].tail == cr) {

    tangle_contents_worker(pd,t,transverse_edge[1],depth+1);

  } else {

    pd_error(SRCLOC,
	     "something unexpected has happened. Edges %EDGE and %EDGE\n"
	     "are supposed to be head-to-tail and meet at crossing %CROSS\n"
	     "but neither seems to be leaving this crossing. Quitting.\n",
	     pd,transverse_edge[0],transverse_edge[1],cr);
    exit(1);

  }
  
}


pd_code_t *pd_flype(pd_code_t *pd,pd_idx_t e[4], pd_idx_t f[4])

  /*  This is a classical "flype" move, which we define by giving input
      and output edges. These must be connected in the (topological) 
      situation below.

                         f[0]

                           +---------+          
      +--e[0]--+    +------+         +---e[3]--+
               |    |      |         |          
               |    |      |         |          
    f[1]    +--------+     |  VVVVV  |    f[3]        
            |  |           |         |          
            |  |           |         |          
      +-e[1]+  +-----------+         +----e[2]--+
                           +---------+          
                      f[2]

		      |
		      v
         +---------+                            
         |         |                            
      +--+         +------+  +-----------------+
         |         |      |  |                  
         |  ^^^^^  |      |  |                  
         |         |      +------+              
         |         |         |   |              
      +--+         +---------+   +-------------+
         +---------+                            

  */

{
  pd_idx_t i;
  assert(pd != NULL);
  for(i=0;i<4;i++) { pd_check_edge(SRCLOC,pd,e[i]); pd_check_face(SRCLOC,pd,f[i]); }

  return NULL;
}

  /* A flype operates on a tangle in the form

                        f[0]

                           +---------+          
      +-e[0]---+    +------+         +----e[3]--+
               |    |      |         |          
               |    |      |         |          
   f[1]     +--------+     |  VVVVV  |        f[3]         
            |  | cr        |         |          
            |  |           |         |          
      +-e[1]+  +-----------+         +----e[2]--+
                           +---------+          

                        f[2]

     We start by collecting all of the crossings and edges
     involved in the tangle.

  */
  

  /* Ok, that's all the self-checking that we can do... 
     We have now identified the crossing, and are ready to do the actual sewing. 
    

                          f[0]

                           +---------+          
      +--e[0]--+    +-ae[1]+         +---e[3]--+
               |    |      |         |          
               |    |      |         |          
    f[1]    +-------+      |  VVVVV  |    f[3]        
            |  | cr        |         |          
            |  |           |         |          
      +-e[1]+  +---ae[0]---+         +----e[2]--+
                           +---------+          
                      f[2]

		      |
		      v
         +---------+                            
         |         |                            
    e[0]-+         +--ae[1]  +--------------e[3]
         |         |      |  |                  
         |  ^^^^^  |      |  | cr                 
         |         |      +------+              
         |         |         |   |              
    e[1]-+         +--ae[0]--+   +----------e[2]
         +---------+                            
  
    We're going to identify "adjacent" edges ae[0] and ae[1] 
    across from e[0] and e[1] at cr. (This is going to depend on 
    the boundary orientation of e[0] and e[1], obviously.) 

    Then we're going to sew e[0] and ae[0] together (eliminating ae[0])
    and likewise for e[1] and ae[1] (likewise eliminating ae[1]).
    Afterwards, we're going to splice in these edges on the other 
    side of the tangle. 

    This is going to require us to walk the tangle strands again
    in order to shift edge indices, and hence to HAVE the strands.
    
  */

/*************************************************************************
   Tangle Slide Code. This is broken into several subprocedures so that 
   we can write unit tests for the individual pieces. 
  **********************************************************************/


bool pdint_check_tslide_data_ok_and_find_te(pd_code_t *pd,pd_tangle_t *t,
					    pd_idx_t n,
					    pd_idx_t *overstrand_edges, 
					    pd_idx_t *border_faces,
					    pd_idx_t **tangle_slide_edges)

/* We make sure that the data is appropriate for a 
   tangle slide. 

		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
             | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |
       |       tse[n-1]    tse[0]     |
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]---<----e[1]-----e[0]-+	      
       	       	 |   | 	     |		   

 
   We need to make several checks:

   0) The indices in overstrand_edges and border_faces are
      in-range edge and face indices for the pd code pd.

   1) The border faces are actually border faces of the 
      tangle, either in ccw or cw order.

   2) The overstrand edges are incident to the appropriate 
      faces, in order, and not interior to the tangle itself.

   3) The crossings on the tangle edges have corresponding 
      signs.

   The procedure returns only if all this is true; otherwise we quit
   with a call to pd_error. The output buffer tangle_slide_edges (of 
   size n-1) is allocated in this function and must be disposed of 
   externally.

*/
{
  pd_idx_t i;

  /* 0. Check indices of everything. */

  for(i=0;i<n;i++) { 

    pd_check_edge(SRCLOC,pd,overstrand_edges[i]);
    pd_check_face(SRCLOC,pd,border_faces[i]);

  }  

  /* 1. We start by looking for the sequence of faces occuring 
     along the boundary of the tangle. Notice that the same face
     may occur multiple times on the boundary of a tangle:
		  
	       +--+    f   +-----------+
	       |  |        |  	       |
       	     +-------------------+     |	  
	     |	     	      	 |--+  |	  
             | 	    Tangle     	 |  |  |	       	  
         f   | 	       	       	 |--+  |     	    
             | 	     	      	 |     |   	     
             +-------------------+     |     
                 |   |       | 	       |     
                 |   |       |         |      
       	         +---+ 	f    +---------+   
                                  
     and we can even have a sequence of faces occur more than
     once   

		     |	   |
  +----------------+ |	   |
  |   +----------+ | |	   |
  |   |   +----+ | | |     |
  |   |	  |    | | | |     |   	       	 
  |   |	  |  +-------------------+      	  
  |   |	  |  |	      	      	 |      	  
  |   |   |  | 	    Tangle     	 |      	       	  
  |f2 |f1 |f0| 	       	       	 |           	    
  |   |   |  | 	      	      	 |         	     
  |   |   |  +-------------------+           
  |   |   |    | | | |       | 	             
  |   |   +----+ | | |       |                 	 
  |   +----------+ | |       |                     
  +----------------+ |	     |   

      
  */

  bool found_sequence = false;
  pd_idx_t j;

  /* First we look for the sequence in the forward orientation. */
	      	       			      
  for(i=0;i<t->nedges;i++) {  
			      
    if (t->face[i] == border_faces[0]) { 

      /* This is a possible start-- we check for the 
	 border_faces sequence in positive orientation. */

      bool match_continues = true;

      for (j=1;j<n && match_continues;j++) { 

	match_continues = (t->face[(i+j)%t->nedges] == border_faces[j]);

      }

      if (match_continues) { found_sequence = true; }

    }

  }

  /* We now look for the sequence in the backward orientation. */

  for(i=0;i<t->nedges;i++) {  
			      
    if (t->face[i] == border_faces[n-1]) { 

      /* This is a possible start-- we check for the 
	 border_faces sequence in positive orientation. */

      bool match_continues = true;

      for (j=1;j<n && match_continues;j++) { 

	match_continues = (t->face[(i+j)%t->nedges] == border_faces[n-1-j]);

      }

      if (match_continues) { found_sequence = true; }

    }

  }

  if (!found_sequence) { 

    if (n >= 4) { 

      pd_error(SRCLOC,
	       "couldn't find %d border_faces sequence %d %d %d ... %d on the border\n"
	       "of %TANGLE.\n",pd,n,border_faces[0],border_faces[1],border_faces[2],border_faces[n-1],t);

    } else { 

      pd_error(SRCLOC,
	       "couldn't find %d border_faces sequence %d ... %d on the border\n"
	       "of %TANGLE.\n",pd,n,border_faces[0],border_faces[n-1],t);

    }

    return false;

  }

  /* 2. We now check to see whether the overstrand edges are actually on 
     the faces in question.

		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
             | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]---<----e[1]-----e[0]-+	      
       	       	 |   | 	     |		   

     Notice that however this goes, either the overstrand edges should
     all occur positively or all negatively along the faces in question. 

  */

  pd_idx_t posface, pfp, negface, nfp;

  pd_face_and_pos(pd,overstrand_edges[0],
		  &posface, &pfp, &negface, &nfp);

  if (posface == border_faces[0]) { /* We're go in the positive direction. */

    pd_idx_t j;

    for(j=1;j<n;j++) { 

      pd_face_and_pos(pd,overstrand_edges[j],
		      &posface, &pfp, &negface, &nfp);

      if (posface != border_faces[j]) { 

	pd_error(SRCLOC,
		 "overstrand_edges[0] = %EDGE occurs on border_face[0] = %FACE positively, but\n"
		 "overstrand_edges[%d] = %EDGE does not occur on border_face[%d] = %FACE positively\n"
		 "meaning that overstrand and border_face aren't consistent.\n",
		 pd,overstrand_edges[0],border_faces[0],j,overstrand_edges[j],border_faces[j]);

	return false;

      }

    }

  } else if (negface == border_faces[0]) { /* We're go in the negative direction. */

    pd_idx_t j;

    for(j=1;j<n;j++) { 

      pd_face_and_pos(pd,overstrand_edges[j],
		      &posface, &pfp, &negface, &nfp);

      if (negface != border_faces[j]) { 

	pd_error(SRCLOC,
		 "overstrand_edges[0] = %EDGE occurs on border_face[0] = %FACE negatively, but\n"
		 "overstrand_edges[%d] = %EDGE does not occur on border_face[%d] = %FACE negatively\n"
		 "meaning that overstrand and border_face aren't consistent.\n",
		 pd,overstrand_edges[0],border_faces[0],j,overstrand_edges[j],border_faces[j]);

	return false;

      }

    }

  } else { 

    pd_error(SRCLOC,
	     "overstrand_edges[0] = %EDGE does not occur on border_face[0] = %FACE\n"
	     "meaning that overstrand and border_face aren't consistent.\n",
	     pd,overstrand_edges[0],border_faces[0]);
    
  }

  /* We now check that the overstrand edges are not interior to 
     the tangle itself and that overstrand edges 1..n-1 are not 
     border strands either. */

  for(i=0;i<n;i++) { 

    pd_idx_t j;

    for(j=0;j<t->ninterior_edges;j++) { 

      if (t->interior_edge[j] == overstrand_edges[i]) { 

	pd_error(SRCLOC,"overstrand_edges[%d] == %EDGE is interior to %TANGLE\n",
		 pd,i,overstrand_edges[i],t);
	return false;

      }

    }

  }

  for(i=1;i<n-1;i++) { 

    pd_idx_t j;

    for(j=0;j<t->nedges;j++) { 

      if (t->edge[j] == overstrand_edges[i]) { 

	pd_error(SRCLOC,"overstrand_edges[%d] == %EDGE is an incident edge of %TANGLE\n",
		 pd,i,overstrand_edges[i],t);
	return false;

      }

    }

  }

  /* Now make sure that the overstrand edges are actually in orientation order. */

  for(i=0;i<n-1;i++) { 

    if (pd_next_edge(pd,overstrand_edges[i]) != overstrand_edges[i+1]) { 

      pd_error(SRCLOC,"overstrand edges %d and %d (%EDGE and %EDGE) are not\n"
	       "in orientation order along the their component in %PD",pd,
	       i,i+1,overstrand_edges[i],overstrand_edges[i+1]);
      return false;

    }

  }
  
  /* Finally, we're ready to check for tangle edges at each crossing. This is 
     complicated by the fact that BOTH of the edges at the crossing may be tangle
     edges, as in:

       	+----------+
       	|	   |	     +----e[0]
       	|    +--------+	     |	      
       	| e[n-1]      |------|----    
	| +--|   T    |      | 	      
	| |  |	      |------|----    
	| |  +--------+	     | 	      
	| |  	  | 	     |
     	| +------------------+
	+---------+

     However, we know that only ONE of these two tangle edges should be incident
     to the corresponding border face:

     		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
             | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]---<----e[1]-----e[0]-+	      
       	       	 |   | 	     |		   

     So our strategy is to figure out which of the crossing edges is incident
     to the correct face, then look for it in the tangle's edges.

  */

  pd_idx_t *te;
  te = calloc(n-1,sizeof(pd_idx_t));
  assert(te != NULL);

  for(i=0;i<n-1;i++) {

    pd_idx_t candidate_edges[2];
    pd_idx_t oe = overstrand_edges[i];

    candidate_edges[0] = pd->cross[pd->edge[oe].head].edge[(pd->edge[oe].headpos+1)%4];
    candidate_edges[1] = pd->cross[pd->edge[oe].head].edge[(pd->edge[oe].headpos+3)%4];

    if (pd_edge_on_face(pd,candidate_edges[0],border_faces[i]) && 
	pd_edge_on_face(pd,candidate_edges[0],border_faces[i+1])) { 

      te[i] = candidate_edges[0];

    } else if (pd_edge_on_face(pd,candidate_edges[1],border_faces[i]) && 
	       pd_edge_on_face(pd,candidate_edges[1],border_faces[i+1])) { 

      te[i] = candidate_edges[1];

    } else {

      pd_error(SRCLOC,
	       "at overstrand_edges[%d] == %EDGE, neither of the\n"
	       "candidate edges %EDGE and %EDGE are incident to\n"
	       "both border_faces[%d] == %FACE and border_faces[%d] == %FACE\n",
	       pd,i,overstrand_edges[i],candidate_edges[0],candidate_edges[1],
	       i,border_faces[i],i+1,border_faces[i+1]);

      return false;

    }

    /* We now check that te[i] is indeed a tangle edge. */

    pd_idx_t j;
    bool found = false;

    for(j=0;j<t->nedges && !found;j++) {

      if (t->edge[j] == te[i]) { found = true; }

    }

    if (!found) { 

      pd_error(SRCLOC,
	       "at %CROSS joining overstrand_edges[%d] and [%d] (%EDGE and %EDGE),\n"
	       "%EDGE is incident to both border faces %FACE and %FACE, but not \n"
	       "part of %TANGLE\n",
	       pd,pd->edge[overstrand_edges[i]].head,i,i+1,overstrand_edges[i],overstrand_edges[i+1],
	       te[i],border_faces[i],border_faces[i+1],t);
      return false;

    }

  }

  /* Last, we check that the run of overstrand edges is actually over (or under) */
  /* We start by seeing what the first crossing does, then loop to check the */
  /* remaining crossings are the same. */

  /* Since we know the overstrand_edges are in orientation order, if the overstrand
     is going OVER, it is the incoming_edgenum in pd_overstrand:

  void pd_overstrand(pd_code_t *pd,pd_idx_t cr,pd_idx_t *incoming_edgenum, pd_idx_t *outgoing_edgenum);

  */
  
  pd_idx_t incoming_e, outgoing_e;

  pd_overstrand(pd,pd->edge[overstrand_edges[0]].head,&incoming_e,&outgoing_e);

  if (incoming_e == overstrand_edges[0]) { /* We're supposed to be always OVER. */

    pd_idx_t j;

    for(j=1;j<n-1;j++) { 

      pd_overstrand(pd,pd->edge[overstrand_edges[j]].head,&incoming_e,&outgoing_e);

      if (incoming_e != overstrand_edges[j]) { 

	pd_error(SRCLOC,
		 "sequence of overstrand_edges in tangle starts with OVERcrossing at\n"
		 "head of overstrand_edges[0] == %EDGE (%CROSS), but overstrand_edges[%d]\n"
		 "== %EDGE goes UNDER at its head (%CROSS)\n",
		 pd,overstrand_edges[0],pd->edge[overstrand_edges[0]].head,j,
		 overstrand_edges[j],pd->edge[overstrand_edges[j]].head);
	return false;

      }

    }

  } else { /* Presumably, we're going UNDER, but we check that the orientation isn't
	      screwed up somehow just to be sure. */

    pd_understrand(pd,pd->edge[overstrand_edges[0]].head,&incoming_e,&outgoing_e);

    if (incoming_e == overstrand_edges[0]) { /* We're going UNDER. */

      pd_idx_t j;

      for(j=1;j<n-1;j++) { 
	
	pd_understrand(pd,pd->edge[overstrand_edges[j]].head,&incoming_e,&outgoing_e);

	if (incoming_e != overstrand_edges[j]) { 

	  pd_error(SRCLOC,
		 "sequence of overstrand_edges in tangle starts with UNDERcrossing at\n"
		 "head of overstrand_edges[0] == %EDGE (%CROSS), but overstrand_edge[%d]\n"
		 "== %EDGE goes OVER at its head (%CROSS)\n",
		 pd,overstrand_edges[0],pd->edge[overstrand_edges[0]].head,j,
		 overstrand_edges[j],pd->edge[overstrand_edges[j]].head);
	  return false;

	}

      }

    } else {

      pd_idx_t ieo,oeo,ieu,oeu;

      pd_overstrand(pd,pd->edge[overstrand_edges[0]].head,&ieo,&oeo);
      pd_understrand(pd,pd->edge[overstrand_edges[0]].head,&ieu,&oeu);
    
      pd_error(SRCLOC,
	       "orientation or other internal error: edge %EDGE is neither the\n"
	       "incoming overstrand %d or incoming understrand %d at its head\n"
	       "(%CROSS).\n",
	       pd, overstrand_edges[0],ieo,ieu,pd->edge[overstrand_edges[0]].head);
      
      return false;

    }

  }

  /* We've now checked everything. 
     We assign the tangle_slide_edges buffer address and quit. */

  *tangle_slide_edges = te;
  return true;

}

void pd_tangle_slide(pd_code_t *pd,pd_tangle_t *t,
		     pd_idx_t n,
		     pd_idx_t *overstrand_edges, 
		     pd_idx_t *border_faces,
		     pd_idx_t **npieces,
		     pd_code_t **pd_pieces)

 /* Given a list of edges overstrand_edges (e[0]...e[n-1], below) and
    corresponding faces bordering the tangle (f[0]...f[n-1], below),
    slide the strand over the tangle to cross the remaining edges
    of the tangle, as below. The edges e[0]..e[n-1] are supposed to 
    occur in orientation order along their component.
    

		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
             | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]---<----e[1]-----e[0]-+	      
       	       	 |   | 	     |		   


    becomes

		  |    	   |
       +------------------------------+
       |	  |  	   |	      |
       |     +-------------------+    |		  
       |     |		    	 |    |		  
       |     | 	    Tangle     	 |    | 	       	  
    ---+     | 	       	       	 |    +---   	    
             | 	    	       	 |      	     
             +-------------------+     	     
                 |   |       | 	       	     
                 |   |       |         	      
       	       	 |   | 	     |		   
    			     
    We handle correctly the case where the initial and/or final
    edges of the strand are tangle edges themselves. We also note
    that while we call the strand the "overstrand", we also handle
    the case where the strand goes UNDER all of the tangle strands.

    This operation can disconnect the diagram, potentially into many
    pieces. We return the number of connected components of the
    diagram in "npieces" and the components themselves in
    "pd_pieces". Both are allocated internally and are the caller's
    responsibility to dispose of.
     	       	       	    
   */
{

  /* We start by checking the input, and coming up with a list of 
     tangle_slide_edges which are crossed by the overstrand. This 
     list is in orientation order along the overstrand, as below:
  
		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
             | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]---<----e[1]-----e[0]-+	      
       	       	 |   | 	     |		   
               tse[n-1]    tse[0] 



  */

  pd_idx_t *tangle_slide_edges;

  if (!pdint_check_tslide_data_ok_and_find_te(pd,t,n,
					      overstrand_edges, 
					      border_faces,
					      &tangle_slide_edges)) { 

    pd_error(SRCLOC,"tangle slide input data is invalid");
    exit(1);

  }

  /* Step 1. There are number of very ugly special cases in the forms
					         
     	   +--------------+		       +--------------+
           |		  |	       	       |	      |
	   |   Tangle T   |    	       	       |   Tangle T   |
	   |	          |----+       	       |  	      |
       	   |         	  |    |       ---+    | 	      |
  ----+    +--------------+    |       	  |    +--------------+
      |	      |	   |   |       |      	  |	    |  	   |   
      |	      |	   |   |       |      	  +----------------------+
      +-------|----|---|-------+       	       	    |  	   |   	 | 
	      |	   |   |   e[0]	      		    |  	   |	 |e[0]
			   				   +-----+     

     The problem is that these cases lack an "anchor" crossing on the
     other end of e[0] which is neither part of the tangle nor part of
     the overstrand. (Of course, the situation is just as bad if
     e[n-1] lacks an anchor vertex.) 

     For the moment, we're going to assume that there's an anchor vertex.
     Later, we'll do an R1 loop insert to create an anchor vertex, but we 
     haven't coded that yet. This test simply checks for presence of an 
     anchor vertex.

     We'll call e[0] and e[n-1]'s anchor crossing the "start_anchor" and 
     "end_anchor".

  */

  pd_idx_t start_anchor, end_anchor;

  start_anchor = pd->edge[overstrand_edges[0]].tail;
  end_anchor   = pd->edge[overstrand_edges[n-1]].head;

  pd_idx_t i;
  bool start_anchor_ok = true,end_anchor_ok = true;

  /* First, we check the crossings involved in the overstrand. */

  for(i=1;i<n-1;i++) { 

    if (pd->edge[overstrand_edges[i]].tail == start_anchor ||
	pd->edge[overstrand_edges[i]].head == start_anchor ) {

      start_anchor_ok = false;

    } 

    if (pd->edge[overstrand_edges[i]].tail == end_anchor ||
	pd->edge[overstrand_edges[i]].head == end_anchor) { 

      end_anchor_ok = false;

    }

  }

  /* Next, we check the crossings internal to the tangle T. */

  for(i=0;i<t->ninterior_cross;i++) { 

    if (t->interior_cross[i] == start_anchor) { start_anchor_ok = false; }
    if (t->interior_cross[i] == end_anchor) { end_anchor_ok = false; }

  }

  /* Now we act on the information we've got. */

  if (!start_anchor_ok) {

    pd_error(SRCLOC,
	     "this version of pd_tangle_slide requires that the overstrand start\n"
	     "at an anchor crossing which is not part of the remainder of the overstrand\n"
	     "or the interior of the tangle.\n"
	     "overstrand_edges[0] = %EDGE\n"
	     "starts at %CROSS \n"
	     "which is part of overstrand edges or interior to %TANGLE\n",
	     pd,overstrand_edges[0],overstrand_edges[0].tail,t);
    exit(1);

  }

  if (!end_anchor_ok) {

    pd_error(SRCLOC,
	     "this version of pd_tangle_slide requires that the overstrand end\n"
	     "at an anchor crossing which is not part of the remainder of the overstrand\n"
	     "or the interior of the tangle.\n"
	     "overstrand_edges[n-1] = %EDGE\n"
	     "ends at %CROSS \n"
	     "which is part of overstrand edges or interior to %TANGLE\n",
	     pd,overstrand_edges[n-1],overstrand_edges[n-1].tail,t);
    exit(1);

  }
			   
  /* We are now going to start working on the edges. The idea is that
     we're going to delete n-1 crossings on the "overstrand" side of
     the tangle and add (at most) (t->nedges - (n-1)) crossings on the other
     side of the tangle.   

     Note that if the first and last overstrand edges are themselves
     tangle edges, the overstrand doesn't cross them after the tangle slide,
     so the number of crossings added is in the interval

     [t->nedges - n - 1, t->nedges - n + 1]

     This means that the net change in number of crossings is in the interval

     [t->nedges - 2n,t->nedges - 2n + 2]

     We don't actually know that the interval is sorted this way,
     because the net change in crossings may be negative. The output
     pd_code is going to have 2*crossings edges in the end, so a first
     idea would be to simply prepare for this many crossings in the
     new pd. 

     However, this causes some problems because we then have to be
     very careful about order of operations (all deletions from the
     edge list have to happen first) and keeping track of where in the
     edge array to add edges (if we are to add them in order).

     The strategy here is going to be to trade space for time: we're going 
     to expand the edge array to make space for 

     2(t->nedges - n + 1)

     NEW edges. We're then going to add the edges with new edge
     numbers, and delete edges in place as we go, replacing deleted
     edges with pd_unset_edge(). Once all edge inserts and deletions
     are complete, we'll run the components and renumber everything.

     We'll do the same thing with crossings: expanding the buffer of 
     crossings to make room for 

     t->nedges - n + 1

     new crossings (even though the net change may not be this large,
     and may even be negative!).

  */

  /* We start by making a working copy of the pd code with room for the 
     new edges. We do this by copying and then realloc'ing the edge buffer
     because we may later add stuff to the pd_code_t which will be copied
     by copy, and don't want to duplicate the functionality of that code 
     here (as it probably won't get updated!). */

  pd_code_t *pd_working = pd_copy(pd); 

  pd_working->MAXEDGES = pd->nedges + 2*(t->nedges - n + 1);
  pd_working->edge = realloc(pd_working->edge,pd_working->MAXEDGES*sizeof(pd_edge_t));
  assert(pd_working->edge != NULL); 

  pd_working->MAXVERTS = pd->ncross + (t->nedges - n + 1);
  pd_working->cross = realloc(pd_working->cross,pd_working->MAXVERTS*sizeof(pd_crossing_t));
  assert(pd_working->cross != NULL);


  /* Step 2. Delete the internal edges of the overstrand. 

		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
             | 	    Tangle     	 |      	       	  
end_anchor   |                   |    +---start_anchor   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       |      tse[n-1]    tse[0]      |
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]-+   +       +---e[0]-+	      
       	       	 |   | 	     |		   
               
  */

  pd_idx_t i;

  for(i=1;i<n-1;i++) { 

    pd_delete_edge(pd_working,overstrand_edge[i]);

  }

  /* Step 3. We're now in the state: 

		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
 end_anchor  | 	    Tangle     	 |      	       	  
       |     | 	       	       	 |    +---start_anchor   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       |      tse[n-1]    tse[0]      |
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]-+   +       +---e[0]-+	      
       	       	 |   | 	     |		   
	
     We want to delete the crossings at the outward end of the 
     tangle_slide_edges. We've marked them-- each crossing we 
     want to get rid of is now adjacent to an unset edge. 
	
     We're going to delete them in chains-- if we have something
     like
	
		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
 end_anchor  | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor   	    
       |     | 	 +---+	     |	 |    | 	     
       |     +---------------|---+    |	     
       |      tse[n-1]    tse[0]      |
       | f[n-1]  |   |  f[1] v 	 f[0] |	     
       +--e[n-1]-+   +       +---e[0]-+	      
       	       	 |   | 	     |		  
                     +-------+
	   		      
     we'll keep going forward, deleting crossings as we go, until we 
     reach a crossing that's NOT marked for deletion, leaving us with 
     something like
	   
	   	  |  	   |
	   	  |  	   |
       	     +-------------------+		  
       	     | 	       	       	 |     	       	       	
end_anchor   | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor 
       |     | 	      	     |	 |    |     	     
       |     +---------------|---+    |	     
       |                  +--+ 	      |	    
       |                  |    	      |	     
       +--e[n-1]-         |   ---e[0]-+	      
       	       	       	  |   		  
                 +--------+   
	   	 |	      
	   
     The crossing deletion strategy will always be to go in
     orientation order, so we either start on tangle_slide_edge[i] (if
     that's oriented OUT), or on the predecessor edge to that (if
     tse[i] is oriented IN). This has some dangers:

     1) We could delete a component entirely. 


	   	  |  	   |
	   	  |  	   |
       	     +-------------------+		  
       	     | 	       	       	 |     	       	       	
end_anchor   | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor 
       |     | 	  +----------+	 |    |
       |     +-------------------+    |
       |       |  |          | |      |
       |       |  |          | |      |
       +-e[n-1]+  +          + +-e[0]-+
       	       |  |    	     | |
               |  +----------+ |
	   	  	      
                       |
                       v

	   	  |  	   |
	   	  |  	   |
       	     +-------------------+		  
       	     | 	       	       	 |     	       	       	
end_anchor   | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor 
       |     | 	              	 |    |
       |     +-------------------+    |
       |       |               |      |
       |       |               |      |
       +-e[n-1]+  +          + +-e[0]-+
       	       |      	       |
               |               | 

     2) We could separate one or more components from the rest of the link without 
     deleting them entirely.

	   	  |  	   |
	   	  |  	   |
       	     +-------------------+		  
       	     | 	       	       	 |     	       	       	
end_anchor   | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor 
       |     | 	  +----------+	 |    |
       |     +-------------------+    |
       |       |  |          | |      |
       |       |  |          | |      |
       +-e[n-1]+  +          + +-e[0]-+
       	       |  | +--+     | |
               |  | |  |     | |
		  | +--------+
		  -----+     

                       |
                       v
			      
	   	  |  	   |  
	   	  |  	   |  
       	     +-------------------+		  
       	     | 	       	       	 |     	       	       	
end_anchor   | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor 
       |     | 	              	 |    |
       |     +-------------------+    |
       |       |               |      |
       |       |               |      |
       +-e[n-1]+               +-e[0]-+
       	       |  +----------+ |
      	       |  | +--+     | |
               |  | |  |     | |
	      	  | +--------+
	      	  -----+      

     We're not going to code around these specifically, but we need to 
     stay aware that they can happen.
  */	   		      
	   		      
  pd_pos_t j;		      
	   		      
  for(i=0;i<n-1;i++) { 
	   
    /* First, we check if this edge has already been deleted. */
	   
    if (pd_working->edge[tangle_slide_edges[i]] != pd_unset_edge()) { 
	   
      pd_idx_t edgechain_start;
    	   
      if (pd_tangle_bdy_or(t,tangle_slide_edges[i]) == out) { 
	
	edgechain_start = tangle_slide_edges[i];

      } else {

	edgechain_start = pd_previous_edge(pd_working,tangle_slide_edges[i]);

      }

      /* Now follow the edgechain, deleting as we go. */

      pd_idx_t cur_edge_idx;
      pd_edge_t cur_edge; 
      for(cur_edge_idx = edgechain_start, cur_edge = pd_working->edge[cur_edge_idx];
	  /* Until we reach an edge whose head is not going to be deleted */
	  pd_working->cross[cur_edge.head].edge[(cur_edge.headpos+1)%4] 
	    != PD_WORKING_UNSET_IDX && 
	    pd_working->cross[cur_edge.head].edge[(cur_edge.headpos+3)%4] 
	    != PD_WORKING_UNSET_IDX;
	  /* Advance to the next edge in the chain */
	  cur_edge_idx = pd_next_edge(pd_working,cur_edge_idx), 
	    cur_edge = pd_working->edge[cur_edge]) { 

	pd_edge_delete(pd_working,cur_edge_idx);
      
      }

      /* We've now reached a cur_edge which has a head which is not to be deleted 
	 as a crossing. We're going to modify the start of the chain to end here. */

      pd_working->edge[edgechain_start].head = pd_working->edge[cur_edge_idx].head;
      pd_working->edge[edgechain_start].headpos = pd_working->edge[cur_edge_idx].headpos;

      pd_delete_edge(pd_working,cur_edge_idx); 
      
      /* Now this call to pd_delete_edge set the crossing incident to
	 cur_edge_idx at the head to have a PD_UNSET_IDX in the
	 position where cur_edge_idx came in. We want to set this
	 crossing to have the start of the edgechain coming in
	 anyway. */

      pd_working->cross[pd_working->edge[edgechain_start].head].
	edge[pd_working->edge[edgechain_start].headpos] = edgechain_start;

      /* We've killed the entire edgechain. It's time to look for the
         start of the next edgechain. */

    }

  }

   /* Step 4. We're now (hopefully) in the state: 

		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
 end_anchor  | 	    Tangle     	 |      	       	  
       |     | 	       	       	 |    +---start_anchor   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       |         |   |       |        |
       |         |   |       | 	      |	     
       +--e[n-1]+|  +|       |+--e[0]-+	      
       	       	 |   | 	     |		   
	
     where edges have been deleted and crossings set to 
     the state where they have all incident edges unset.
     At this point, we should have exactly two crossings
     which have three edges unset (the head of e[0] and
     the tail of e[n-1] and various crossings with four 
     edges unset (everything in the original overstrand).
     
     We're now going to take a moment to check that this
     state holds. */

  pd_idx_t UNSET_count;

  /* We start by checking the head of e[0]. */

  for(UNSET_count=0,i=0;i<4;i++) { 

    if (pd_working->cross[pd_working->edge[overstrand_edges[0]].head].edge[i] 
	== PD_UNSET_IDX) { UNSET_count++; }

  }

  if (UNSET_count != 3 && n != 2) { 

    pd_error(SRCLOC,
	     "after deleting edges on a %d-edge overstrand\n"
	     "head of overstrand_edges[0] = %EDGE (%CROSS) does not \n"
	     "have exactly three incident UNSET edges\n",pd_working,n,overstrand_edges[0],
	     pd_working->edge[overstrand_edges[0]].head);
    exit(1);

  }
  
  if (UNSET_count != 2 && n == 2) { 

    pd_error(SRCLOC,
	     "after deleting edges on a %d-edge overstrand\n"
	     "head of overstrand_edges[0] = %EDGE (%CROSS) does not \n"
	     "have exactly two incident UNSET edges\n",pd_working,n,overstrand_edges[0],
	     pd_working->edge[overstrand_edges[0]].head);
    exit(1);

  }

  /* Now we check the tail of e[n-1] */

  for(UNSET_count=0,i=0;i<4;i++) { 

    if (pd_working->cross[pd_working->edge[overstrand_edges[n-1]].tail].edge[i] 
	== PD_UNSET_IDX) { UNSET_count++; }

  }

  if (UNSET_count != 3 && n != 2) { 

    pd_error(SRCLOC,
	     "after deleting edges on a %d-edge overstrand\n"
	     "tail of overstrand_edges[n-1] = %EDGE (%CROSS) does not \n"
	     "have exactly three incident UNSET edges\n",pd_working,n,overstrand_edges[n-1],
	     pd_working->edge[overstrand_edges[n-1]].tail);
    exit(1);

  }
  
  if (UNSET_count != 2 && n == 2) { 

    pd_error(SRCLOC,
	     "after deleting edges on a %d-edge overstrand\n"
	     "head of overstrand_edges[0] = %EDGE (%CROSS) does not \n"
	     "have exactly two incident UNSET edges\n",pd_working,n,overstrand_edges[n-1],
	     pd_working->edge[overstrand_edges[n-1]].tail);
    exit(1);

  }

  /* Now we look at the middle of the overstrand. This is a little
     complicated, since we just deleted all those edges! However, we
     can look back at pd itself (still unchanged) and recover the
     crossing numbers. */

  for(i=1;i<n-1;i++) { 

    pd_idx_t cr = pd->edge[overstrand_edges[i]].head;

    for(UNSET_count=0,i=0;i<4;i++) { 

      if (pd_working->cross[cr].edge[i] 
	  == PD_UNSET_IDX) { UNSET_count++; }

    }
    
    if (UNSET_count != 4) { 

      pd_error(SRCLOC,
	       "after deleting edges on a %d-edge overstrand\n"
	       "head of overstrand_edge[%d] = %EDGE (%CROSS) does not \n"
	       "have all 4 incident edges UNSET\n",pd_working,n,overstrand_edges[i],
	       pd_working->edge[overstrand_edges[i]].head);
      exit(1);

    }
    
  }

  /* Step 5. Having checked our deletions, we're now going to insert new crossings
     and edges on the other side of the tangle. This involves several steps. First, 
     we're going to identify a list of "complementary edges". These are going be 
     in the order that they will be encountered in the new overstrand passage, which 
     could be reversed from their order in the t->edge buffer. Going back to the 
     original diagram, we HAD something like this:

       	       ce[k-1]    ce[0]
		  |  	   |
		  |  	   |
       	     +-------------------+		  
	     |		    	 |		  
 end_anchor  | 	    Tangle     	 |      	       	  
    ---+     | 	       	       	 |    +---start_anchor   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       | f[n-1]  |   |  f[1] | 	 f[0] |	     
       +--e[n-1]---<----e[1]-----e[0]-+	      
       	       	 |   | 	     |		   
               tse[n-1]    tse[0] 
  
      Once we have the complementary edges, we'll split them by introducing new 
      crossings just outside the tangle, and string these crossings together with 
      new edges (if needed) joining e[0] to e[n-1]. 

      The order of the complementary edges is important, but it can
      only be determined using the order of the border_faces, since
      there are guaranteed to be two border faces (but there may be
      only one tangle_slide_edge).

  */

  pd_idx_t step;
  bool found_f0 = false;
  pd_idx_t f0_tanglepos;

  for(i=0;i<t->nedges;i++) { 

    if (t->face[i] == border_faces[0]) { 

      found_f0 = true;
      f0_tanglepos = i;

    }

  }

  if (!found_f0) { 

    pd_error(SRCLOC,
	     "couldn't find border_faces[0] = %d in the "
	     "list of faces in %TANGLE",
	     pd_working,border_faces[0],t);
    
    exit(0);

  }

  if (t->face[(f0_tanglepos+1)%t->nedges] == border_faces[1]) { 

    step = -1; /* Go backwards. */

  } else { 

    step = +1; /* Go forwards. */

  }

  pd_idx_t *complementary_edges = calloc(t->nedges /* too big */,sizeof(pd_idx_t));
  assert(complementary_edges != NULL);
  pd_boundary_or_t *complementary_or = calloc(t->nedges, sizeof(pd_boundary_or_t));
  assert(complementary_or != NULL);
  pd_idx_t ncomplementary_edges = 0;

  pd_idx_t tpos;
  for(tpos = f0_tanglepos; /* We could be counting down in mod arithmetic */
      t->face[tpos] != border_faces[n-1];
      tpos = (tpos + step + t->nedges)%t->nedges) /* Move around face buffer */ { 

    complementary_edges[ncomplementary_edges] = t->edge[tpos];
    complementary_or[ncomplementary_edges] = t->edge_bdy_or[tpos];
    ncomplementary_edges++;

  }

  /* We now have a list of 0 or more complementary edges to work
     with. We're going to split each one in the middle with a new
     crossing, and add a new edge.

     Here, we remember that pd_working->MAXVERTS has already
     been increased to give us room, but we'll still check the size of 
     the buffer as we go in case that code is somehow screwy. */

  for(i=0;i<ncomplementary_edges;i++) { 

    if (pd_working->ncross > pd_working->MAXVERTS-1) { 

      pd_error(SRCLOC,
	       "we should have increased MAXVERTS (=%d) in pd %PD to make room for\n"
	       "ncomplementary_edges = %d new crossings, but we're out of\n"
	       "room after only %d new crossings\n",pd_working,pd_working->MAXVERTS,
	       ncomplementary_edges,i);
      exit(1);

    }

    if (pd_working->nedges > pd_working->MAXEDGES-1) { 

      pd_error(SRCLOC,
	       "we should have increased MAXEDGES (=%d) in pd %PD to make room for\n"
	       "ncomplementary_edges = %d new edges, but we're out of\n"
	       "room after only %d new edges\n",pd_working,pd_working->MAXEDGES,
	       ncomplementary_edges,i);
      exit(1);

    }

    /* We now have a picture like:
 						       |
	    |tail crossing			       |tail crossing
       	----|-----                                 ----|----	     
       	    |				     new_edge->|
       	    |<-complementary_edges[i] 		       0        	     	
	    |                             ---->	       +<-new_cross    	       	   
       .....|.....(tangle boundary)	      .........2..........(tangle boundary)
	    V	 			               |<-complementary_edges[i]    
      	----------				       V       	      
	    |head crossing			   ---------          
						       | head crossing
							 
      or 						 
							 
      	    |head crossing			       |head crossing
       	----|-----                                 ----|---- 	     
       	    ^	 			     new_edge->^     
       	    |<-complementary_edges[i] 		       0       	     		
            |                             ---->	       +<-new_cross    	       	   
       .....|.....(tangle boundary)	      .........2..........(tangle boundary)
	    |	 			               |       	       	       	
	----------				       |<-complementary_edges[i]   
            |tail crossing			   ---------          
		       	       	       	       	       | tail crossing
							
     We are choosing to do the crossing inserts in this
     orientation-sensitive way in order to make sure that the new
     edges are consistenly OUTSIDE the original tangle.  This is going
     to help us keep the crossings oriented consistently as we sew in
     the overstrand.
						       	
    */

    pd_idx_t new_cross, new_edge;
    new_cross = pd_working->ncross; 
    new_edge = pd_working->nedges;

    pd_working->cross[new_cross].edge[0] = new_edge;
    pd_working->cross[new_cross].edge[1] = PD_UNSET_EDGE;
    pd_working->cross[new_cross].edge[2] = complementary_edges[i];
    pd_working->cross[new_cross].edge[3] = PD_UNSET_EDGE;

    if (complementary_bdy_or == in) { 

      pd_working->edge[new_edge].head = new_cross;
      pd_working->edge[new_edge].headpos = 0;
      pd_working->edge[new_edge].tail = pd_working->edge[complementary_edges[i]].head;
      pd_working->edge[new_edge].tailpos = pd_working->edge[complementary_edges[i]].headpos;

    } else if (complementary_bdy_or == out) { 

      pd_working->edge[new_edge].head = pd_working->edge[complementary_edges[i]].head;
      pd_working->edge[new_edge].headpos = pd_working->edge[complementary_edges[i]].headpos;
      pd_working->edge[new_edge].tail = new_cross;
      pd_working->edge[new_edge].tailpos = 0;

    } else { 

      pd_error(SRCLOC,"complementary_bdy_or[%d] = %d is neither 'in' nor 'out'\n",pd_working,
	       i,(int)(complementary_bdy_or[i]));
      exit(1);

    }

    pd_working->ncross++;
    pd_working->nedges++;

  }

  /* We now deal with a special case: there are NO complementary edges. */

  if (ncomplementary_edges == 0) { 

    /* We have the following picture: 
     		       	        	    
       	     +-------------------+		  
	     |		    	 |		  
 end_anchor  | 	    Tangle     	 |      	       	  
       |     | 	       	       	 |    +---start_anchor   	    
       |     | 	    	    	 |    | 	     
       |     +-------------------+    |	     
       |         |   |       |        |
       |         |   |       | 	      |	     
       +--e[n-1]+|  +|       |+--e[0]-+	      
       	       	 |   | 	     |		   
  			    
    		     | 
	       	     | 	       	     
	       	     V	    	     
      +-------------------e[0]--------+
      |	     +-------------------+    |
      +      |        	      	 |    ^
 end_anchor  | 	    Tangle       |    |
             | 	       	       	 |    +---start_anchor   	    
             | 	    	      	 |      	     
             +-------------------+     	     
  e[n-1]         |   |       |         
 deleted         |   |       | 	       	     	       	
       	     +   |  +|       |  +      	       	       	
 former tail   	 |   | 	     |	  former head of e[0]  
   of e[n-1]   	      				       
    
    */

    pd_edge_t e0,enm1;
    e0 = pd->working->edge[overstrand_edges[0]];
    enm1 = pd->working->edge[overstrand_edges[n-1]];

    /* Disconnect from the former head and tail crossings. */
    pd_working->cross[e0.head].edge[e0.headpos] = PD_UNSET_IDX;
    pd_working->cross[enm1.tail].edge[enm1.tailpos] = PD_UNSET_IDX;

    /* Connect e0 to the end_anchor crossing record. */
    pd_working->cross[end_anchor].edge[enm1.headpos] = overstrand_edges[0];
    
    /* Change the head and headpos of e0 to join the edge to end_anchor. */
    pd_working->edge[overstrand_edges[0]].head = end_anchor;
    pd_working->edge[overstrand_edges[0]].headpos = enm1.headpos;

    /* Delete the edge e[n-1] */
    pd_working->edge[overstand_edges[n-1]] = pd_unset_edge(); 

  } else {

       /* We have the following picture: 

		  |   |	   |  
		  +   +	   +  <- new crossings, one for each complementary edge 
     		  |   |	   |   
       	     +-------------------+		  
	     |		       	 |		  
 end_anchor  | 	    Tangle     	 |      	       	  
       |     | 	       	       	 |    +---start_anchor   	    
       |     | 	    	       	 |    | 	     
       |     +-------------------+    |	     
       |         |   |       |        |
       |         |   |       | 	      |	     
       +--e[n-1]+|  +|       |+--e[0]-+	      
       	       	 |   | 	     | 		   
  			       
    		     | 	       
	       	     V 	       	     
	       	       	       
		|    |	      |
    e[n-1]	|    |	      |
      +---------+----+--------+-------+
      |         |    |        |       |e[0]
      |	     +-------------------+    |
      +      |        	      	 |    ^
 end_anchor  | 	    Tangle       |    |
             | 	       	       	 |    +---start_anchor   	    
             | 	    	      	 |      	     
             +-------------------+     	     
                 |   |       |         
                 |   |       | 	       	     	       	
       	     +   |  +|       |  +      	       	       	
 former tail   	 |   | 	     |	  former head of e[0]  
   of e[n-1]   	      				       
    	  
    */	  
      	  
    pd_edge_t e0,enm1;
    e0 = pd->working->edge[overstrand_edges[0]];
    enm1 = pd->working->edge[overstrand_edges[n-1]];

    /* We still start by disconnecting from the former head and tail crossings. */
    pd_working->cross[e0.head].edge[e0.headpos] = PD_UNSET_IDX;
    pd_working->cross[enm1.tail].edge[enm1.tailpos] = PD_UNSET_IDX;

    /* Now we need to connect e[0] and e[n-1] correctly. We already know 
       that the complementary edges are going in the order in which the 
       new strand will encounter them along the boundary-- we just need to
       know whether this order is clockwise or counterclockwise. 

       This was determined above by the variable "step", which is +1 if we're
       going counterclockwise (positively oriented, same as the direction of 
       tangle_edges and the like) and -1 if we're going clockwise. 

       In the first case, we're coming IN at position 3 on the new crossings
       and leaving at position 1, the second, vice-versa.*/

    pd_pos_t enter_pos, exit_pos;
    if (step == +1) { enter_pos = 1; exit_pos = 3; }
    else if (step == -1) { enter_pos = 3; exit_pos = 1; }
    else { pd_error(SRCLOC,"step (%d) is neither +1 nor -1. Something is corrupted.\n",
		    pd_working); exit(1); }

    /* Change the head and headpos of e[0] to join the edge to first new crossing. */
    pd_working->edge[overstrand_edges[0]].head    = pd->ncross;
    pd_working->edge[overstrand_edges[0]].headpos = enter_pos;

    /* Change the tail and tailpos of e[n-1] to join the edge to last new crossing.*/
    pd_working->edge[overstrand_edges[n-1]].tail = pd->ncross + (ncomplementary_edges-1);
    pd_working->edge[overstrand_edges[n-1]].tailpos = exit_pos;

    /* Now add new edges joining new crossings in pairs. */

    for(i=1;i<ncomplementary_edges;i++) { 

      /* Check to make sure we're ok on memory */

      if (pd_working->nedges > pd_working->MAXEDGES-1) { 

	pd_error(SRCLOC,
		 "should have room to add another %d edges to %PD\n"
		 "but we already have %d and MAXEDGES is %d\n",
		 pd_working,ncomplementary_edges-i,pd->nedges,pd->MAXEDGES);
	exit(1);

      }

      /* Now add the new edge. */
      
      pd_idx_t new_edge = pd->nedges;

      pd_working->edge[new_edge].head = pd->ncross+i;
      pd_working->edge[new_edge].headpos = enter_pos;
      pd_working->edge[new_edge].tail = pd->ncross+i-1;
      pd_working->edge[new_edge].tailpos = exit_pos;

      pd_working->nedges++;

    }

  }

  /* Step 6. Cleanup. At this point, we should have valid crossing data for a pd
     code, and edges giving orientations and so forth. However, there are several
     problems.

     #1. There are "dead" crossings sprinkled throughout pd_working->cross, and 
     "dead" edges sprinkled throughout pd_working->edges.

     #2. The new crossings don't have sign information. 

     #3. Neither edges nor crossings are ordered. 

     #4. Components don't make sense-- some components may consist
     ONLY of dead edges and crossings, which others may simply have
     some dead edges and components.

     #5. The diagram may be disconnected. 

     #6. Faces need to be regenerated. */

}
