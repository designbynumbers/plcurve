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
    
      pd_printf("edge t->edge[%d] = %EDGE is not on face t->face[%d] %FACE exactly once: tangle test fails\n",pd,i,t->edge[i],i,t->face[i]);
      return false; 

    }

    if (!pd_xor(pos_face[i] == t->face[(i+1)%t->nedges],neg_face[i] == t->face[(i+1)%t->nedges])) { 

      pd_printf("edge t->e[%d] = %EDGE is not on face t->face[%d] %FACE exactly once: tangle test fails\n",pd,i,t->edge[i],(i+1)%(t->nedges),t->face[(i+1)%t->nedges]);
      return false; 

    }
 
  }

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
     Let's check the strand data itself for each strand. */

  for(i=0;i<t->nstrands;i++) { 

    pd_check_edge(SRCLOC,pd,t->strand[i].start_edge);
    pd_check_edge(SRCLOC,pd,t->strand[i].end_edge);
    pd_check_cmp(SRCLOC,pd,t->strand[i].comp);

    if (t->strand[i].nedges == PD_UNSET_IDX) { 

      pd_error(SRCLOC,"strand[%d] of tangle has number of edges == PD_UNSET_IDX (maybe tangle not regenerated?)\n",
	       pd,i);
      return false;

    }

    if (t->strand[i].nedges > pd->comp[t->strand[i].comp].nedges) {

      pd_error(SRCLOC,"strand[%d] of tangle claims to have %d edges on component %COMP, which only has %d edges\n",
	       pd,i,t->strand[i].nedges,t->strand[i].comp,pd->comp[t->strand[i].comp].nedges);
      return false;

    }

  }

  /* We're now going to try to test-walk each strand to verify that the
     number of edges and the start and end data is ok. */

  for(i=0;i<t->nstrands;i++) { 

    pd_idx_t edge;
    for(edge = t->strand[i].start_edge, j=0; 
	edge != t->strand[i].end_edge && j != t->strand[i].nedges;
	j++,edge = pd_next_edge(pd,edge)) { }

    if (edge != t->strand[i].end_edge) { 

      pd_error(SRCLOC,
	       "strand[%d] of tangle claims to have %d edges\n"
	       "but counting %d edges from start edge %EDGE\n"
	       "leaves us at %EDGE != the end_edge of %EDGE\n",
	       pd,i,t->strand[i].nedges,j,t->strand[i].start_edge,
	       edge,t->strand[i].end_edge);

      return false;

    }

    if (j != t->strand[i].nedges) { 

      pd_error(SRCLOC,
	       "strand[%d] of tangle claims to have %d edges\n"
	       "but counting edges from start_edge %EDGE\n"
	       "finds the end_edge %EDGE after %d edges\n",
	       pd,i,t->strand[i].nedges,t->strand[i].start_edge,
	       t->strand[i].end_edge,j);

      return false;

    }

  }

  /* Now we check that comp makes sense. */

  for(i=0;i<t->nstrands;i++) { 

    pd_idx_t comp,pos;
    pd_component_and_pos(pd,t->strand[i].start_edge,&comp,&pos);
    
    if (t->strand[i].comp != comp) { 

      pd_error(SRCLOC,
	       "tangle strand %d claims to be on component %COMP, but this\n"
	       "component does not include start_edge %EDGE\n",pd,i,
	       t->strand[i].comp,t->strand[i].start_edge);
      return false;

    }

    pd_component_and_pos(pd,t->strand[i].end_edge,&comp,&pos);
    
    if (t->strand[i].comp != comp) { 

      pd_error(SRCLOC,
	       "tangle strand %d claims to be on component %COMP, but this\n"
	       "component does not include end_edge %EDGE\n",pd,i,
	       t->strand[i].comp,t->strand[i].end_edge);
      return false;

    }

  }
 
  /* There's no really practical way to check that the interior crossings
     are, in fact, interior (same with the edges) without actually regenerating
     the data here as part of the check, which seems like too much redundant 
     code. 

     We can check that they are allocated, however, and in fact that they refer 
     to valid crossing and edge numbers, and we do this in the spirit of checking
     for data corruption wherever it might crop up. */
     
  pd_check_notnull(SRCLOC,"t->interior_cross",t->interior_cross);
  pd_check_notnull(SRCLOC,"t->interior_edge",t->interior_edge);

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

  if (t->ninterior_edges > pd->ncross) { 

    pd_error(SRCLOC,"t->ninterior_edges = %d > # of crossings in entire pd code (%d)\n",pd,
	     t->ninterior_edges,pd->ncross);
    return false;

  }
  
  /* Now we run through the buffers and check them: */

  for(i=0;i<t->ninterior_cross;i++) { 

    pd_check_cr(SRCLOC,pd,t->interior_cross[i]);
    
  }

  for(i=0;i<t->ninterior_edges;i++) { 

    pd_check_cr(SRCLOC,pd,t->interior_edge[i]);

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

      pd_error(SRCLOC,"edge %EDGE and face %FACE are paired on tangle, but this edge is not on this face in %PD",pd,t->edge[i],t->face[i]);
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
      pd_component_and_pos(pd,t->edge[i],&(t->strand[s].comp),&pos);

      pd_idx_t j;
      bool found_end;

      for(j=pd_next_edge(pd,t->strand[s].start_edge),found_end = false;
	  !found_end && j != t->strand[s].start_edge; // failsafe
	  j=pd_next_edge(pd,j),
	    t->strand[s].nedges++) {

	/* Now we have to search the list of edges of the tangle 
	   to try to find this particular edge j. */

	pd_idx_t bdypos;

	if (pd_buf_contains(t->edge,t->nedges,j,&bdypos)) { 

	  if (!(t->edge_bdy_or[bdypos] == out)) { /* Just to be safe */
	    
	    pd_error(SRCLOC,"tangle strand which started (in) with %EDGE "
		     "has encountered tangle edge %d (%EDGE) with boundary "
		     "orientation not 'out' \n",pd,
		     t->strand[s].start_edge,bdypos,t->edge[bdypos]);
	    exit(1);

	  }

	  found_end = true;
	  t->strand[s].end_edge = j;

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

  
