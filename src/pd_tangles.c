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
  pd_check_notnull(SRCLOC,pd,"input pd code");
  pd_check_notnull(SRCLOC,t,"input tangle");
  
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
	     t->nstrands,t->nedges,t->nedges/2);
    return false;

  }

  pd_check_notnull(SRCLOC,t->edge,"t->edge");
  pd_check_notnull(SRCLOC,t->face,"t->face");

  pd_idx_t *pos_face,*pfp,*neg_face,*nfp;
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

      pd_printf("edge t->e[%d] = %EDGE is not on face t->face[%d] %FACE exactly once: tangle test fails\n",pd,i,t->edge[i],(i+1)%(t->nedges),t->f[(i+1)%t->nedges]);
      return false; 

    }
 
  }

  /* We now check edge_bdy_or */

  pd_check_notnull(SRCLOC,t->edge_bdy_or,"t->edge_bdy_or");
  pd_idx_t in_count = 0, out_count = 0;

  for(i=0;i<t->nedges;i++) { 
   
    if (edge_bdy_or[i] == in) { 
      
      if (in_count >= t->nstrands) { 

	pd_error(SRCLOC,"tangle has %d strands, but more than %d flow into the tangle\n",
		 pd,t->nstrands,in_count);
	return false;
      }

      in_count++; 

    } else if (edge_bdy_or[i] == out) {

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
    pd_check_comp(SRCLOC,pd,t->strand[i].comp);

    if (t->strand[i].nedges == PD_UNSET_IDX) { 

      pd_error(SRCLOC,"strand[%d] of tangle has number of edges == PD_UNSET_IDX (maybe tangle not regenerated?)\n",
	       pd,i);
      return false;

    }

    if (t->strand[i].nedges > pd->comps[t->strand[i].comp].nedges) {

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
     
  pd_check_notnull(SRCLOC,t->interior_cross,"t->interior_cross");
  pd_check_notnull(SRCLOC,t->interior_edge,"t->interior_edge");

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

  for(i=0;i<t->ninterior_edge;i++) { 

    pd_check_cr(SRCLOC,pd,t->interior_edge[i]);

  }

  /* We've checked everything that we can! */

  return true;
   
}
	       
void pd_regenerate_tangle() 

{
 pd_idx_t in_count = 0, out_count = 0;

  for(i=0;i<4;i++) { 
    if (edge_bdy_or[i] == in) { 
      
      if (in_count >= 2) { 
	pd_error(SRCLOC,"edges e[0] = %EDGE, e[1] = %EDGE, e[2] = %EDGE, e[3] = %EDGE have incorrect orientations to be part of tangle\n",
		 pd,e[0],e[1],e[2],e[3]);
	exit(1);
      }

      in_edges[in_count] = i;
      in_count++; 

    } else {

      if (out_count >= 2) { 
	pd_error(SRCLOC,"edges e[0] = %EDGE, e[1] = %EDGE, e[2] = %EDGE, e[3] = %EDGE have incorrect orientations to be part of tangle\n",
		 pd,e[0],e[1],e[2],e[3]);
	exit(1);
      }

      out_edges[out_count] = i;
      out_count++; 
    } 
  }

}

 void pd_tangle_bdy_orientation(pd_code_t *pd, pd_idx_t e[4], pd_idx_t f[4],
				pd_boundary_or_t edge_bdy_or[4],
				pd_idx_t in_edges[2], pd_idx_t out_edges[2])

  /* We are now going to start by deducing "innie/outie" orientations
     for each of the four boundary edges of the tangle. */
{
  pd_idx_t i;

  for(i=0;i<4;i++) { 
    pdx_idx_t pos_face,pfp,neg_face,nfp;
    pd_face_and_pos(pd,e[i],&pos_face,&pfp,&neg_face,&nfp);
    if (pos_face == f[i]) { 
      edge_bdy_or[i] = in;
    } else if (neg_face == f[i]) { 
      edge_bdy_or[i] = out;
    } else {
      pd_error(SRCLOC,"edge e[%d] = %EDGE in tangle is not on face f[%d] = %FACE in pd\n",
	       pd,i,e[i],i,f[i]);
      exit(1);
    }
  }

 
}

    
void pd_tangle_contents(pd_code_t *pd, pd_idx_t e[4], pd_idx_t f[4],
			pd_idx_t *ncross,pd_idx_t **tangle_cross,  
			pd_idx_t *nedges,pd_idx_t **tangle_edge)
    
{
  assert(pd_tangle_ok(pd,e,f));
  
  pd_boundary_or_t edge_bdy_or[4];
  pd_idx_t in_edges[2], out_edges[2];

  pd_tangle_bdy_orientation(pd,e,f,edge_bdy_or,
			    in_edges,out_edges);

  /* Now we need to construct the two strands of the tangle. We'll call these
     strand[0] and strand[1] and observe that each has a start_edge and an end_edge. */



  /* The first thing we do is make a pass along each component 
     until we find an end edge. */

  bool found_end;

  for(i=0;i<2;i++) { 
    
    strand[i].start_edge = in_edges[i];
    strand[i].num_edges = 1;

    for(j=pd_next_edge(pd,strand[i].start_edge),found_end = false;
	j != strand[i].start_edge;
	j=pd_next_edge(pd,j),
	strand[i].num_edges++) {

      if (j == out_edges[0] || j == out_edges[1]) { 

	strand[i].end_edge = j;
        found_end = true;

      }
    
    }

    if (!found_end) { pd_error(SRCLOC,"couldn't find either out edge %EDGE or %EDGE along component containing in edge %EDGE in %PD",
			       pd,out_edges[0],out_edges[1],in_edges[i]); exit(1); }
    
  }

  /* Now that we have isolated both strands of the tangle, we need to fill 
     in the remainder of the tangle. Roughly, the problem is this: we have
     counted the edges in each strand of the tangle. But other components 
     might cross the these strands (in which case, the entire component 
     has to be added to the tangle's total). */
  
  pd_idx_t pos;
  for(i=0;i<2;i++) { pd_component_and_pos(pd,strand[i].start_edge,&(strand[i].comp),&pos); }

  /* So now we're going to rewalk each strand, from start to end, marking components for 
     inclusion in the tangle. */

  bool *component_in_tangle = calloc(pd->ncomps,sizeof(bool));
  assert(component_in_tangle != NULL); for(i=0;i<pd->ncomps;i++) { component_in_tangle = false; }

  for(i=0;i<2;i++) { 

    for(j=strand[i].start_edge;
	j != strand[i].end_edge;
	j = pd_next_edge(pd,j)) { 

      /* At each crossing in the tangle, look at the OTHER strand through the crossing, and 
	 add that component (if it's not either of the strand components) */

      pd_idx_t comp, pos;
      pd_component_and_pos(pd,pd->cross[pd->edge[j].head].edge[(pd->edge[j].headpos+1)%4],&comp,&pos);
      if (comp != strand[i].comp && comp != strand[1-i].comp) { component_in_tangle[comp] = true; }

    }
    
  }
      
  /* We now have a list of components which are entirely contained in
     the tangle.  Our job is to figure out a list of edges in the
     tangle and a list of crossings, then allocate space for them
     and copy them into the output buffers. 
  */

  bool *edge_in_tangle = calloc(pd->nedges,sizeof(bool));
  assert(edge_in_tangle != NULL); for(i=0;i<pd->nedges;i++) { edge_in_tangle = false; }

  bool *crossing_in_tangle = calloc(pd->ncross,sizeof(bool));
  assert(crossing_in_tangle != NULL); for(i=0;i<pd->ncross;i++) { crossing_in_tangle = false; }
  
  /* We first walk the strands again, adding crossings and edges as we go. */

  for(i=0;i<2;i++) { 

    for(j=strand[i].start_edge;
	j != strand[i].end_edge;
	j = pd_next_edge(pd,j)) { 

      edge_in_tangle[j] = true;
      crossing_in_tangle[pd->edge[j].head] = true;

    }
    
  }
     
  for(i=0;i<4;i++) { edge_in_tangle[e[i]] = true; }

  /* Now we walk the list of components adding entire components as we go. */

  for(i=0;i<pd->ncomps;i++) { 

    if (component_in_tangle[i]) { 

      for(j=0;j<pd->comp[i].nedges;j++) { 

	edge_in_tangle[j] = true;
	crossing_in_tangle[pd->edge[j].head] = true;

      }

    }

  }

  /* Finally, we count the crossings and edges that we've found and copy them 
     to output buffers. */

  for(i=0,*ncross=0;i<pd->ncross;i++) { if (crossing_in_tangle[i]) { (*ncross)++; } }
  for(i=0,*nedges=0;i<pd->nedges;i++) { if (edge_in_tangle[i]) { (*nedges)++; } }

  (*tangle_cross) = calloc(*ncross,sizeof(pd_idx_t)); assert(tangle_cross != NULL);
  for(j=0,i=0;i<pd->ncross && j < *ncross;i++) { if (crossing_in_tangle[i]) { (*tangle_cross)[j++] = i; } }

  (*tangle_edge) = calloc(*nedges,sizeof(pd_idx_t)); assert(tangle_edge != NULL);
  for(j=0,i=0;i<pd->nedges && j < *nedge;i++) { if (edge_in_tangle[i]) { (*tangle_edge)[j++] = i; } }

  /* Last, we do some housekeeping */

  free(component_in_tangle);
  free(crossing_in_tangle);
  free(edge_in_tangle);

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
  pd_idx_t i,j;
  assert(pd != NULL);
  for(i=0;i<4;i++) { pd_check_edge(SRCLOC,pd,e[i]); pd_check_face(SRCLOC,pd,f[i]); }

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
  
  pd_idx_t ntangle_cross, *tangle_cross;
  pd_idx_t ntangle_edge, *tangle_edge;

  pd_tangle_contents(pd,e,f,
		     &ntangle_cross,&tangle_cross,
		     &ntangle_edge,&tangle_edge);

  /* We'll need access to the orientation data, too. */

  pd_boundary_or_t edge_bdy_or[4];
  pd_idx_t in_edges[2], out_edges[2];

  pd_tangle_bdy_orientation(pd,e,f,edge_bdy_or,
			    in_edges,out_edges);

  /* We can find the crossing cr by taking the end of the 
     edge e[0] inside the tangle. */

  if (edge_bdy_or[0] == in) { cr = pd->edge[e[0]].head; }
  else { cr = pd->edge[e[0]].tail; }

  /* But it's wise to check that this agrees with the same
     test done on e[1]. */

  if (edge_bdy_or[1] == in) { 
    assert(cr == pd->edge[e[1]].head); 
  } else { 
    assert(cr == pd->edge[e[1]].tail); 
  }
      
  /* We also should check that this crossing is listed in-tangle. */    

  bool found_in_tangle = false;
  for(j=0;j<ntangle_cross;j++) { if (tangle_cross[j] == cr) { found_in_tangle = true; } }

  if (!found_in_tangle) { 
    
    pd_error(SRCLOC,"edges e[0] = %EDGE and e[1] = %EDGE share crossing %CROSS, but this isn't listed in-tangle\n"
	     "in %PD.\n",
	     pd,e[0],e[1],cr);
    exit(1);

  }

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

  
