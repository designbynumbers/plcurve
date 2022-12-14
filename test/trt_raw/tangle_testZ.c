/*                                                          */
/*              /----C--7---<-----\          (1)            */
/*           c1 |       (5)       | c3                      */
/*      /-----------<---3---<----------2-\                  */
/*      |       |                 |      |                  */
/*      |   T   6                 6      |                  */
/*      |       v                 ^      |                  */
/*      v   +------------------------+   ^                  */
/*      |   |   |       (4)       |  |   0                  */
/*      |   |   \-->--C--6--->----/  |   A                  */
/*     e0   |     (0)          (0)   |   0                  */
/*      |   |                        |   2                  */
/*      |   |   /----<-B-5---<----\  |   0                  */
/*      A   |   |                 |  |   |                  */
/*      |   |   5       (3)       5  |   ^                  */
/*      v   |   v                 ^  |   0                  */
/*      0   +------------------------+   0                  */
/*      0    c0 |                 | c2   |                  */
/*      \--0--000000->--0-1--00->-0000-2-/                  */
/*              |                 |                         */
/*              4       (2)       4         (1)             */
/*              \-->-B---4---->---/                         */
/*                                                          */
 						
pd_idx_t nedges = 4;	      			
pd_idx_t tangle_faces[4] = {4,0,3,0};	
pd_idx_t tangle_edges[4] = {6,5,5,6};	

pd_boundary_or_t edge_bdy_or[4] = {in,out,in,out};
pd_idx_t ninterior_cross = 0;
pd_idx_t interior_cross[1] = {PD_UNSET_IDX};

pd_idx_t ninterior_edges = 0;
pd_idx_t interior_edge[1] = {PD_UNSET_IDX};

pd_idx_t nstrands = 2;
pd_tangle_strand_t strand_data[2] = {{0,3,1,2},
				  {2,1,1,1}};
