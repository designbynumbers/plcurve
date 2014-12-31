/*              _________________________________                             */
/*             /                                 \                            */
/*    (6)     /         _________________________________3>_______            */
/*           /         /                           3              \           */
/*          /   +-----/---------------------------|---------+      \          */
/*         /    |    /                            6         |      |          */
/*        /     |   /             _____9>________ v ___     |      |          */
/*        |     |  /             /                |0   \    |      |          */
/*        |     | /             /                 |     \   |<-T   |          */
/*        |     | |            /        (7)       |      \  |      |          */
/*        |     | |           /                   7       \ |      |          */
/*        |     | |          /                    v       | |      |          */
/*        |     | |          |      __<10__       |       | | (0)  |          */
/*        |     | |     (3)  |     /       \      |       0 |      |          */
/*        |     | |          |    /   (9)   \     |       v |      |          */
/*        |     | |          |    |7        |6    |  (4)  | |      |          */
/*        |     | |          \_________<8___|_____/       | |      |          */
/*        |     | |                         |             / |      |          */
/*        |     | \               |         ^            /  |      |          */
/*        |     |  \              11   (2) 15           /   |      |          */
/*        |     |   \             v         |          /    |      |          */
/*        |(5)  |    \_______<2__ | ___<1_____________/     |      |          */
/*        \     |                 |2         1              |      /          */
/*         \    |                 12  (1)   |               |     /           */
/*          \   +-----------------|---------14--------------+    /            */
/*           \____________<5___________<4__ | __________________/             */
/*                                 5        |4                                */
/*                                |  (8)    |                                 */
/*                                \         /                                 */
/*                                 \_13>___/                                  */
pd_idx_t nedges = 4 ;
pd_idx_t tangle_faces[4] = {0,3,5,1} ;
pd_idx_t tangle_edges[4] = {6,2,12,14} ;
pd_boundary_or_t edge_bdy_or[4] = {in,out,out,in} ;

pd_idx_t ninterior_cross = 5 ;
pd_idx_t interior_cross[5] = {0,1,2,6,7} ;

pd_idx_t ninterior_edges = 8 ;
pd_idx_t interior_edge[8] = {0,1,7,8,9,10,11,15} ;

pd_idx_t nstrands = 2;
pd_tangle_strand_t strand_data[2] = {{0,1,7,0},
				     {3,2,5,1}};
