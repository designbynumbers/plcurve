/*            ____________________________________<4______________            */
/*    (1)    /                                                    \           */
/*          /             +--------------------------------------- \ -------+ */
/*         /              |                      (0)                \___    | */
/*        /               |                                             \   | */
/*        |               |                                              \  | */
/*        |            T->|             __6>___                          |  | */
/*        |               |            /  (9)  \                         |  | */
/*        \    __________9>___________/__10>_______11>_                  |  | */
/*         \  /           |   (8)    / 6         7     \    _2>_   _13_  |  | */
/*           / _______5>____________/           |       \  /    \ /    \ |  | */
/*          / 5           |           (4)       |         / (6)  /  (7)  /  | */
/*          \             |                     /        /      /3\     /   | */
/*           \______________________<8____     7        /2\_12_/   \3>_/4\  | */
/*                        |               \   v        /                 |  | */
/*                        |                \ /        /                  |  | */
/*        ________________|                 /1       /                   |  | */
/*       /                |\    _____<0____/  \     /                    |  | */
/*      /                 | \  /              |    /                     |  | */
/*     /                  |  \         (5)    |   /                      |  | */
/*    /                   |   \               |  /                       |  | */
/*   /                    | / 0\___15>________/ /             (3)        |  | */
/*  /                     |/                   /                         |  | */
/*  |                     /           (2)     /                          |  | */
/*  |                    /+------------------/---------------------------|--+ */
/*  |                    |                  /                            |    */
/*  |                    \_______1>________/                             /    */
/*  |                                                                   /     */
/*  \                                                                  /      */
/*   \_________________________________________<14____________________/       */
pd_idx_t nedges = 8 ;
pd_idx_t tangle_faces[8] = {0,8,4,1,3,2,3,1} ;
pd_idx_t        tangle_edges[8] = {9 ,5 ,8  ,14,1  ,1 ,14 ,4} ;
pd_boundary_or_t edge_bdy_or[8] = {in,in,out,in,out,in,out,out} ;
/*                                 0  1  2   3  4   5  6   7    */
pd_idx_t ninterior_cross = 7  ;
pd_idx_t interior_cross[7]  = {0,1,2,3,4,6,7} ;

pd_idx_t ninterior_edges = 10 ;
pd_idx_t interior_edge[10] = {0,2,3,6,7,10,11,12,13,15} ;

pd_idx_t nstrands = 4;
pd_tangle_strand_t strand_data[4] = {{0,6,6,0},
				     {1,5,5,0},
				     {3,2,3,1},
				     {5,7,4,0}};
