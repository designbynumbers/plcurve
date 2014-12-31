/*           +------------------------------+                                 */
/*           |                              | <-T                             */
/*        ____<11_________      _______16>_____________________________       */
/*       /   |            \    /            |                          \      */
/*      /   _|__<2____     \  /             |                           \     */
/*     /   / |        \     \ 8             |              (11)          \    */
/*    /    | |         \    /\              |                            |    */
/*   /     | | (3)      \ 15  \             |                            |    */
/*  /      | |     _     /2    \____________|________<10_____            |    */
/*  |      | |    / \  14 \                 |                \           |    */
/*  | (7)  | |  19   \ /   \       _______<5|___              \          |    */
/*  |      | |  v     \10   ^     /         |   \              \         |    */
/*  |      | | / (12)  \     1   /          |    \              \        |    */
/*  |      | ||     /   \( 2) \ /           |     \              \       |    */
/*  |      | | \  13     \     \ 1      (9 )|     |    (1)        \      |    */
/*  |      | |  \ /       ^   | \           |     |                \     |    */
/*  \      | |   / 9      18  |  \          |     |                |     /    */
/*   \____ |_12>/          \  6   \         |     |                |    /     */
/*        3| |     \        \ V    \___<0___|____ | _______<17_________/      */
/*         | |      \  (6)   \|             |     |0                7         */
/*          \|       \        |5            |     |                |          */
/*           \       20       |\      (0)   |     |                |          */
/*           |\         v     7 \________   |     |                |          */
/*           | \         \    v (10)     \  |     |                |          */
/*           |  \         \    \         |  |     /                |          */
/*           |   \         \_________21>_/  |    /                 |          */
/*       (5) |    \   (4)         6         |   /    (8)           |          */
/*           |     \              \         |  /                   |          */
/*           |      \              8        | /                    |          */
/*           |       \_3>__________ \ ___4>__/                     |          */
/*           +-----------------------\4-----+                      /          */
/*                                    \___________9>______________/           */

pd_idx_t nedges = 10 ;
pd_idx_t tangle_faces[10] = {0,9,1,11,5,7,3,4,5,8} ;
pd_idx_t tangle_edges[10]        = {0 ,5 ,10,16 , 11,2  ,12,3 ,9  ,4} ;
pd_boundary_or_t edge_bdy_or[10] = {in,in,in,out,out,out,in,in,out,out} ;
/*                                  0   1  2  3    4   5  6  7   8   9 */

pd_idx_t ninterior_cross = 8 ;
pd_idx_t interior_cross[8] = {1,2,4,5,6,8,9,10} ;

pd_idx_t ninterior_edges = 11 ;
pd_idx_t interior_edge[11] = {1,6,7,8,13,14,15,18,19,20,21} ;

pd_idx_t nstrands = 5;
pd_tangle_strand_t strand_data[5] = {{0,5,3,0},
				     {1,8,5,0},
				     {2,4,2,0},
				     {6,3,5,0},
				     {7,9,2,0}};