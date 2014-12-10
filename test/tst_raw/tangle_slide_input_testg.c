/*                     __________20>___________________________               */
/*                    /                                        \              */
/*                   /                                          \             */
/*         (0)      /                             (1)           |             */
/*                 /     +--------------------+                 |             */
/*                 |     |                    |                 /             */
/*              __ | _1>_____                 |                /              */
/*             /   |1    |   \             ______________________             */
/*            /    |     |    \           ^   |                9 \            */
/*           /     \     |(11) \         9    |             /     \           */
/*          /       \    |      \        |    |   (9)      /       \          */
/*         /         \____<19__ | __<18_____  |           /        |          */
/*         |             |      |2       |8 \ |         21         |          */
/*         |    (6)      |      2        ^   ^|         v          /          */
/*         |             |      v  (2)   8    17       /          /           */
/*         \             |      |        |(8) | 0     /          /            */
/*          \_______________<0_____<23__ | _<22_ 0 __/   (10)  10             */
/*                       |      |0       ^ 7  |  011          v               */
/*                       |      3   (3)  7 (5)|  ^           /                */
/*                       |      v        |    |  16         /                 */
/*          ______________<13__ | __<12____<11_____________/                  */
/*         /             |      |3       |6   |  0 10                         */
/*        /              |      4   (4)  ^    |  0                            */
/*        \     (12)     +------v--------6----+  0                            */
/*         \                    |        |  (7)  0                            */
/*          \0000000000000000000|0014>000000015>0                             */
/*                              |4       |5                                   */
/*                              \  (13)  /                                    */
/*                               \__5>__/                                     */
pd_idx_t nedges = 10 ;            
pd_idx_t tangle_faces[10] = {0,12,4,7,5,8,9,1,11,6} ;
pd_idx_t tangle_edges[10] = {13,4,6,11,22,17,9,1,19,0} ;

pd_idx_t noverstrand_edges = 5;
pd_idx_t overstrand_edges[5] = {13,14,15,16,17};
pd_idx_t border_faces[5] = {12,4,7,5,8};

bool valid_ts_input = false ; /* overstrand alternates under/over */

pd_idx_t tangle_slide_edges[1] = {PD_UNSET_IDX};
pd_idx_t complementary_edges[1] = {PD_UNSET_IDX};
pd_boundary_or_t complementary_or[1] = {unset};
bool overstrand_goes_OVER = false;
pd_or_t overstrand_orientation = PD_UNSET_ORIENTATION;
