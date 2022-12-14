/*                                _____<1___                                  */
/*                               /          \                                 */
/*                               |           \                                */
/*                               |    (9)     \                               */
/*                               |            |                               */
/*                  ___<11_____________<10____|_______________                */
/*      (0)        /             |2           |1              \               */
/*                /              |            |                \              */
/*               /               2            /                 \             */
/*               |     (4)       v    (3)    /                   \            */
/*               |               |          /                     \           */
/*               |    +----------|---------/------+ (1)            ^          */
/*               \    |      17>_|____0>__/       |                 9         */
/*                \   |     /    3 0              |                  \        */
/*                 \  |    / (5) v                |<-T                \       */
/*                  \_____/__12>_____________13>_________________     |       */
/*                    |  / 8     4 3              |              \    |       */
/*             ___16>___/(2)     v     (6)        |               \   |       */
/*            /       +----------|----------------+               |   /       */
/*       0000/00000007>0000000000|000000000008>0000000000000000000|00/        */
/*      /   /6                   |4                               |7          */
/*     /   /                     |                                /           */
/*    /   /                      5                               /            */
/*    |   |          (7)         v            (8)               /             */
/*    |   |                      |                             /              */
/*    |   \                      |                            /               */
/*    |    \_______<15____________________________<14________/                */
/*    |                          /5                                           */
/*    |         (10)            /                                             */
/*    \___________________<6___/                                              */

pd_idx_t nedges = 6 ;            
pd_idx_t tangle_faces[6] = {0,2,6,1,3,4} ;
pd_idx_t tangle_edges[6] = {16,4,13,0,2,11} ;

pd_idx_t noverstrand_edges = 4;
pd_idx_t overstrand_edges[4] = {4,7,8,9};
pd_idx_t border_faces[4] = {0,2,6,1};

bool valid_ts_input = false ; /* incorrect input overstrand edges */

pd_idx_t tangle_slide_edges[1] = {PD_UNSET_IDX};
pd_idx_t complementary_edges[1] = {PD_UNSET_IDX};
pd_boundary_or_t complementary_or[1] = {unset};
bool overstrand_goes_OVER = false;
pd_or_t overstrand_orientation = PD_UNSET_ORIENTATION;
