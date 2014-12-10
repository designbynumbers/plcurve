/*                                                         */
/*                                                         */
/*                   /--3--\                               */
/*                   A     ^                               */
/*                   | (6) |c0                             */
/*                   \-3>--|-->-\          (0)             */
/*                         2    0                          */
/*                      c2 | (2)|c1                        */
/*                /----7-----<6----5----\                  */
/*                |        |    |       |                  */
/*                |  T     ^ (3)v       |                  */
/*                B  +---------------+  ^                  */
/*                |  |     \--1-/    |  |                  */
/*                |  |       (1)     |  |                  */
/*                7  |               |  |                  */
/*                |  |  /-C--9----\  |  |                  */
/*                v  +---------------+  |                  */
/*                |     v   (5)   ^     |                  */
/*                |     |         |     |                  */
/*                \7---------4-->-----5-/                  */
/*                   c4 |         |c3                      */
/*                      C   (4)   |                        */
/*                      \----8-->-/                        */
		       
pd_idx_t nedges = 4 ;            
pd_idx_t tangle_faces[4] = {3,1,5,1};
pd_idx_t tangle_edges[4] = {1,9,9,1};

pd_idx_t noverstrand_edges = 3;
pd_idx_t overstrand_edges[3] = {7,4,5};
pd_idx_t border_faces[3]     = {1,5,1};

bool valid_ts_input = true ; 

pd_idx_t tangle_slide_edges[2] = {9,9};
pd_idx_t complementary_edges[3] = {1,1};
pd_boundary_or_t complementary_or[3] = {out,in};
bool overstrand_goes_OVER = true;
pd_or_t overstrand_orientation = PD_POS_ORIENTATION;
