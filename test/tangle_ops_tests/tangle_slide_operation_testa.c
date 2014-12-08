/*                                                           */
/*                                                           */
/*                         +---------e6----------+           */
/*                         o                     |           */
/*              (2)        o <-overstrand        |           */
/*                         o                     |           */
/*   T-> +------------+    v                     |           */
/*       |            |    |c3                   |           */
/*       |    +-------|------------+             |           */
/*       |    |       |    o       |             |           */
/*       |   e4       |    e7     e3             |           */
/*       |    |   (3) |    v  (5)  |             |           */
/*       |  c2|       |    o       |             |           */
/*    +--|----|--->---|-e2------>--+             |           */
/*    |  |    v       |    oc0                   |           */
/*   e1  +------------+    v           (1)       |           */
/*    | (4)  e5    (0)     o                     |           */
/*    +---<---------e0--ooo+                     |           */
/*          c1|                                  |           */
/*            +----->-------e6-------------------+           */

pd_idx_t nedges = 4;
pd_idx_t tangle_faces[4] = {0,3,2,4};
pd_idx_t tangle_edges[4] = {2,4,1,5};

pd_idx_t noverstrand_edges = 3;
pd_idx_t overstrand_edges[3] = {6,7,0};
pd_idx_t border_faces[3] = {2,3,0};

/*                                                           */
/*                                                           */
/*     +-------2--------+                                    */
/*     |                |                                    */
/*     |      (4)       v      (0)                           */
/*     |                |                                    */
/*     2                |c0                                  */
/*     |      +---5--------->---6-----+                      */
/*     |      ^   (2)   3     (1)     |                      */
/*     |    c2|         |c0           |c1                    */
/*     +--<---|----1-<---------0--<---|----7----+            */
/*            |         |             |         |            */
/*            ^         v             v         |            */
/*            4   (3)   |             |   (5)   7            */
/*            |         4             7         |            */
/*            |         |             |         |            */
/*            +---------+             +---------+            */
/*                                                           */
/*                                                           */
/*                                                           */
