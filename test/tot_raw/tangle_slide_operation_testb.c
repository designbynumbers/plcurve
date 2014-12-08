/*                                                          */
/*                                        (1)               */
/*                      /->--0-A--\                         */
/*                      |c1+ (4)  |c0-                      */
/*         /000-11-0->-------8-->----0000-B-000\            */
/*         0            ^         v            |            */
/*         0   (0)      7   (2)   1     (0)    B            */
/*         0   T->+-------------------+        |            */
/*         0      |     |c4-      vc2-|        0            */
/*         0   /--|->3--|----4->------|-5-\    0            */
/*         0   |  |     |   (3)   |   |   |    |            */
/*         0   ^  |     \--6--\   2   |   v    V            */
/*         0   3  |(6)        ^   v  (5)  5    |            */
/*         0   |  |      /<-A-|---/   |   |    0            */
/*         |   A  |      | c3+|       |   A    0            */
/*         11  \--|---3--/    \--<-5--|---/    |            */
/*         |      |       (0)         |        9            */
/*         0   /--|---C------13-------|---\    |            */
/*         ^   |  |                   |   |    0            */
/*         0   v  +-------------------+   ^    v            */
/*         0   |          (8)             |    0            */
/*         \000000-B-000000-10-000000-<-0000000/            */
/*          c6-|           (7)            |c5-              */
/*             \--------C----12---->------/                 */
/*                                                          */
/*    (1)                                                   */
       	       	  
pd_idx_t nedges = 8;
pd_idx_t tangle_faces[8] = {2,0,6,0,8,0,5,0};
pd_idx_t tangle_edges[8] = {7,3,3,13,5,5,1};

pd_idx_t noverstrand_edges = 3;
pd_idx_t overstrand_edges[3] = {9,10,11};
pd_idx_t border_faces[3] = {0,8,0};

/*                                                          */
/*                                                          */
/*                     /-A--0-->--\                         */
/*     	       	   c0+ |       	  | c1-                     */
/*            /--11--------12>--B--------\                  */
/*            B        9          1      B                  */
/*            ^    c7- |          |      v                  */
/*            \--11--------10<--B-----13-/                  */
/*                     ^          2 c2+                     */
/*                     8          v                         */
/*            /---4->--|--->-5--A----->6---\                */
/*            |    c5- |          | c3-    |                */
/*            |        \--7-<--\  3        v                */
/*            ^                |  v        |                */
/*            |            /---|--/        6                */
/*            4            |   | c4+       |                */
/*            |            v   |           A                */
/*            A            4   \--<--6-----/                */
/*            |            |                                */
/*            \-----A------/                                */
/*                                                          */
/*                                                          */
/*                                                          */
/*                                                          */
/*                                                          */
/*                                                          */
/*                                                          */
