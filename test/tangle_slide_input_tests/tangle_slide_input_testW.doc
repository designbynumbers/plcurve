/*                                                         */
/*                                                         */
/*                                                         */
/*                                                         */
/*             /-----0-->--\                               */
/*             ^           |                               */
/*             |c1-        |c0+                            */
/*     /--9--------<-8----------<----7-----------\         */
/*     |       13          |                     ^         */
/*     |       |c9+        1                     |c7+      */
/*     |  /----------16->-------------17--->----------\    */
/*     v  |    ^    (5)    v c2-                 |    |    */
/*     |  ^  +-------------+------+<-tangle      |    |    */
/*     |  |  | |   /--<-2--/      |              |    |    */
/*     |  |  | 12  |              |              |    |    */
/*     |  |  | |   |c3+           |              |    |    */
/*     9  |  | \-------\          | (1)          6    18   */
/*     |  |  | (6) |(7)|          |              |    |    */
/*     |  ^  +--------------------+              |    v    */
/*     |  |        3  11                         ^    |    */
/*     |  |     c4+v  ^c8-                       |    |    */
/*     |  \---15----14----<-------19--\          |    |    */
/*     |           |  10              \-----<---------/    */
/*     |           4>--|-->\                     |c6-      */
/*     v           c5- |   |                     |         */
/*     |               |   5                     |         */
/*     \------>--9-----/   \-------->---5--------/         */
/*                                                         */
/*                                                         */
/*                                                         */
		       	
pd_idx_t nedges = 4;            
pd_idx_t tangle_faces[4] = {5 ,6,7 ,1};
pd_idx_t tangle_edges[4] = {12,3,11,2};

pd_idx_t noverstrand_edges = 3;
pd_idx_t overstrand_edges[3] = {14,15,16};
pd_idx_t border_faces[3]     = {7,6,5};

bool valid_ts_input = true; 

pd_idx_t tangle_slide_edges[2] = {3,12};
pd_idx_t complementary_edges[3] = {11,2};
pd_boundary_or_t complementary_or[3] = {in,out};
bool overstrand_goes_OVER = true;
pd_or_t overstrand_orientation = PD_NEG_ORIENTATION;
