/*
 * This procedure closes all open components of plCurve by distributing the
 * necessary change over all the vertices of each such component. It also
 * removes the now duplicate vertex, changes the "open" flag and calls
 * fix_wrap.
 *
 */
void plc_force_closed(plCurve * const L) {
  int i, cmp, nv;
  plc_vector diff;
  double half;
  plc_constraint *pcst,*cst;

  for (cmp=0;cmp < L->nc;cmp++) {
    if (L->cp[cmp].open == true) {  /* Isolate the open components. */
      nv = L->cp[cmp].nv;
      half = (nv-1.0)/2.0;

      /* Compute the error in closure */
      diff = plc_vect_diff(L->cp[cmp].vt[nv-1],L->cp[cmp].vt[0]);

      for (i=0; i < nv; i++) {
        /* add half of diff to vt[0] and subtract half of diff from vt[0] and
         * prorate the others :-) */
        plc_M_vmadd(L->cp[cmp].vt[i],(half-i)/(nv-1.0),diff);
      }

      /* We claim to have moved the last vertex on top of the first. */
      diff = plc_vect_diff(L->cp[cmp].vt[0], L->cp[cmp].vt[nv-1]);
      assert(plc_M_dot(diff,diff) < DBL_EPSILON);

      /* Thus we eliminate the last vertex. */
      L->cp[cmp].nv--;
      L->cp[cmp].open = false;

      /* And perhaps the last color. */
      L->cp[cmp].cc = intmin(L->cp[cmp].cc,L->cp[cmp].nv);

      pcst = NULL;
      cst = L->cst;
      while (cst != NULL && (cst->cmp < cmp ||
          (cst->cmp == cmp && cst->vert + cst->num_verts <= L->cp[cmp].nv))) {
        pcst = cst;
        cst = cst->next;
      }
      if (cst != NULL && cst->cmp == cmp) {
        cst->num_verts = L->cp[cmp].nv - cst->vert;
        if (cst->num_verts == 0) {
          /*@-kepttrans@*/
          if (pcst != NULL) {
            /*@-mustfreeonly@*/
            pcst->next = cst->next;
            /*@=mustfreeonly@*/
            free(cst);
            cst = NULL;
          } else {
            /*@-nullderef@*/
            L->cst = cst->next;
            /*@=nullderef@*/
            free(cst);
            cst = NULL;
          /*@-branchstate@*/
          }
          /*@=kepttrans@*/
        }
      }
    }
    /*@=branchstate@*/
  }
  
  /* Someday we'll deal with quantifiers here */

  plc_fix_wrap(L);
/*@-compdef -usereleased@*/
} /* plc_force_closed */
/*@=compdef =usereleased@*/
