/* 

   test_isomorphism.c : Unit tests for the code in pd_isomorphisms.c.


*/

#ifdef HAVE_CONFIG_H
  #include"config.h"
#endif

#ifdef HAVE_STDIO_H
   #include<stdio.h>
#endif 

#ifdef HAVE_STRING_H
   #include<string.h>
#endif

#ifdef HAVE_STDINT_H
   #include<stdint.h>
#endif

#ifdef HAVE_STDLIB_H
   #include<stdlib.h>
#endif

#ifdef HAVE_GSL_GSL_PERMUTATION_H
   #include<gsl/gsl_permutation.h>
#endif

#ifdef HAVE_STDBOOL_H
   #include<stdbool.h>
#endif

#ifdef HAVE_ASSERT_H
   #include<assert.h>
#endif

#include<plcTopology.h>

#include<pd_multidx.h>
#include<pd_perm.h>
#include<pd_dihedral.h>

#include<pd_isomorphisms.h>

#include"argtable2.h" /* We use a local copy of argtable */

int PD_VERBOSE=50;

bool test_compgrps()

{

  printf("pd_build_compgrps test suite\n"
	 "-------------------------------\n");

  pd_code_t pdA,pdB;

  printf("different number of components...");

  pdA.ncomps = 5;
  pdB.ncomps = 4;

  pd_compgrp_t *grps = NULL;
  pd_idx_t     ngrps = 255;

  grps = pd_build_compgrps(&pdA,&pdB,&ngrps);

  if (grps == NULL && ngrps == 0) {

    printf("pass.\n");

  } else {

    printf("FAIL (grps,ngrps == %p, %d instead of NULL, 0).\n",grps,ngrps);
    return false;

  }
  
  printf("edge vectors (2,3,4) and (2,3,5) don't match...");

  pdA.ncomps = 3;
  pdA.comp[0].nedges = 2; pdA.comp[1].nedges = 3; pdA.comp[2].nedges = 4;
  
  pdB.ncomps = 3;
  pdB.comp[0].nedges = 2; pdB.comp[1].nedges = 3; pdB.comp[2].nedges = 5;

  grps = pd_build_compgrps(&pdA,&pdB,&ngrps);

  if (grps == NULL && ngrps == 0) {

    printf("pass.\n");

  } else {

    printf("FAIL (grps,ngrps == %p, %d instead of NULL, 0).\n",grps,ngrps);
    return false;

  }
  
  printf("\n"
	 "specific group tests\n"
	 "***********************\n"
	 "nedges vec (2 3 5) -> ");

  pdA.ncomps = 3;
  pdA.comp[0].nedges = 2; pdA.comp[1].nedges = 3; pdA.comp[2].nedges = 5;

  pdB.ncomps = 3;
  pdB.comp[0].nedges = 2; pdB.comp[1].nedges = 3; pdB.comp[2].nedges = 5;
 
  grps = pd_build_compgrps(&pdA,&pdB,&ngrps);
  
  if (grps == NULL || ngrps != 3) {

    printf("FAIL. ngrps = %d != 3.\n",ngrps);
    return false;

  }

  if (grps[0].ncomps != 1 || grps[0].comp[0] != 0) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 0 ) \n",&pdA,&(grps[0]));
    return false;

  } 

  if (grps[1].ncomps != 1 || grps[1].comp[0] != 1) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 1 ) \n",&pdA,&(grps[1]));
    return false;

  } 

  if (grps[2].ncomps != 1 || grps[2].comp[0] != 2) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 2 ) \n",&pdA,&(grps[2]));
    return false;

  } 

  pd_printf("%COMPGRP   %COMPGRP  %COMPGRP    pass.\n",&pdA,&(grps[0]),&(grps[1]),&(grps[2]));
  free(grps);

  printf("nedges vec (3 3 5) -> ");

  pdA.ncomps = 3;
  pdA.comp[0].nedges = 3; pdA.comp[1].nedges = 3; pdA.comp[2].nedges = 5;

  pdB.ncomps = 3;
  pdB.comp[0].nedges = 3; pdB.comp[1].nedges = 3; pdB.comp[2].nedges = 5;
 
  grps = pd_build_compgrps(&pdA,&pdB,&ngrps);
  
  if (grps == NULL || ngrps != 2) {

    printf("FAIL. ngrps = %d != 2.\n",ngrps);
    return false;

  }

  if (grps[0].ncomps != 2 || grps[0].comp[0] != 0 || grps[0].comp[1] != 1) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 0 1 ) \n",&pdA,&(grps[0]));
    return false;

  } 

  if (grps[1].ncomps != 1 || grps[1].comp[0] != 2) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 2 ) \n",&pdA,&(grps[1]));
    return false;

  } 

  pd_printf("%COMPGRP %COMPGRP                   pass.\n",&pdA,&(grps[0]),&(grps[1]));
  free(grps);

  printf("nedges vec (3 5 5) -> ");

  pdA.ncomps = 3;
  pdA.comp[0].nedges = 3; pdA.comp[1].nedges = 5; pdA.comp[2].nedges = 5;

  pdB.ncomps = 3;
  pdB.comp[0].nedges = 3; pdB.comp[1].nedges = 5; pdB.comp[2].nedges = 5;
 
  grps = pd_build_compgrps(&pdA,&pdB,&ngrps);
  
  if (grps == NULL || ngrps != 2) {

    printf("FAIL. ngrps = %d != 2.\n",ngrps);
    return false;

  }

  if (grps[0].ncomps != 1 || grps[0].comp[0] != 0) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 0 ) \n",&pdA,&(grps[0]));
    return false;

  } 

  if (grps[1].ncomps != 2 || grps[1].comp[0] != 1 || grps[1].comp[1] != 2) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 1 2 ) \n",&pdA,&(grps[1]));
    return false;

  }  

  pd_printf("%COMPGRP   %COMPGRP                 pass.\n",&pdA,&(grps[0]),&(grps[1]));
  free(grps);

  printf("nedges vec (5 5 5) -> ");

  pdA.ncomps = 3;
  pdA.comp[0].nedges = 5; pdA.comp[1].nedges = 5; pdA.comp[2].nedges = 5;

  pdB.ncomps = 3;
  pdB.comp[0].nedges = 5; pdB.comp[1].nedges = 5; pdB.comp[2].nedges = 5;
 
  grps = pd_build_compgrps(&pdA,&pdB,&ngrps);
  
  if (grps == NULL || ngrps != 1) {

    printf("FAIL. ngrps = %d != 1.\n",ngrps);
    return false;

  }

  if (grps[0].ncomps != 3 || grps[0].comp[0] != 0 || grps[0].comp[1] != 1 || grps[0].comp[2] != 2) {

    pd_printf("FAIL. group 0 = %COMPGRP != ( 0 1 2 ) \n",&pdA,&(grps[0]));
    return false;

  } 

  pd_printf("%COMPGRP                               pass.\n",&pdA,&(grps[0]));
  free(grps);

  printf("-------------------------------------------\n"
	 "pd_build_compgrps test suite        PASS\n\n");

  return true;

}

bool evec_test(pd_idx_t ncomps,pd_idx_t *evec,bool printperms,bool printemaps)

{
  if (printemaps) { printperms = true; }

  if (ncomps > PD_MAXCOMPONENTS) { ncomps = PD_MAXCOMPONENTS; }

  pd_code_t pdA, pdB;
  pd_idx_t comp, compedge,edge=0;

  printf("Checking nedges vec ...                  ( ");

  pdA.ncomps = pdB.ncomps = ncomps;
  pdA.nedges = pdB.nedges = 0;
  
  for(comp=0;comp<ncomps;comp++) {
    
    pdA.comp[comp].nedges = evec[comp];
    pdB.comp[comp].nedges = evec[comp];
    
    for(compedge=0;compedge < evec[comp];compedge++,edge++) {
      
      pdA.comp[comp].edge[compedge] = edge;
      pdB.comp[comp].edge[compedge] = edge;
      
    }
    
    pdA.nedges += evec[comp];
    pdB.nedges += evec[comp];
    
    printf("%d ",evec[comp]);

  }

  printf(")\n");

  /* First, we need to build comp_grps */

  pd_compgrp_t *compgrps;
  pd_idx_t      ngrps,grp;

  
  printf("Building compgrps ...                   ");
  
  compgrps = pd_build_compgrps(&pdA,&pdB,&ngrps);
  for(grp=0;grp<ngrps;grp++) { 

    if (pd_compgrp_ok(&pdA,&(compgrps[grp])) == false) {
      
      printf("FAIL. grp %d doesn't pass pd_compgrp_ok.\n",grp);
      return false;

    }

  }

  printf(" pass (groups ok)\n");
  printf("Groups are ... \n");

  for(grp=0;grp<ngrps;grp++) { 

    pd_printf("\t %COMPGRP \n",NULL,&(compgrps[grp])); 

  }

  printf("\n");

  printf("Building comp_perms from these groups... ");

  pd_perm_t **comp_perms;
  unsigned int i,ncomp_perms;
  
  comp_perms = pd_build_compperms(ngrps,compgrps,&ncomp_perms);

  printf("\n");

  for(i=0;i<ncomp_perms;i++) {

    if (printperms) {

      pd_printf("\t %PERM ... ",NULL,comp_perms[i]);
      if (pd_perm_ok(comp_perms[i])) { printf("ok\n"); }
      else {printf("FAIL (not ok).\n"); return false; }

    } else {

      printf("\t %d \n",i);

    }

    pd_edgemap_t **edgemap_buf;
    unsigned int nedgemaps;
    
    edgemap_buf = pd_build_edgemaps(&pdA,&pdB,comp_perms[i],&nedgemaps);
      
    printf("\t\t %d edgemaps built ... ",nedgemaps);

    if (printemaps) { 

      printf("\n");
      
      unsigned int j;
      for(j=0;j<nedgemaps;j++) {

	pd_printf("\t\t %EDGEMAP ... ",NULL,edgemap_buf[j]);
	if (pd_edgemap_ok(edgemap_buf[j])) { printf("ok.\n"); }
	else { printf("FAIL (not ok).\n"); return false; }

      }
    
    }
    
    if (!pd_edgemaps_ok(nedgemaps,edgemap_buf)) {

      printf("FAIL. (edgemaps list not ok or contains dups).\n");
      return false;

    }

    unsigned int j;
    for(j=0;j<nedgemaps;j++) {

      if (!pd_edgemap_consistent(&pdA,&pdB,edgemap_buf[j])) {

	pd_printf("FAIL. %EDGEMAP %d not consistent with pdA, pdB.\n",NULL,
		  edgemap_buf[j],j);
	return false;

      }

    }

    printf("pass (edgemaps ok && consistent).\n");

    pd_free_edgemaps(nedgemaps,&edgemap_buf);
    pd_free_edgemaps(nedgemaps,&edgemap_buf);

  }

  printf("\n"
	 "pd_compperms_ok   ...                    ");


  if (!pd_compperms_ok(ncomp_perms,comp_perms)) {

    printf("FAIL. List does not pass (pd_compperms_ok).\n");
    return false;

  } else {

    printf("pass.\n");

    printf("                                         ");

    printf("perms are ok and unique.\n");

  }

  printf("Checking # perms generated ...           ");

  unsigned int ncomp_perms_expected=1;
  
  for(grp=0;grp<ngrps;grp++) { 

    for(comp=1;comp<=compgrps[grp].ncomps;comp++) {

      ncomp_perms_expected *= comp;

    }

  }

  if (ncomp_perms != ncomp_perms_expected) { 

    printf("FAIL. # comp_perms %d != expected value %d.\n",ncomp_perms,ncomp_perms_expected);
    return false;

  } else {

    printf("pass. (%d == %d).\n",ncomp_perms,ncomp_perms_expected);

  }

  printf("edgevec ( ");

  for(comp=0;comp<ncomps;comp++) {

    printf("%d ",evec[comp]);

  }

  printf(") pd_build_compperms tests pass.\n\n");

  free(compgrps);

  pd_free_compperms(ncomp_perms,&comp_perms);
  pd_free_compperms(ncomp_perms,&comp_perms);

  return true;

}  

bool test_compperms_and_edgemaps() {

  int store_verb = PD_VERBOSE;
  PD_VERBOSE = 0;

  printf("pd_build_compperms and pd_build_edgemaps test suite\n"
	 "---------------------------------------------------\n");

  pd_idx_t evec[10];

  evec[0] = 2; evec[1] = 3; evec[2] = 3; evec[3] = 5;
  if (!evec_test(4,evec,true,false)) { return false; } 

  evec[0] = 2; evec[1] = 3; evec[2] = 3; evec[3] = 3;
  if (!evec_test(4,evec,true,false)) { return false; }

  evec[0] = 3; evec[1] = 3; evec[2] = 3; evec[3] = 3;
  if (!evec_test(4,evec,true,false)) { return false; } 

  /* evec[0] = 3; evec[1] = 3; evec[2] = 3; evec[3] = 3; evec[4] = 3; */
  /* if (!evec_test(5,evec,false,false)) { return false; }  */ 

  /* evec[0] = 2; evec[1] = 2; evec[2] = 2;  */
  /* evec[3] = 2; evec[4] = 2; evec[5] = 2; */
  /* evec[6] = 2; evec[7] = 2; */
  
  /* if (!evec_test(7,evec,false,false)) { return false; }   */

  printf("--------------------------------------------------------\n"
	 "pd_build_compperms and pd_build_edgemaps test suite PASS\n\n");

  PD_VERBOSE = store_verb;

  return true;

}

bool unknot_iso_test(pd_idx_t ncross,bool print)

{
  pd_code_t *pdpA, *pdpB;
  pd_code_t pdA,pdB;
  
  printf("%d-crossing unknot automorphism group test  \n",ncross);
  printf("--------------------------------------------\n");

  pdpA = pd_build_unknot(ncross);
  pdpB = pd_build_unknot(ncross);

  printf("pd is ");
  pd_write(stdout,pdpA);

  pdA = *pdpB; pdB = *pdpB;

  printf("\t building (auto)morphisms...");

  unsigned int nisos;
  pd_iso_t **isos = pd_build_isos(pdpA,pdpB,&nisos);

  if (print) {

    unsigned int i;

    printf("\n");

    for(i=0;i<nisos;i++) {

      pd_printf("\t %ISO \n",NULL,isos[i]);

    }

    printf("\t checking # isos generated...");

  }

  if (nisos != 4) {

    printf("FAIL. (# isos generated == %d != 4).\n",nisos);
    return false;

  } else {

    printf("pass (# isos generated == %d == 4).\n",nisos);

  }

  printf("\t checking isomorphism effect on plane...");

  pd_idx_t plane_neg[4];
  pd_idx_t plane_pos[4];
  pd_idx_t iso,npn=0,npp=0;

  for(iso=0;iso<nisos;iso++) {

    if (isos[iso]->crossmap->or == PD_NEG_ORIENTATION) {

      plane_neg[npn++] = iso;

    } else {

      plane_pos[npp++] = iso;
      
    }

  }

  if (npn != 2 || npp != 2) {

    printf("FAIL. (# (plane) or reversing, # (plane) or preserving) == (%d,%d) != (2,2).\n",
	   npn,npp);
    return false;

  } else {

    printf("pass. (# (plane) or -, # (plane) or +) == (%d,%d) == (2,2).\n",
	   npn,npp);
  
  }

  /* The n-crossing unknot should have 4 automorphisms: e, rot, hflip, and vflip. */
  
  pd_idx_t e, rot, hflip, vflip;

  printf("\t checking for e, rot, hflip, vflip...");

  if (ncross % 2 != 0) { /* EVEN # of crossings. */

    if (isos[plane_pos[0]]->edgemap->or[0] == isos[plane_pos[1]]->edgemap->or[0]) { /* Same orientation */

      printf("FAIL.\n(# crossings == %d == odd, but BOTH plane+ isos (rot and e) are component %c)",
	     ncross,pd_print_or(isos[plane_pos[0]]->edgemap->or[0]));
      return false;

    }

    if (isos[plane_pos[0]]->edgemap->or[0] == PD_POS_ORIENTATION &&
	isos[plane_pos[1]]->edgemap->or[0] == PD_NEG_ORIENTATION) { 

      e = plane_pos[0]; 
      rot = plane_pos[1];

    } else {

      e = plane_pos[1];
      rot = plane_pos[0];

    }

    if (isos[plane_neg[0]]->edgemap->or[0] == isos[plane_neg[1]]->edgemap->or[1]) {

       printf("FAIL.\n(# crossings == %d == odd, but BOTH plane- isos (hflip, vflip) are component %c)",
	     ncross,pd_print_or(isos[plane_neg[0]]->edgemap->or[0]));
      return false;

    }

    if (isos[plane_neg[0]]->edgemap->or[0] == PD_POS_ORIENTATION &&
	isos[plane_neg[1]]->edgemap->or[0] == PD_NEG_ORIENTATION) { 

      hflip = plane_neg[0]; 
      vflip = plane_neg[1];

    } else {

      hflip = plane_neg[1]; 
      vflip = plane_neg[0];

    }

  } else {  /* number of crossings is odd */

    if (isos[plane_pos[0]]->edgemap->or[0] != isos[plane_pos[1]]->edgemap->or[0]) { /* diff orientation */

      printf("FAIL. (# crossings == %d == odd, but plane+ isos are component %c and %c)",
	     ncross,
	     pd_print_or(isos[plane_pos[0]]->edgemap->or[0]),
	     pd_print_or(isos[plane_pos[1]]->edgemap->or[0]));
      return false;

    }

    /* We now know the plane+ isomorphisms should match e
       and rot.  Since these are component+ as well, to
       tell them apart, we will check their effect on the
       1-edge faces. */

    assert(pdpA->face[pdpA->nfaces-1].nedges == 1 && pdpA->face[pdpA->nfaces-2].nedges == 1);

    bool swaps_endfaces[2];

    for(iso=0;iso<2;iso++) { 

      swaps_endfaces[iso] = (isos[plane_pos[iso]]->facemap->perm->map[pdpA->nfaces-1] == pdpA->nfaces-2 &&
			     isos[plane_pos[iso]]->facemap->perm->map[pdpA->nfaces-2] == pdpA->nfaces-1);

    }

    if (swaps_endfaces[0] == swaps_endfaces[1]) { 

      printf("FAIL. (# crossings == %d == odd, but both plane+ maps",ncross);

      if (swaps_endfaces[0]) { 

	printf("swap endfaces.)\n");

      } else {

	printf("DON'T swap endfaces.)\n");

      }

      return false;

    }

    if (swaps_endfaces[0] && !swaps_endfaces[1]) { 

      e = plane_pos[1];
      rot = plane_pos[0];

    } else if (swaps_endfaces[1] && !swaps_endfaces[0]) { 

      e = plane_pos[0];
      rot = plane_pos[1];

    }

    /* We have now found both e and rot. We try to identify hflip and vflip. */ 

    for(iso=0;iso<2;iso++) { 

      swaps_endfaces[iso] = (isos[plane_neg[iso]]->facemap->perm->map[pdpA->nfaces-1] == pdpA->nfaces-2 &&
			     isos[plane_neg[iso]]->facemap->perm->map[pdpA->nfaces-2] == pdpA->nfaces-1);

    }

    if (swaps_endfaces[0] == swaps_endfaces[1]) { 

      printf("FAIL. (# crossings == %d == odd, but both plane- maps",ncross);

      if (swaps_endfaces[0]) { 

	printf("swap endfaces.)\n");

      } else {

	printf("DON'T swap endfaces.)\n");

      }

      return false;

    }

    if (swaps_endfaces[0] && !swaps_endfaces[1]) { 

      vflip = plane_neg[1];
      hflip = plane_neg[0];

    } else if (swaps_endfaces[1] && !swaps_endfaces[0]) { 

      vflip = plane_neg[0];
      hflip = plane_neg[1];

    } else {

      printf("FAIL. (couldn't identify hflip and vflip)\n");
      return false;

    }

  }

  printf("pass. (tentative id: (e,rot,hflip,vflip) = isos (%d,%d,%d,%d))\n",
	 e,rot,hflip,vflip);

  printf("\t checking that iso e == isos[%d] is identity... ",e);

  if (!pd_iso_is_e(isos[e])) { 

    printf("FAIL. (isos[e] == isos[%d] != identity).",e);
    return false;

  } 

  printf("pass (isos[e] == identity).\n");

  printf("\t checking that rot, hflip, vflip have period 2... ");

  if (pd_iso_period(isos[rot]) != 2) {

    printf("FAIL. (isos[rot] == isos[%d] has period %d != 2).\n",
	   rot,pd_iso_period(isos[rot]));
    return false;

  } 

  if (pd_iso_period(isos[hflip]) != 2) {

    printf("FAIL. (isos[hflip] == isos[%d] has period %d != 2).\n",
	   hflip,pd_iso_period(isos[hflip]));
    return false;

  } 

  if (pd_iso_period(isos[vflip]) != 2) {

    printf("FAIL. (isos[vflip] == isos[%d] has period %d != 2).\n",
	   vflip,pd_iso_period(isos[vflip]));
    return false;

  } 

  printf("pass.\n");

  printf("\t checking that hflip * vflip == rot ...");

  pd_iso_t *hv = pd_compose_isos(isos[hflip],isos[vflip]);
  
  if (pd_iso_cmp(&hv,&(isos[rot])) != 0) { 

    printf("FAIL.\n");
    return false;

  } 

  pd_free_iso(&hv);
  printf("pass.\n");

  printf("\t checking that hflip * rot == vflip ...");

  pd_iso_t *hr = pd_compose_isos(isos[hflip],isos[rot]);
  
  if (pd_iso_cmp(&hr,&(isos[vflip])) != 0) { 

    printf("FAIL.\n");
    return false;

  } 

  pd_free_iso(&hr);
  printf("pass.\n");

  printf("\t checking that vflip * rot == hflip ...");

  pd_iso_t *vr = pd_compose_isos(isos[vflip],isos[rot]);
  
  if (pd_iso_cmp(&vr,&(isos[hflip])) != 0) { 

    printf("FAIL.\n");
    return false;

  } 

  pd_free_iso(&vr);
  printf("pass.\n");

  printf("\t free data test... ");

  pd_free_isos(&nisos,&isos);
  pd_free_isos(&nisos,&isos);
  free(pdpA); free(pdpB);

  printf("pass.\n");
 
  printf("-------------------------------------\n");
  printf("%d-crossing unknot tests ... PASS \n\n",ncross);

  return true;

}  

bool twist_iso_test(pd_idx_t ntwist,bool print)

{
  pd_code_t *pdpA, *pdpB;
  pd_code_t pdA,pdB;

  if (ntwist < 3) {

    printf("FAIL. n-twist knot test not implemented for n < 3.\n");
    return false;

  }
  
  printf("%2d-twist knot automorphism group test \n",ntwist);
  printf("--------------------------------------------\n");
  
  pdpA = pd_build_twist_knot(ntwist);
  pdpB = pd_build_twist_knot(ntwist);
  
  printf("pd is ");
  pd_write(stdout,pdpA);
  
  pdA = *pdpB; pdB = *pdpB;
  
  printf("\t building (auto)morphisms...");

  unsigned int nisos;
  pd_iso_t **isos = pd_build_isos(pdpA,pdpB,&nisos);

  if (print) {

    unsigned int i;

    printf("\n");

    for(i=0;i<nisos;i++) {

      pd_printf("\t %ISO \n",NULL,isos[i]);

    }

    printf("\t checking # isos generated...");

  }
  
  if (nisos != 4) {
    
    printf("FAIL. (# isos generated == %d != 4).\n",nisos);
    return false;
    
  } else {

    printf("pass (# isos generated == %d == 4).\n",nisos);
    
  }
  
  printf("\t identifying 3-edge face preserving isos...");
  
  /* Since the faces are reverse-sorted, these should be faces 2 and 3 */

  pd_idx_t nfound=0;
  pd_idx_t evflip[4];
  pd_idx_t i;

  for(i=0;i<nisos;i++) {

    if (isos[i]->facemap->perm->map[2] == 2 &&
	isos[i]->facemap->perm->map[3] == 3) {

      evflip[nfound++] = i;

    }

  }

  if (nfound != 2) { 

    printf("FAIL. Found %d != 2 isomorphisms which preserve the 3 edge faces.\n",nfound);
    return false;

  }

  printf("pass. (isos %d and %d preserve the 3 edge faces)\n",evflip[0],evflip[1]);
  
  printf("\t distinguishing e from vflip by effect on orientation...");

  for(i=0;i<2;i++) { 

    if (isos[evflip[i]]->edgemap->or[0] != isos[evflip[i]]->crossmap->or) { 

      printf("FAIL.\n isos[%d] doesn't preserve or reverse BOTH plane and curve orientation.\n",evflip[i]);
      return false;

    }

  }

  if (isos[evflip[0]]->crossmap->or == isos[evflip[1]]->crossmap->or) {

    printf("FAIL.\n isos[%d] AND isos[%d] preserve or reverse plane orientation.\n",evflip[0],evflip[1]);
    return false;

  }

  pd_idx_t e, vflip;

  if (isos[evflip[0]]->crossmap->or == PD_POS_ORIENTATION) { 

    e = evflip[0]; vflip = evflip[1];
    printf("pass. (e == isos[%d], vflip == isos[%d])\n",evflip[0],evflip[1]);

  } else {

    e = evflip[1]; vflip = evflip[0];
    printf("pass. (e == isos[%d], vflip == isos[%d])\n",evflip[1],evflip[0]);

  }

  printf("\t checking that iso e == isos[%d] is identity... ",e);
  
  if (!pd_iso_is_e(isos[e])) { 
    
    printf("FAIL. (isos[e] == isos[%d] != identity).",e);
    return false;
    
  } 
  
  printf("pass (isos[e] == identity).\n");
  
  printf("\t checking that vflip has period 2... ");
  
  if (pd_iso_period(isos[vflip]) != 2) {
    
    printf("FAIL. (isos[vflip] == isos[%d] has period %d != 2).\n",
	   vflip,pd_iso_period(isos[vflip]));
    return false;
    
  } 
  
  printf("pass.\n");

  pd_idx_t hflip;
  pd_idx_t hbuf[4];

  printf("\t looking for hflip, which should preserve %d edge faces, swap 3 edge faces...",ntwist+1);

  for(i=0,nfound=0;i<nisos;i++) {

    if (isos[i]->facemap->perm->map[0] == 0 &&
	isos[i]->facemap->perm->map[1] == 1 &&
	isos[i]->facemap->perm->map[2] == 3 &&
	isos[i]->facemap->perm->map[3] == 2) {

      hbuf[nfound++] = i;

    }

  }

  if (nfound != 1) {

    printf("FAIL. Found more than one iso with this property.\n");
    return false;

  }

  hflip = hbuf[0];

  printf("pass. (hflip = isos[%d])\n",hflip);

  printf("\t checking that hflip has period 2...");

  if (pd_iso_period(isos[hflip]) != 2) {
    
    printf("FAIL. (isos[hflip] == isos[%d] has period %d != 2).\n",
	   hflip,pd_iso_period(isos[hflip]));
    return false;
    
  } 
  
  printf("pass.\n");
  
  printf("\t checking that hflip reverses plane orientation...");

  if (isos[hflip]->crossmap->or != PD_NEG_ORIENTATION) { 

    printf("FAIL. hflip PRESERVES plane orientation.\n");
    return false;

  } 

  printf("pass.\n");

  if (ntwist % 2 == 0) {

    printf("\t checking that hflip reverses curve orientation...");

    if (isos[hflip]->edgemap->or[0] != PD_NEG_ORIENTATION) {

      printf("FAIL. (hflip PRESERVES curve orientation, but ntwist == %d is even).\n",ntwist);
      return false;

    }

    printf("pass.\n");

  } else {

     printf("\t checking that hflip preserves curve orientation...");

    if (isos[hflip]->edgemap->or[0] != PD_POS_ORIENTATION) {

      printf("FAIL. (hflip REVERSES curve orientation, but ntwist == %d is odd).\n",ntwist);
      return false;

    }

    printf("pass.\n");

  }
  
  /* We have now isolated the hflip, vflip, and e operations and tested them. */
  /* The remaining operation is rot. This swaps BOTH pairs of faces. */

  pd_idx_t rot;
  pd_idx_t rotbuf[4];

  printf("\t looking for hflip, which should swap %d edge faces, swap 3 edge faces...",ntwist+1);

  for(i=0,nfound=0;i<nisos;i++) {

    if (isos[i]->facemap->perm->map[0] == 1 &&
	isos[i]->facemap->perm->map[1] == 0 &&
	isos[i]->facemap->perm->map[2] == 3 &&
	isos[i]->facemap->perm->map[3] == 2) {

      rotbuf[nfound++] = i;

    }

  }

  if (nfound != 1) {

    printf("FAIL. Found more than one iso with this property.\n");
    return false;

  }

  rot = rotbuf[0];

  printf("pass. (rot = isos[%d])\n",rot);

  printf("\t checking that rot has period 2...");

  if (pd_iso_period(isos[rot]) != 2) {
    
    printf("FAIL. (isos[rot] == isos[%d] has period %d != 2).\n",
	   rot,pd_iso_period(isos[rot]));
    return false;
    
  } 
  
  printf("pass.\n");
  
  printf("\t checking that rot preserves plane orientation...");

  if (isos[rot]->crossmap->or != PD_POS_ORIENTATION) { 

    printf("FAIL. rot REVERSES plane orientation.\n");
    return false;

  } 

  printf("pass.\n");

  if (ntwist % 2 != 0) {

    printf("\t checking that rot reverses curve orientation...");

    if (isos[rot]->edgemap->or[0] != PD_NEG_ORIENTATION) {

      printf("FAIL. (rot PRESERVES curve orientation, but ntwist == %d is odd).\n",ntwist);
      return false;

    }

    printf("pass.\n");

  } else {

    printf("\t checking that rot preserves curve orientation...");

    if (isos[rot]->edgemap->or[0] != PD_POS_ORIENTATION) {

      printf("FAIL. (rot REVERSES curve orientation, but ntwist == %d is even).\n",ntwist);
      return false;

    }

    printf("pass.\n");

  }
 
  printf("\t free data test... ");

  pd_free_isos(&nisos,&isos);
  pd_free_isos(&nisos,&isos);
  free(pdpA); free(pdpB);

  printf("pass.\n");
 
  printf("-------------------------------------\n");
  printf("%d-twist twist knot tests ... PASS \n\n",ntwist);

  return true;

}  


bool torusknot_iso_test(pd_idx_t p,pd_idx_t q,bool prime,bool print)

{
  pd_code_t *pdpA, *pdpB;
  pd_code_t pdA,pdB;

#define MAXQ 64

  if (p != 2 || q > MAXQ) {

    printf("FAIL. (p,q)-torus knot test not implemented for p != 2, q > %d\n",MAXQ);
    return false;

  }
  
  printf("(%2d,%2d)-torus knot automorphism group test \n",p,q);
  printf("---------------------------------------------\n");
  
  pdpA = pd_build_torus_knot(p,q);
  pdpB = pd_build_torus_knot(p,q);
  
  printf("pd is ");
  pd_write(stdout,pdpA);
  
  pdA = *pdpB; pdB = *pdpB;
  
  printf("\t building (auto)morphisms...");

  unsigned int nisos;
  pd_iso_t **isos = pd_build_isos(pdpA,pdpB,&nisos);

  if (print) {

    unsigned int i;

    printf("\n");

    for(i=0;i<nisos;i++) {

      pd_printf("\t %ISO \n",NULL,isos[i]);

    }

    printf("\t checking # isos generated...");

  }
  
  if (nisos != 4*q) {
    
    printf("FAIL. (# isos generated == %d != %d).\n",nisos,4*q);
    return false;
    
  } else {

    printf("pass (# isos generated == %d == %d).\n",nisos,4*q);
    
  }

  printf("\t finding rots and reflections...");

  pd_idx_t nrots=0,nrefs=0;
  pd_idx_t rotbuf[2*MAXQ],refbuf[2*MAXQ];
  pd_idx_t i;

  for(i=0;i<nisos;i++) {

    if (isos[i]->facemap->perm->map[0] == 0 &&   /* Not a "vflip" */
	isos[i]->facemap->perm->map[1] == 1) {

      if (isos[i]->crossmap->or == PD_POS_ORIENTATION) {
	
	if (isos[i]->edgemap->or[0] != PD_POS_ORIENTATION) {
	  
	  printf("FAIL. isos[%d] preserves plane, but reverses curve orientation.\n",i);
	  return false;
	  
	}
	
	rotbuf[nrots++] = i;
	
      } else {
	
	if (isos[i]->edgemap->or[0] != PD_NEG_ORIENTATION) {
	  
	  printf("FAIL. isos[%d] reverses plane, but preserves curve orientation.\n",i);
	  return false;
	  
	}
	
	refbuf[nrefs++] = i;
	
      }

    }

  }
    
  if (nrots != q) {

    printf("FAIL. Found %d != q == %d rotations.\n",nrots,q);
    return false;

  } 

  if (nrefs != q) {

    printf("FAIL. Found %d != q == %d reflections.\n",nrots,q);
    return false;

  }

  printf("pass (#rots == #refs == q == %d).\n",q);

  printf("\t checking period for refs...");
  
  for(i=0;i<nrefs;i++) { 

    if (pd_iso_period(isos[refbuf[i]]) != 2) {
      
      printf("FAIL. (refbuf[%d] == isos[%d] has period %d != 2).\n",
	     i,refbuf[i],pd_iso_period(isos[refbuf[i]]));
      return false;
      
    } 

  }

  printf("pass. (all reflections have period 2).\n");

  printf("\t looking for e...");

  pd_idx_t e;
  pd_idx_t ebuf[2*MAXQ];
  pd_idx_t efound=0;

  for(i=0;i<nrots;i++) {

    if (pd_iso_is_e(isos[rotbuf[i]])) { 
    
      ebuf[efound++] = i;

    }
    
  }

  if (efound != 1) {

    printf("FAIL. Found %d identity isos.\n",efound);
    return false;

  }

  e = ebuf[0];

  printf("pass. (e == rotbuf[%d] == isos[%d]).\n",e,rotbuf[e]);

  if (prime) {

    printf("\t checking period of (non-e) rotations...");

    for(i=0;i<nrots;i++) {

      if (i != e) {
    
	if (pd_iso_period(isos[rotbuf[i]]) != q) {
	
	printf("FAIL. (rotbuf[%d] == isos[%d] has period %d != q == %d).\n",
	       i,rotbuf[i],pd_iso_period(isos[rotbuf[i]]),q);
	return false;
	
	} 
	
      }
      
    }

    printf("pass (all rots have period q == %d).\n",q);

  }
    
 
  printf("\t free data test... ");

  pd_free_isos(&nisos,&isos);
  pd_free_isos(&nisos,&isos);
  free(pdpA); free(pdpB);

  printf("pass.\n");
 
  printf("-------------------------------------\n");
  printf("(%2d,%2d)-torus knot tests ... PASS \n\n",p,q);

  return true;

}  


bool test_isos() {

  int store_verb = PD_VERBOSE;
  PD_VERBOSE = 0;

  printf("pd_build_isos test suite\n"
	 "-------------------------------\n");

  if (!unknot_iso_test(2,true)) { return false; }
  if (!unknot_iso_test(3,true)) { return false; }
  if (!unknot_iso_test(4,false)) { return false; }
  if (!unknot_iso_test(5,false)) { return false; }
  if (!unknot_iso_test(6,false)) { return false; }
 
  if (!twist_iso_test(3,true)) { return false; }
  if (!twist_iso_test(4,false)) { return false; }
  if (!twist_iso_test(5,false)) { return false; }
  if (!twist_iso_test(6,false)) { return false; }

  if (!torusknot_iso_test(2,3,true,true)) { return false; }
  if (!torusknot_iso_test(2,4,false,false)) { return false; }
  if (!torusknot_iso_test(2,5,true,false)) { return false; }
  if (!torusknot_iso_test(2,6,false,false)) { return false; }
  if (!torusknot_iso_test(2,7,true,false)) { return false; }
  if (!torusknot_iso_test(2,13,true,false)) { return false; }

  printf("------------------------------------\n"
	 "pd_build_isos test suite PASS\n\n");

  PD_VERBOSE = store_verb;

  return true;

}

bool test_isomorphic()
{
  int store_verb = PD_VERBOSE;
  PD_VERBOSE = 0;
  
  printf("pd_isomorphic test suite\n"
	 "-------------------------------\n");

  printf("testing n-crossing unknots, for n in 1..5 ... ");

  pd_idx_t i,j;
  pd_code_t *pdA, *pdB;

  for(i=2;i<=5;i++) { 

    pdA = pd_build_unknot(i);
    
    for(j=2;j<=5;j++) {
      
      pdB = pd_build_unknot(j);

      if (i == j && !pd_isomorphic(pdA,pdB)) { 

	printf("FAIL. %d crossing unknot != %d crossing unknot.\n",i,j);
	return false;

      } else if (i != j && pd_isomorphic(pdA,pdB)) {

	printf("FAIL. %d crossing unknot == %d crossing unknot.\n",i,j);
	return false;

      }

      free(pdB);

    }

    free(pdA);

  }

  printf("PASS. (distinguishes n-crossing unknots 2<=n<=5).\n");   

  printf("testing n-twist twist knots, for n in 1..5 ... ");

  for(i=2;i<=5;i++) { 

    pdA = pd_build_twist_knot(i);
    
    for(j=2;j<=5;j++) {
      
      pdB = pd_build_twist_knot(j);

      if (i == j && !pd_isomorphic(pdA,pdB)) { 

	printf("FAIL. %d-twist twist knot != %d-twist twist knot.\n",i,j);
	return false;

      } else if (i != j && pd_isomorphic(pdA,pdB)) {

	printf("FAIL. %d-twist twist knot == %d-twist twist knot.\n",i,j);
	return false;

      }

      free(pdB);

    }

    free(pdA);

  }

  printf("PASS. (distinguishes n-twist twist knots 2<=n<=5).\n");    

  printf("unknots, torus knots, and twist knots are distinguished...");

  pd_code_t *pdC;
  pd_idx_t k;

  for(i=2;i<=5;i++) { 

    pdA = pd_build_twist_knot(i);
    
    for(j=2;j<=5;j++) {
      
      pdB = pd_build_unknot(j);

      for(k=3;k<=5;k++) {

	pdC = pd_build_torus_knot(2,k);

	if (pd_isomorphic(pdA,pdB)) {
	  
	  printf("FAIL. %d-twist twist knot == %d-crossing unknot.\n",i,j);
	  return false;
	  
	}

	if (pd_isomorphic(pdA,pdC)) {

	  printf("FAIL. %d-twist twist knot == (2,%d)-torus knot.\n",i,k);
	  return false;
	  
	}

	if (pd_isomorphic(pdB,pdC)) {

	  printf("FAIL. %d crossing unknot == (2,%d)-torus knot",j,k);
	  return false;
	  
	}
	
	free(pdC);

      }
	  
      free(pdB);

    }

    free(pdA);

  }

  printf("PASS. (twist knots/unknots/torus knots different).\n"); 
   
  printf("------------------------------------\n"
	 "pd_isomorphic test suite PASS\n\n");
  
  PD_VERBOSE = store_verb;
  
  return true;

}

int main() {

  printf("test_isomorphism (%s)\n",PACKAGE_STRING);
  printf("---------------------------------------\n"
	 "Unit tests for pd_isomorphisms.c\n"
	 "=======================================\n");

  if (!test_isos() || !test_compgrps() || !test_compperms_and_edgemaps() || !test_isomorphic()) {

    printf("=====================================\n");
    printf("test_isomorphism:  FAIL.\n");
    exit(1);

  } else {

    printf("=====================================\n");
    printf("test_isomorphism:  PASS.\n");
    exit(0);

  }

  return 0;

}
  
  

  

  

  
