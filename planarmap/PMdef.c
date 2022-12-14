
/* ici on definit les types et leurs gestions dont on a besoin. */


#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

#include "PMdef.h"

/****************************************/
/* data utility functions */
/****************************************/

/* gestion memoire */
void pmMemoryFault(void){
  printf("\nAllocation failure\n");
  exit(3);
}
/* gestion du mot */
void pmCreateWrd(long n, char **Wrd)
{
  *Wrd=(char *)calloc(n+1,sizeof(char));
  if (*Wrd == NULL) pmMemoryFault();
}
void pmFreeWrd(char *Wrd)
{
  free(Wrd);
}

/* gestion des sommets */
pm_vertex *pmVtxSet; long pmNxtVtxNbr=-1;           

void pmCreateVtx(long n)
{
  pmVtxSet=(pm_vertex *)calloc(n,sizeof(pm_vertex));
  if (pmVtxSet == NULL) pmMemoryFault();
  pmNxtVtxNbr=0;
}

void pmFreeVtx(void)
{
  free(pmVtxSet);
  pmNxtVtxNbr=-1;
}

pm_vertex *pmNewVtx(pm_edge *Edge)
{
  pm_vertex *NxtVtx;
  NxtVtx=pmVtxSet+pmNxtVtxNbr++;
  NxtVtx->root=Edge;
  if (Edge != NULL){
    Edge->from=NxtVtx;
    //    NxtVtx->degree=1;
  }else{
    //    NxtVtx->degree=0;
  }return NxtVtx;
}

pm_vertex *pmNewFace(pm_edge *Edge)
{
  pm_vertex *NxtVtx;
  NxtVtx=pmVtxSet+pmNxtVtxNbr++;
  NxtVtx->root=Edge;
  if (Edge != NULL){
    Edge->face=NxtVtx;
    //    NxtVtx->degree=1;
  }else{
    //    NxtVtx->degree=0;
  }return NxtVtx;
}

/* gestion des aretes */
pm_edge *pmEdgeSet;  long pmNxtEdgeNbr=-1;

void pmCreateEdge(long n)
{
  pmEdgeSet=(pm_edge *)calloc(n,sizeof(pm_edge));
  if (pmEdgeSet == NULL) pmMemoryFault();
  pmNxtEdgeNbr=0;
}

void pmFreeEdge(void)
{
  free(pmEdgeSet);
  pmNxtEdgeNbr=-1;
}

pm_edge *pmEmptyEdge(void)
{
  return (pmEdgeSet+pmNxtEdgeNbr++);
}

pm_edge *pmNewEdge(pm_vertex *from, 
	      pm_edge *prev, pm_edge *next, pm_edge *oppo,
	      short type)
{
  pm_edge *NxtEdge;
  NxtEdge=pmEdgeSet+pmNxtEdgeNbr++;
  NxtEdge->from=from;
  NxtEdge->prev=prev;
  NxtEdge->next=next;
  NxtEdge->oppo=oppo;
  NxtEdge->type=type;
  return NxtEdge;
}

/* gestion de la pile */
void pmCreateStck(long n, pmStck *Stack)
{
  Stack->s=(pm_edge **)calloc(n,sizeof(pm_edge *));
  if (Stack->s == NULL) pmMemoryFault();
  Stack->pos=0;
}

void pmFreeStck(pmStck Stack)
{
  free(Stack.s);
}

void pmStckIn(pm_edge *e, pmStck *Stack)
{
  Stack->s[Stack->pos++]=e;
}

pm_edge *pmStckOut(pmStck *Stack)
{
  if (Stack->pos == 0)
    return NULL;
  else
    return Stack->s[--(Stack->pos)];
}

/* gestion des marques */
long pmAbsMark=1;
long pmAbsLabel=1;

long pmNewMark(void) 
{
  pmAbsMark++;
  if (pmAbsMark == 0) pmAbsMark++;
  return pmAbsMark;
}
long pmCurMark(void)
{
  return pmAbsMark;
}
long pmNewLabel(void) 
{
  pmAbsLabel++;
  if (pmAbsLabel == 0) pmAbsLabel++;
  return pmAbsLabel;
}


/* gestion des sous composantes */
pm_edge **pmBloc;
pm_edge **pmPost;
pm_edge **pmSeed;
pm_edge **pmComp;

long pmBlocBeg, pmBlocEnd, pmCompBeg, pmCompEnd, 
  pmPostBeg, pmPostEnd, pmSeedBeg, pmSeedEnd;
  
void pmCreateBloc(long n)
{
  pmBloc=(pm_edge **)calloc(n,sizeof(pm_edge *));
  if (pmBloc == NULL) pmMemoryFault();
  pmBlocBeg=-1;
  pmBlocEnd=-1;
}
int pmIsBloc(void){
  if (pmBlocBeg == pmBlocEnd) return(FALSE);
  else return(TRUE);
}
void pmFreeBloc(void)
{
  free(pmBloc);
  pmBlocBeg=-1;
  pmBlocEnd=-1;
}
void pmNewBloc(pm_edge *e)
{
  pmBloc[++pmBlocEnd]=e;
}
pm_edge *pmNextBloc(void)
{
  if (pmBlocBeg == pmBlocEnd)
    return NULL;
  else
    return pmBloc[++pmBlocBeg];
}
  
void pmCreateComp(long n)
{
  pmComp=(pm_edge **)calloc(n,sizeof(pm_edge *));
  if (pmComp == NULL) pmMemoryFault();
  pmCompBeg=-1;
  pmCompEnd=-1;
}
int pmIsComp(void){
  if (pmCompBeg == pmCompEnd) return(FALSE);
  else return(TRUE);
}
void pmFreeComp(void)
{
  free(pmComp);
  pmCompBeg=-1;
  pmCompEnd=-1;
}
void pmFirstComp(void)
{
  pmCompBeg = -1;
}
void pmNewComp(pm_edge *e)
{
  pmComp[++pmCompEnd]=e;
}
pm_edge *pmNextComp(void)
{
  if (pmCompBeg == pmCompEnd)
    return NULL;
  else
    return pmComp[++pmCompBeg];
}



void pmCreatePost(long n)
{
  pmPost=(pm_edge **)calloc(n,sizeof(pm_edge *));
  if (pmPost == NULL) pmMemoryFault();
  pmPostBeg=-1;
  pmPostEnd=-1;
}
int pmIsPost(void){
  if (pmPostBeg == pmPostEnd) return(FALSE);
  else return(TRUE);
}
void pmResetPost(void)
{
  pmPostBeg=-1;
  pmPostEnd=-1;
}
void pmFreePost(void)
{
  free(pmPost);
  pmPostBeg=-1;
  pmPostEnd=-1;
}
void pmNewPost(pm_edge *e)
{
  pmPost[++pmPostEnd]=e;
}
pm_edge *pmNextPost(void) 
{
  if (pmPostBeg == pmPostEnd)
    return NULL;
  else
    return pmPost[++pmPostBeg];
}

void pmCopyPostSeed(void)
{
  pmSeedBeg = -1;
  pmSeedEnd = -1;
  while(pmPostBeg < pmPostEnd)
    pmSeed[++pmSeedEnd] = pmPost[++pmPostBeg];
  pmPostBeg = -1;
  pmPostEnd = -1;
}
  
  
void pmCreateSeed(long n)
{
  pmSeed=(pm_edge **)calloc(n,sizeof(pm_edge *));
  if (pmSeed == NULL) pmMemoryFault();
  pmSeedBeg=-1;
  pmSeedEnd=-1;
}
int pmIsSeed(void){
  if (pmSeedBeg == pmSeedEnd) return(FALSE);
  else return(TRUE);
}
void pmFreeSeed(void)
{
  free(pmSeed);
  pmSeedBeg=-1;
  pmSeedEnd=-1;
}
void pmFirstSeed(void)
{
  pmSeedBeg = -1;
}
void pmNewSeed(pm_edge *e)
{
  pmSeed[++pmSeedEnd]=e;
}
pm_edge *pmNextSeed(void)
{
  if (pmSeedBeg == pmSeedEnd)
    return NULL;
  else
    return pmSeed[++pmSeedBeg];
}




/******************************/
/* this function initialize the random generator */
/******************************/
int pmInitRND(pmMethod *Meth){
  srand48(Meth->seed); printvf("# Seed: %ld\n",Meth->seed);
  return(TRUE);
}

/******************************/
/* this function generate random integers between 1,n */
/******************************/
long pmRandom(long n){
  if (n>0) 
    return lrand48()%n+1;
  else return 1; 
}
