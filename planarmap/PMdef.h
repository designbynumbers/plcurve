#define TRUE -1
#define FALSE 0

#define MSK 'a'
#define BLACK  2
#define BLACK2 3
#define WHITE  4
#define ROOT   6
#define FREE   8
#define FREE2  9
#define CLSD   10
#define COFREE 12

#define VIRT  14
#define INNER 16
#define OUTER 18
#define DLTD  20

#define FORCED 40
#define DONE   42

#define printvf if (Meth->verbose) printf
#define printwf if (Meth->verbose == TRUE+TRUE) printf

/* les aretes sont toutes pareilles, aussi bien pour les arbres
   que pour les cartes */

typedef struct e /* The data type used for edges */
{ 
  struct v *from;          /* vertex where the edge starts */
  struct v *face;          /* face on the right side of the edge
		      note: only valid if addFace() called */
  struct e *prev;    /* previous edge in clockwise direction */
  struct e *next;    /* previous edge in clockwise direction */
  struct e *oppo;    /* the edge that is inverse to this one */
  long mark;
  short type; 
  long label;
} pm_edge;

#define Next oppo->next   // Pour faciliter les deplacements
#define Prev prev->oppo   // dans la carte duale.
#define Oppo oppo

typedef struct v /* data type for vertices and faces*/
{
  pm_edge *root;
  struct v *next;
  long mark;
  short type;
  long label;
} pm_vertex;
    

typedef struct pmmap /* data type for root and info */
{
  pm_edge *root;
  long e, v, f, i;
  long geodesicDistance;
}pmMap;

/* une carte ou un arbre sont constitues de sommets, d'aretes et de faces
   prises dans un des tableaux correspondant, qui servent de "reservoir":
   la gestion est faite manuellement dans ces tableaux, en particulier
   parce qu'on sait toujours d'avance de quelle place on aura besoin. */


typedef struct st /* data type for stack */
{
  pm_edge **s;
  long pos;
}pmStck;


/******************************/
/* These data types serve for the command flow */
/******************************/

#define PM_MAP_TYPE_QUART_2C 4
#define PM_MAP_TYPE_QUART_2C_2LEG 40

#define PM_BASIC_TYPE_QUART_2C 4

typedef struct pmsize /* Data type for map type and size */
{
  char m, b;     /* map and basic map type */
  long e, v, f;  /* edges, vertices, faces */
  long r, g, d;  /* red and black vertices, max degree */
  long t;        /* tolerence on e */

  // Min number of loop components (inclusive)
  long minLoopComps;
  // Max number of loop components (inclusive)
  long maxLoopComps;

  long *dgArr;   /* pt on vertex list */
} pmSize;

typedef struct pmmethod /* Data type for methods */
{
  char core, pic;
  long seed;
  char verbose;
} pmMethod;

typedef struct pmmemory  /* Data type for memory requirements */
{ 
  char dTree;
  long sTree, rTree, gTree, sWrd, sEdge, sVtx, sLeaf;
} pmMemory;

typedef struct pmoutput /* Data type for output */
{
  char format, transform;
  char map, dual;
} pmOutput;

#define PM_STAT_GEODESICDISTANCE 12

typedef struct pmstats  /* Data type for stats */
{
  long nb;        /* number of maps generated */
  char stats, loop, core, dist, facedeg, gauss, gaussmax, geodist;
} pmStats;


#if 0

#include "PMdef.c"
#include "PMconjugation.c"
#include "PMenlight.c"
#include "PMextract.c"
#include "PMdisplay.c"
#include "PMplanmap.c"

#endif

extern void pmMemoryFault(void);
extern void pmCreateWrd(long n, char **Wrd);
extern void pmFreeWrd(char *Wrd);
extern void pmCreateVtx(long n);
extern void pmFreeVtx(void);
extern pm_vertex *pmNewVtx(pm_edge *Edge);
extern pm_vertex *pmNewFace(pm_edge *Edge);
extern void pmCreateEdge(long n);
extern void pmFreeEdge(void);
extern pm_edge *pmEmptyEdge(void);
extern pm_edge *pmNewEdge(pm_vertex *from, pm_edge *prev, pm_edge *next, pm_edge *oppo, short type);
extern void pmCreateStck(long n, pmStck *Stack);
extern void pmFreeStck(pmStck Stack);
extern void pmStckIn(pm_edge *e, pmStck *Stack);
extern pm_edge *pmStckOut(pmStck *Stack);
extern long pmNewMark(void) ;
extern long pmCurMark(void);
extern long pmNewLabel(void) ;
extern void pmCreateBloc(long n);
extern void pmFreeBloc(void);
extern void pmNewBloc(pm_edge *e);
extern pm_edge *pmNextBloc(void);
extern void pmCreateComp(long n);
extern int pmIsComp(void);
extern void pmFreeComp(void);
extern void pmFirstComp(void);
extern void pmNewComp(pm_edge *e);
extern pm_edge *pmNextComp(void);
extern void pmCreatePost(long n);
extern int pmIsPost(void);
extern void pmResetPost(void);
extern void pmFreePost(void);
extern void pmNewPost(pm_edge *e);
extern pm_edge *pmNextPost(void) ;
extern void pmCopyPostSeed(void);
extern void pmCreateSeed(long n);
extern int pmIsSeed(void);
extern void pmFreeSeed(void);
extern void pmFirstSeed(void);
extern void pmNewSeed(pm_edge *e);
extern pm_edge *pmNextSeed(void);
extern int pmInitRND(pmMethod *Meth);
extern long pmRandom(long n);
extern int pmIsBloc(void);
extern long pmStatGauss(pmMap *Map);
extern long pmStatGeodesicDistance(pmMap *Map);
