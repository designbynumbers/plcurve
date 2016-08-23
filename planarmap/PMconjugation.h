extern long pmLuka1 (long n, long k, char *LkWrd );
extern long pmLuka2 (long l, long DgArr[], char LkWrd[] );
extern long pmLuka3 (long i, long j, char LkWrd[] );
extern pm_edge *pmLuka2tree(long st, char LkWrd[]);
extern pm_edge *pmChottin2tree(long st, char LkWrd[]);
extern void pmSpring1(pm_edge *Root);
extern void pmSpring2(pm_edge *Root);
extern void pmSpring3(pm_edge *Root);
extern void pmSpring4(pm_edge *Root);
extern void pmSpring5(pm_edge *Root);
extern pm_edge *pmBalance(pm_edge *Root);
extern pm_edge *pmClosure(pm_edge *Free, pmStck *Stk);
extern pm_edge *pmTwoLegClosure(pm_edge *Free, pmStck *Stk);
extern pm_edge *pmSuppress(pm_edge *Root);
