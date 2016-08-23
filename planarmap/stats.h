typedef struct cumul{
  long *allDist, *maxDist, *gauss, *maxgauss, *geodist;
}pmCumul;

extern long pmStatGauss(pmMap *Map);
extern void pmStatistic(pmMap *Map, pmStats *Stat, pmCumul *Cumul);

