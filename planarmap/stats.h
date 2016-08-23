typedef struct cumul{
  long *allDist, *maxDist, *gauss, *maxgauss, *geodist;
}pmCumul;

extern void pmStatistic(pmMap *Map, pmStats *Stat, pmCumul *Cumul);

