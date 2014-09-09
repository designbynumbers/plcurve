/*          FAST VERSION       */
/* A load balanced calculation of the oriented polynomial  Ver 1.1f */
/* Bruce Ewing and Kenneth C. Millett  (author: Ewing)    Feb 5, 01 */
/* Department of Mathematics */
/* University of California */
/* Santa Barbara, CA 93106 */
/* there is a hard limit of 255 crossings in this program */

#include <config.h>

#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#ifdef HAVE_ASSERT_H
#include <assert.h>
#endif

  /* values the USER needs to define!! */
#define XCNT       250     /* maximum crossings in knot (at most 255) */
#define XCNTSQ   62500     /* YOU MUST CALCULATE XCNT SQUARED! */
#define XCSQTR  187500     /* YOU MUST CALCULATE XCNTSQ times THREE! */


/* system values -- it's best if these are not touched */
#define EOFCHR       0     /* integer value of my "end of file" character */
#define MAXBIL   32767     /* positive value where overflow warning triggers */

#define MAXPOLY  12000     /* Maximum # of characters for output poly */

/* globals */

typedef struct state_struct {

  int readpos;
  int writepos;

  long plybuf[XCNTSQ], *poly[XCNT], notbeg, count[XCNT], b[XCNT+1];
  unsigned char sign[XCNT], donlnk[XCNT+1], t[XCNT+2], crsbuf[XCNT+1][8];
  unsigned char buf[XCSQTR], cbuf[10242], clist[XCNT+2], stc[XCNT*2];
  short numcrs, numlps, poslnk, neglnk, lowx, restrt, gapsto[65];
  short tt[XCNT+2], bstlst[XCNT], bilbuf[XCNTSQ], *bilion[XCNT], suplng;

} lmpoly_state_t;

/* Prototypes for later functions */

void mypause(int sig);
int  ntc(long i,char *buf,lmpoly_state_t *state);
int  conchk(lmpoly_state_t *state);
void mrecon(unsigned char *c,unsigned char *p,lmpoly_state_t *state);
void squish(short sqush,lmpoly_state_t *state);
void sqush2(short n,short g, lmpoly_state_t *state);
void triple(short m,short i,unsigned char *p,lmpoly_state_t *state);
void untwst(short l,short n,short i,short ex,lmpoly_state_t *state);
void rmcir(short top,short vnum,short vbrnch,unsigned char *ovundr,int flag,lmpoly_state_t *state);
int  tstcir(short h,short *bstlst,short *dspair,short *skflag,lmpoly_state_t *state);
void delpow(lmpoly_state_t *state);
void addin(long kcoeff,short ypow,short xpow,int mult,int reach,lmpoly_state_t *state);
int  bust(short tobust,short xpow,short ypow,lmpoly_state_t *state);
int  skinny(int twist,short lencir,short lngpos,lmpoly_state_t *state);

int strread(char *str,void *buffer,size_t size,lmpoly_state_t *state)

/* Reads up to size characters from string, copying the results to buffer. */
/* We keep a (global) count of our position in the string internally.      */
/* Return the number of characters read, or 0 at the end of the string.    */
/* We don't copy the terminating NULL from str. */

{
  int ctried,cread;

  for(ctried=0,cread=0;ctried < size;ctried++) {

    if (str[(state->readpos)] != 0) {

      ((char *) buffer)[cread++] = str[(state->readpos)++];

    }

  }

  return cread;
}

int strwrite(char *str,char *buffer,size_t size,lmpoly_state_t *state)

/*   Writes size characters from buffer to str. Returns the number of  */
/*   characters written. We keep a static counter of where we are in   */
/*   the string. */

{
  int ctried,cwrote;

  for(ctried=0,cwrote=0;ctried < size;ctried++) {

    str[(state->writepos)++] = buffer[cwrote++];

    if ((state->writepos) > MAXPOLY-1) {

      printf("lmpoly: Output polynomial larger than %d characters.\n",
	     MAXPOLY);

      exit(1);

    }

  }

  return cwrote;
}


char *plc_lmpoly(char *code, int timeout)
/* This version takes a timeout in seconds. Returns NULL pointer if the computation takes too long. */
{
 register short i, j, k, h, m, n, *sp;
 register unsigned char *p, *c;
 int in = 0, out, stats, len, len2, pause();
 long lngi, *lp1, *lp2, kstrt, cmpval;
 short g, xpow, ypow, dspair[4], maxcrs, chksiz, nopro, skflag;
 unsigned char nbuf[82], *q, *s;
 //struct tms hi;
 char *outpoly;

 clock_t start,end;
 double  cpu_time_used;

 lmpoly_state_t *state = calloc(1,sizeof(lmpoly_state_t));
 assert(state != NULL);

 start = clock();

 if (code == NULL) { free(state); return NULL; }

 /* It's possible that we'll be fed a code with only one crossing */
 /* This will crash the rest of the code, so we insert some stuff here */
 /* to catch this case and return [[0]]N for the unknot. --JHC */

 int newlinecount = 0;
 char *codeptr;
 for(codeptr = code;codeptr != NULL;) {

   codeptr = strchr(codeptr,'\n');
   if (codeptr != NULL) {
     newlinecount++;
     codeptr++;
   }

} /* Count newlines */

 if (newlinecount == 3) { /* Exactly one crossing */

   outpoly = calloc(128,sizeof(char));
   sprintf(outpoly,"[[1]]N ");
   free(state);
   return outpoly;

 }

/* Now back to your regularly scheduled programming. */

 outpoly = calloc(MAXPOLY,sizeof(char)); // Space for a large polynomial
 (state->readpos) = 0;  // Reset the globals
 (state->writepos) = 0;

 maxcrs=XCNT;    /* maximum crossings in a knot, including all limits */
 chksiz= 50;     /* knot size trigger where (state->bilion) is checked, stats printed */
 *(state->sign)= 255;      /* do 3 tests on the machine to verify trick operation */
 i= (short) *(state->sign);
 if (i< (lngi=0)){
  write (1,"unsigned char problem - can't calculate knot >127 crossings\n",60);
  if (maxcrs>127) maxcrs=127;
 }
 kstrt= 0x01020304;
 p= (unsigned char *) &kstrt;
 c= (unsigned char *) &lngi;
 *(c++)= *(p++);
 *(c++)= *(p++);
 *(c++)= *(p++);
 *c= *p;
 if (lngi!= 0x01020304){
  write (1,"long int pointer copy problem -- halt/restart won't work\n",57);
 }
 *nbuf= nbuf[1]= nbuf[2]= nbuf[3]= 0;
 (state->t)[0]=1;
 (state->t)[1]=2;
 (state->t)[2]=3;
 (state->t)[3]=4;
 lp1= (long *) (state->t);
 lp2= (long *) nbuf;
 *lp2= *lp1;
 if (*nbuf!=1 || nbuf[1]!=2 || nbuf[2]!=3 || nbuf[3]!=4){
  write (1,"4 byte copy problem -- use PORTABLE PROGRAM version!\n",53);
  exit (1);
 }
 if (sizeof(lngi)== ((state->suplng)=2)){
  write(1,"halt/restart won't work, and watch for coefficient overflows\n",61);
  cmpval= 10000;
  chksiz= 20;
 }
 else cmpval= 1000000000;     /* how big a value fits in a poly memory slot? */
 if (sizeof(lngi)== 8) (state->suplng)= 1;

 // if (argc<2){
 // write (1,"usage:  lmpoly -rs knotfile [knotfile ...]\n",43);
 // exit (1);
 //}

 /* if ((out=open("lmknot.out",1))== -1){ */
/*   if ((out=creat("lmknot.out",0644))== -1){ */
/*    write (1,"fatal: cannot create output file 'lmknot.out'\n",46); */
/*    exit (1); */
/*   } */
/*  } */
/*  else lseek (out,(long) 0,2); */
/*  close (out); */

 *(state->count)= (state->restrt)= i= j= kstrt= nopro= 0;
 stats= -1;

 /* This code appears to detect "restarts" and "stats", neither of which
    we are currently implementing. */

 /* if (*argv[1]=='-'){
  if (argc==2 && argv[1][1]=='r'){
   (state->restrt)=1;
   if (argv[1][2]=='s') stats= 0;
  }
  else if (argv[1][1]=='s'){
   stats= 0;
   ++argv;
   --argc;
  }
  } */


 signal(SIGTERM,mypause); // Not sure what this does.

 lp1= state->plybuf;
 sp= (state->bilbuf);
 while (j<XCNT){
  (state->poly)[j]= lp1;
  (state->bilion)[j]= sp;
  lp1+= 1+ 2* (XCNT- ++j);
  sp+= 1+ 2* (XCNT-j);
 }
NEWFIL:

 /* We comment out some "stats" code, since we won't display stats anyway. */

 /* if (stats>0){     // remove stats info for previous knot file
  i= strlen(*argv);
  c= (unsigned char *) *argv+i-1;
  if (i>11) i=11;
  while (i--!=0 && *c!='/') --c;
  p=t;
  while (*++c!=0) *(p++)= *c;
  *(p++)= '.';
  *(p++)= 'l';
  *p=0;
  close (stats);
  unlink(t);
  } */
 /* if (--argc==0) exit(0); */ /* This was the exit point for the program */

 /* We also comment out some "restart" code since we aren't restarting */

 /*if ((state->restrt)!=0){           // if I am restarting an old calculation
  if ((in=open("lmpoly.(state->restrt)",0))== -1){
   write (1,"couldn't find restart file\n",27);
   exit (0);
  }
  read (in,nbuf,82);
  read (in,(state->cbuf),12);
  *argv= (char *) (state->cbuf);
  read (in,(state->buf),7);
  (state->numcrs)= *(state->buf);
  xpow= (short) (state->buf)[2];
  ypow= (short) (state->buf)[1];
  (state->numlps)= (short) (state->buf)[3];
  (state->poslnk)= (short) (state->buf)[4];
  (state->neglnk)= (short) (state->buf)[5];
  if (((state->buf)[6]&1) !=0) ypow= -ypow;
  if (((state->buf)[6]&2) !=0) xpow|= 512;
  read (in,*(state->crsbuf),(state->numcrs)*8);
  i= 0;
  while (i!=(state->numcrs)){
   (state->sign)[i]= (state->crsbuf)[i][1]>>4;
   (state->crsbuf)[i++][1]&=6;
  }
  p= (unsigned char *) &(state->lowx);
  read (in,p,2);
  p= (unsigned char *) &kstrt;
  read (in,p,4);
  p= (unsigned char *) &notbeg;
  read (in,p,4);
  read (in,(state->buf), (int) notbeg);
  read (in,(char *) count,XCNT*4);
  read (in,(char *) state->plybuf,4*XCNT*XCNT);
  read (in,(char *) (state->bilbuf),2*XCNT*XCNT);
  close (in);
  in=0;
  s= (state->cbuf)+14;
  *s= EOFCHR;
 }
 else if ((in=open(*++argv,0))==-1){
  write (1,"cannot open ",12);
  write (1,*argv,strlen(*argv));
  write (1,", skipping\n",11);
  goto NEWFIL;
 }
 if (stats>=0){
  i= strlen(*argv);
  q= (unsigned char *) *argv+i-1;
  if (i>11) i=11;
  while (i--!=0 && *p!='/') --q;
  p=t;
  while (*++q!=0) *(p++)= *q;
  *(p++)= '.';
  *(p++)= 'l';
  *p=0;
  if ((stats=creat(t,0644)) <0) stats=0;
 }
 if ((state->restrt)!=0){
  if (stats>0){
   write (stats,nbuf,strlen(nbuf));
   write (stats,"\n",1);
  }
  (state->restrt)=0;
  goto STEP1; // restart file has set up all NEWNOT data already
 }
 */

 c = (state->cbuf);

 /* This code fragment copies the input file into a character buffer,
    appending the EOF character to the end of the string. Instead, we
    copy our input "code" buffer into the (state->cbuf) array.

    Original code:

    c[read(in,c,10240)]= EOFCHR;

 */

 int   mypos;
 for(mypos = 0;code[mypos] != 0 && mypos < 10240-1;mypos++) {
   c[mypos] = code[mypos];
 }
 c[mypos]= EOFCHR;

 while (*c=='\n') ++c;

 s=c;
NEWNOT:
 i= XCNT-1;
 while (i!=0 && (state->count)[i]==0) --i;
 if (*(state->count)!= ((state->numcrs)=0)){
  if (stats>=0){
   if ((out=open("lmknot.stats",1))== -1) out=creat("lmknot.stats",0644);
   else lseek (out,(long) 0,2);
  }
/*   if (stats>=0 && out>0){ */
/*    write (out,nbuf,strlen(nbuf)); */
/*    write (out,"\n",1);     /\* write out statistics for knot just completed *\/ */
/*    lngi=0; */
/*    while (i!=0){ */
/*     write (out,t,ntc((long)(i+1),t,state)); */
/*     write (out,"     ",5); */
/*     write (out,t,ntc(count[i],t)+1,state); */
/*     lngi+= count[i--]; */
/*    } */
/*    write (out,"\ntotal: ",8); */
/*    write (out,t,ntc(lngi,t,state)+1); */
/*    write (out,"run: ",5); */
/*    times(&hi); */
/*    lngi= ((hi.tms_utime-kstrt)*5)/3; */
/*    write (out,t,ntc(lngi/100,t,state)); */
/*    write (out,".",1); */
/*    len= ntc(lngi%100,t,state); */
/*    write (out,"0",2-len); */
/*    write (out,t,len); */
/*    write (out," s\n\n",4); */
/*    kstrt= hi.tms_utime; */
/*    if (count[2]<0) write (out,"program error: knot became inconsistent\n",40); */
/*    close (out); */
/*   } */
  k= XCNT;
  i=0;
  while (k!=0 && i==0){
   i= 2+ 2*(XCNT-k);
   lp1= (state->poly)[--k];
   sp= (state->bilion)[k];
   while (--i!=0 && *lp1==0 && *(sp++)==0) ++lp1;
  }
  if (i!=0){

    /* Original code:

       if ((out= open("lmknot.out",1))<0) out=open("temp.out",1);

       lseek (out,(long) 0,2);
       write (out,nbuf,strlen(nbuf));
       write (out,"\n",1); */

    /* write out polynomial for knot just completed */

    /* These lines would copy the first line of the crossing code into
       the output. We don't need this, so we forget it */

    /*strwrite(outpoly,(char *)(nbuf),strlen((char *)(nbuf)));
      strwrite(outpoly," ",1); */

    /*
    if (count[2]<0) write (out,"program error: knot became inconsistent\n",40);
    if (count[1]!=0) write (out,"coefficient overflow error: output BAD\n",39);
    */

    if ((state->count)[2]<0) {

      printf("lmpoly: knot became inconsistent.\n");
      exit(1);

    }

    if ((state->count)[1]!=0) {

      printf("lmpoly: coefficient overflow error.\n");
      exit(1);

    }

    len=m= -1;
    while (m++<k || (state->lowx) <= 0){
        if ((state->lowx)== (i=0)) strwrite (outpoly,"[",1,state);
        n= XCNT-m-1;
        j= n*2;
        while (j!=n && (state->poly)[m][j]==0 && (state->bilion)[m][j]==0) --j;
        while (i!=n && (state->poly)[m][i]==0 && (state->bilion)[m][i]==0) ++i;
        if (len==0 || (state->lowx)>=0 || i!=j || (state->poly)[m][i]!=0 || (state->bilion)[m][i]!=0){
            while (i<=j){
                if (i==n) strwrite (outpoly,"[",1,state);
                h=0;
                lngi= (state->poly)[m][i];
                if (lngi>=cmpval || lngi<= -cmpval){
                    h= lngi/cmpval;
                    lngi-= h* cmpval;
                }
                h+= (state->bilion)[m][i];
                if (h*lngi <0){   /* (state->bilion) and poly are different signs */
                    if (h<0){
                        lngi-= cmpval;
                        ++h;
                    }
                    else {
                        lngi+= cmpval;
                        --h;
                    }
                }
                if (h!=0){
                    if (lngi<0) lngi= -lngi;
                    strwrite (outpoly,(char *)(state->t),ntc((long) h,(char *)(state->t),state),state);
                    len= ntc(lngi,(char *)(state->t),state);
                    if (cmpval==10000) len2= 4-len;
                    else len2= 9-len;
                    strwrite (outpoly,"00000000",len2,state);
                    strwrite (outpoly,(char *)(state->t),len,state);
                }
                else strwrite (outpoly,(char *)(state->t),ntc(lngi,(char *)(state->t),state),state);
                if (i++ ==n) strwrite (outpoly,"]",1,state);
                if (i<=j) strwrite (outpoly," ",1,state);
            }
            if ((state->lowx)== (len=0)) strwrite (outpoly,"]",1,state);
            strwrite (outpoly,"N",1,state);
        }
        ++(state->lowx);
    }
  }
  strwrite (outpoly," ",1,state);
  /* close(out); */ /* We are done writing the polynomial, so return */
  free(state);
  return outpoly;
 }
 i= XCNT;
 while (i!=0) (state->count)[--i]=0; /* Erase the (state->count) buffer */
 c=s;
 *nbuf=0;                   /* Erase the nbuf string? */
 if (*c!='1' || (*(c+1)!='+' && *(c+1)!='-')){
  /* Copy first line of 's' into nbuf */
  p= nbuf;
  while (*c!='\n' && *c!= EOFCHR) *(p++)= *(c++);
  *p= 0;
  /* if (*c != EOFCHR && stats>0){ */
/*    lseek (stats,(long) 0,0); */
/*    write (stats,nbuf,strlen(nbuf)); */
/*    write (stats,"\n",1); */
/*   } */
 }
 if (*c=='\n') ++c;
 p= *(state->crsbuf);
 while (*c!='\n' && *c!=EOFCHR){
  if ((state->numcrs)==maxcrs){
    //write (1,"too many crossings in knot\n",27);
    printf("lmpoly: too many crossings in knot\n");
    free(outpoly);
    free(state);
    return NULL;
  }
  while (*c>='0' && *c<='9') ++c;
  if (*c=='+') (state->sign)[(state->numcrs)]=6;      /* (state->sign)[] says what crossings are + or - */
  else if (*c=='-') (state->sign)[(state->numcrs)]=2;
  else if (*c==EOFCHR) --c;
  else {
    //write (1,"the format of this knot is unreadable, skipping file\n",53);
    //close (in);
    //goto NEWFIL;
    printf("lmpoly: Could not parse input string.\n");
    free(outpoly);
    free(state);
    return NULL;
  }
  ++c;
  j=4;
  while (j--!= (i=0) && *c!=EOFCHR){
   while (*c>='0' && *c<='9'){
    i*=10;
    i+= *(c++)-0x30;
   }
   if (i==0){
     //write (1,"the format of this knot is unreadable, skipping file\n",53);
     //close (in);
     //goto NEWFIL;
     printf("lmpoly: Could not parse input string.\n");
     free(outpoly);
     free(state);
     return NULL;
   }
   *(p++)= i-1;
   if (*c!=EOFCHR){
    if (*c<'a' || *c>'d'){
      //write (1,"the format of this knot is unreadable, skipping file\n",53);
      //close (in);
      //goto NEWFIL;
      printf("lmpoly: Could not parse input string.\n");
      free(outpoly);
      free(state);
      return NULL;
    }
    *(p++)= (*(c++)-'a')*2;
   }
  }
  ++(state->numcrs);
  if (*c=='\n') ++c;
 }
 if (*c==EOFCHR){
  c=(state->cbuf);
  len=10240;
  while (*s!=EOFCHR){
   *(c++)= *(s++);
   --len;
  }
  if (in!=0 && (i=read(in,c,len))!=0){
   c[i]=EOFCHR;
   c=(state->cbuf);
   while (*c=='\n') ++c;
   s=c;
   goto NEWNOT;
  }
  *c= EOFCHR;
  if (in!=0){
   close (in);
   in=0;
  }
  else goto NEWFIL;
  if ((state->numcrs)<2) goto NEWFIL;
 }
 while (*c=='\n') ++c;
 s=c;
 lp1= *(state->poly);
 sp= *(state->bilion);
 lngi= XCNT*XCNT;
 while (--lngi!=0) *(lp1++)= *(sp++)= 0;
 (state->notbeg)=(state->numlps)=(state->poslnk)=(state->neglnk)= xpow= ypow= 0;
 if (conchk(state)!=0) exit(0);
 *(state->count)=1;
STEP1:
 skflag=0;
RESTRT:
 if ((state->restrt)!=0){
   /* We don't include the restart functionality for now. */

  /* if ((out=creat("lmpoly.restrt",0644))!= -1){ */
/*    write (out,nbuf,82); */
/*    i= strlen(*argv)-11; */
/*    if (i<0) i=0; */
/*    write (out,*argv+i,12); */
/*    *(state->cbuf)= (state->numcrs); */
/*    if (ypow<0){ */
/*     ypow= -ypow; */
/*     (state->cbuf)[6]=1; */
/*    } */
/*    else (state->cbuf)[6]=0; */
/*    if ((xpow&512) !=0) (state->cbuf)[6]|= 2; */
/*    (state->cbuf)[1]= ypow; */
/*    (state->cbuf)[2]= xpow; */
/*    (state->cbuf)[3]= (state->numlps); */
/*    (state->cbuf)[4]= (state->poslnk); */
/*    (state->cbuf)[5]= (state->neglnk); */
/*    write (out,(state->cbuf),7); */
/*    i= -1; */
/*    while (++i!=(state->numcrs)) (state->crsbuf)[i][1]|= (state->sign)[i]<<4; */
/*    write (out,*(state->crsbuf),(state->numcrs)*8); */
/*    p= (un(state->sign)ed char *) &(state->lowx); */
/*    write (out,p,2); */
/*    times(&hi); */
/*    kstrt-= hi.tms_utime; */
/*    p= (un(state->sign)ed char *) &kstrt; */
/*    write (out,p,4); */
/*    p= (unsigned char *) &notbeg; */
/*    write (out,p,4); */
/*    write (out,buf, (int) notbeg); */
/*    write (out,(char *) count,XCNT*4); */
/*    write (out,(char *) state->plybufn,4*XCNT*XCNT); */
/*    write (out,(char *) (state->bilbuf),2*XCNT*XCNT); */
/*    close (out); */
/*    close (stats); */
/*    exit(0); */
/*   } */
   /*write (1,"couldn't create restart file 'lmpoly.restrt'\n",45);*/

   printf("lmpoly: Warning! This version does not save"
	  " restart info on quit.\n"
	  "        Terminating run now.\n\n");
   exit(0);
   /*restrt=0; */
 }
/* FIRST -- remove figure 8 loops, 2 cases  "#-#d#c#b#a"  "#+#b#a#d#c" */
/* SECOND -- remove monogons (must recheck previous crossings) */
/* 4 cases  "#+#b#a????"  "#+????#d#c"  "#-#d????#a"  "#-??#c#b??" */
/* THIRD -- remove bigons -- make sure loops don't get forgotten. */
/* 2 cases "k?mcmb???? m???kbka??"  "k?mc????md  m?????kakd" */
/* bigon yanker MUST be completely thorough or loops can get lost later */
 *dspair=g= i= (state->numcrs);
 c=p= (state->crsbuf)[i];
 k=2;
 while (i--!=0){  /* loop through all crossings */
  p-=8;
  if (*p==i){
   if (p[4]==i) ++(state->numlps);
   else mrecon (p+4,p+(state->sign)[i],state);
   squish (i,state);
   goto STEP1;
  }
  else if (p[4]==i){
    mrecon (p,p+((state->sign)[i]^4),state);
    squish (i,state);
   goto STEP1;
  }
   /* this bigon test checks that ka points to some mc, */
    /* and kd to md, or kb to mb */
  if (p[1]==4 && ((*p==p[6] && p[7]==6) || (*p==p[2] && p[3]==2))){
   if ((j= *p) ==p[2]) k=6;
   c= (state->crsbuf)[j];
   if (p[4]==j) ++(state->numlps);
   if (p[2]==p[6]) ++(state->numlps);
   mrecon (p+4,c,state);
   mrecon (p+k,c+k,state);
   sqush2 (i,j,state);
   goto STEP1;
  }
  (state->bstlst)[i]= 0;
 }
 if ((state->numcrs)<6) g=0;
/* find very good and ok triples */
 while (g--!= (lngi=0)){
  c-=8;
  p= (state->crsbuf)[*c]+ ((c[1]^2)&2);
  q= (state->crsbuf)[c[4]]+ ((c[5]^2)&2);
  if ((i= c[2]) != (j= c[6]) && *c!=c[4]){
   if ((m= *q) == (n= q[4])) m=n= -1;
   if ((h= *p) == (k= p[4])) h=k= -1;
   if ((h==i && k==j) || (j==h && k==i)) lngi=1;
   if ((m==i && n==j) || (m==j && n==i)) lngi|=2;
   if ((n==j || m==j) && (k==j || h==j)) lngi|=4;
   if ((m==i || n==i) && (h==i || k==i)) lngi|=8;
  }
  else {
   h=0;       /* there is a ring through this crossing */
   if (i!=j){
    i= *c;   /* on over branch */
    h=2;
   }
   else if (*c==c[4]){
    ++(state->numlps); /* knot's a distant union with a 2-link */
    if ((state->sign)[g]==2) ++(state->neglnk);
    else ++(state->poslnk);
    sqush2 (g,i,state);
    goto STEP1;
   }
   j=h;
   c+=h;
   if ((state->sign)[g]!=(state->sign)[i]) ++(state->numlps); /* loose ring */
   else {
    if ((state->sign)[g]==2) ++(state->neglnk);
    else ++(state->poslnk);
    j= h^2;
   }
   p= (state->crsbuf)[i]+j;  /* only do one mrecon unless I HAVE to do 2 */
   q= c+4;
   if (*c==i){
    if ((c[1]&4) ==0) c= p+4;
    else c=p;
   }
   else if (*q==i){
    if ((q[1]&4) ==0) q= p+4;
    else q=p;
   }
   else mrecon (p,p+4,state);
   mrecon (c,q,state);
   sqush2 (g,i,state);
   goto STEP1;
  }
  if (lngi!=0){
   k= lngi;
   *dspair=1;
   if ((k&3)!=0){
    if ((c[3]&c[7]&2)!=0){
      triple((short) ((k&2)<<1),g,c,state);
     goto STEP1;
    }
    if (((c[3]^c[7])&2)!=0){
     if ((c[((k&2)<<1)|1]&2)==0){
      if ((state->sign)[c[(k&2)<<1]]!=(state->sign)[g]) (state->bstlst)[i]= (state->bstlst)[j]= -19;
      else (state->bstlst)[i]= (state->bstlst)[j]= 8;
     }
     else {
      n=8;
      if ((state->sign)[c[(k&2)<<1]]==(state->sign)[g]) n= -8;
      sp= (state->bstlst)+i;
      if ((c[3]&2)!=0) sp= (state->bstlst)+j;
      if (*sp==0 || (n<0 && n<*sp)) *sp= n;
     }
    }
    else {
     if ((state->bstlst)[g]>=0) (state->bstlst)[g]+= 8;
     else (state->bstlst)[g]-= 8;
    }
   }
   if (k>3){
    if (((c[1]|c[5])&2)==0){
      triple((short) ((k&4)|2),g,c,state);
     goto STEP1;
    }
    if (((c[1]^c[5])&2)!=0){
     if ((c[(k&4)|3]&2)!=0){
      if ((state->sign)[c[(k&4)|2]]!=(state->sign)[g]) (state->bstlst)[*c]= (state->bstlst)[c[4]]= -19;
      else (state->bstlst)[*c]= (state->bstlst)[c[4]]= 8;
     }
     else {
      n=8;
      if ((state->sign)[c[(k&4)|2]]==(state->sign)[g]) n= -8;
      sp= (state->bstlst) + c[4];
      if ((c[1]&2)!=0) sp= (state->bstlst)+ *c;
      if (*sp==0 || (n<0 && n<*sp)) *sp= n;
     }
    }
    else {
     if ((state->bstlst)[g]>=0) (state->bstlst)[g]+= 8;
     else (state->bstlst)[g]-= 8;
    }
   }
  }
 }
 /* circuit remover */
 j= (state->numcrs)-1;
 if (j<3) j= -1;
 else if (j<11){    /* <12 will blow up */
  p= *(state->crsbuf)+1;
  i=(state->numcrs);
  while (i--!=0 && *p!=4) p+=8; /* only non-alternating (a to c) knots */
  if (i< (n=0)) j= -1;
  else p= (state->crsbuf)[j]+5;   /* remove twists */
  while (j>=0){
   n^=4;
   k= *(p--)^4;
   if ((k&2) == 0){
    c= (state->crsbuf)[*p]+k;
    i= *c;
    k= c[1]^4;
    if ((k&2) != 0 && (state->crsbuf)[i][k]==j){
      untwst ((short) *p,i,j,n,state);
      squish (j,state);
     goto STEP1;
    }
   }
   if (n==0) --j;
   p-=3;
  }
 }
 g= (state->donlnk)[XCNT]= m= 0;
 while ((n=j) >= (*(state->tt)=(state->tt)[2]=h= 0)){
  k=g;
  do {
   p= (state->crsbuf)[j]+k;
   j= *p;
   k= p[1]^4;
   i= h;
   p=(state->clist);
   while (i--!=0 && *(p++) != j);
   if (i>=0){
    j= n;              /* circuit crosses itself */
    goto NXTCIR; /* skip circuit, wait for a sub-piece that doesn't */
   }
   else {
    ++(state->tt)[k&2];
    *p= j;
    (state->t)[h++]= (k&2)^2;
   }
  } while (n!=j);
  if (k != g) --(state->tt)[k&2];
  else k= -1;
  if (*(state->tt)==0 || (state->tt)[2]==0){
    rmcir(h,j,k,state->t,0,state);
   goto STEP1;
  }
  if (tstcir(h,(state->bstlst),dspair,&skflag,state) !=0) goto RESTRT;
  if (h==3){
   (state->stc)[m++]=j; /* store vertex of 3 crossing circuit for (state->bstlst) later */
   (state->stc)[m++]=g;
   (state->stc)[m++]=k^4;
  }
NXTCIR:
  if (g==0) g= (state->sign)[j];
  else {
   --j;
   g=0;
  }
 }
/* if knot is tiny, store polynomial coefficients and restart a stored knot */
 if ((state->numcrs)<6){
  lngi= 1;
  if ((state->numcrs)==0) --(state->numlps);
  delpow(state);
  if ((xpow&512)!=0){
   lngi= -1;
   xpow&=511;
  }
  if ((state->numcrs)==4){
   c= *(state->crsbuf);
   p= c+4;
   if (*p==*c || c[2]==p[2]){
    g= -3;
    if (*p==p[2] && *c==c[2]) g= -4;
   }
   else {
    g= 0;  /* 0 is flag for knot */
    i= *(state->sign)+(state->sign)[1]+(state->sign)[2]+(state->sign)[3];
    if (i!=16) g=2;     /* flag for 2 links   */
   }
  }
  else if ((state->numcrs)==5){
   n=10;
   i=j= 0; /* step through knot, if don't get to end, it's a link */
   do {
    p= (state->crsbuf)[i]+j;
    i= *p;
    j= p[1]^4;
    --n;
   } while (i!=0 || j!=0);
   if (n>0){
    i=h= 0; /* 5 crossing 2 link  or twisted 3 link  or tref and 2 links */
    g=m= 1;      /* or tref + link */
    k= *(state->sign)+(state->sign)[1]+(state->sign)[2]+(state->sign)[3]+(state->sign)[4];
    if (k>10 && k<30){
     if (k>20){
      while ((state->sign)[i]==6) ++i;
     }
     else {
      while ((state->sign)[i]==2) ++i;
     }
     if (k==14 || k==26) goto FOURX;
     c= (state->crsbuf)[i];
     p= c+4;
     if (*p==*c || p[2]==c[2]){
      if (*p==p[2] && *c==c[2]) g=3;
      else {
       h=4;
       if ((state->sign)[i]==6) m= -1;
      }
     }
     else {
      g=2;
      j=0;
      if ((state->sign)[*c]!=(state->sign)[i]) j=4;
      k= c[j];
      if ((state->crsbuf)[k][j]!=i){
       if ((state->sign)[i]==6) m= -1;
      }
      else if ((state->sign)[i]==2) m= -1;
     }
    }
    if (g==1){
     if (h==0 && (state->sign)[i]==2) m= -1;
     i=j=k= 5;
     c= state->t;
     while (j!=0) c[--j]=0;
     while (k--!=0){
      j= (state->sign)[k]; /* trefs + circles have fake bigons at all crossings */
      p= (state->crsbuf)[k]; /* so find all fake bigons */
      if (p[1]==(j^4) && (*p==p[2] || *p==p[6])) c[k]=c[*p]= 1;
      if (p[5]==j && p[4]==p[j]) c[k]=c[p[4]]= 1;
     }
     while (--i>=0 && c[i]!=0);
     if (i>=0){
FOURX: /* found twisted 3 link -- "squish" twist before making polynomial */
      while (++i!=6) (state->sign)[i-1]= (state->sign)[i];
      (state->numcrs)=4;
      g= -3;  /* flag for 4 crossing 3 links */
     }
     else {
      p= *(state->crsbuf);
      if (n==4 || (n==8 && *p==p[4] && p[2]==p[6])) g=3;
     }
    }
   }
   else g= -1; /* it's a knot! */
  }
  i= (*(state->sign))/2 -2;
  k= xpow;
  n=ypow;
  p= *(state->crsbuf);
  if ((state->numcrs)==5){
   if (g<0){ /* knots */
    if (i>0){
     if (*p==p[6] && p[2]==p[4] && p[10]==p[12]) g= -2;
    }
    else if (*p==p[2] && p[4]==p[6] && p[8]==p[10]) g= -2;
    if (g== -1){ /* knot five-two */
      addin(lngi,(short)(n-i*6),k,0,0,state);      /* the flags at the end of addin  */
      addin(lngi,(short)(n-i*4),k,-1,2,state);     /* let you do two addins at once, */
      addin((long)-lngi,(short)(n-i*2),k,-1,2,state);    /* for addins with equal ypows.   */
    }
    else { /* knot five-one */
     k+=2;
     addin((long)(-4*lngi),(short)(n-i*4),k,0,0,state);
     addin((long)-lngi,(short)(n-i*6),k,-2,-2,state);
     addin(lngi,(short)(n-i*4),(short)(k+2),3,-4,state);
    }
   }
   else { /* links */
    if (g==1){ /* trefoil + circle */
     if (h!=0){
       addin(lngi,(short)(n-m*3),++k,-1,-2,state);
       addin((long)(3*lngi),(short)(n-m),k,-1,-2,state);
       addin(lngi,(short)(n+m),k,-2,-2,state);
       addin((long)-lngi,(short)(n-m),(short)(k+2),0,0,state);
     }
     else {
       addin((long)-lngi,(short)(n-m*7),--k,0,0,state);
       addin((long)(-3*lngi),(short)(n-m*5),k,0,0,state);
      k+=2;
      addin((long)(3*lngi),(short)(n-m*3),k,0,0,state);
      addin((long)(2*lngi),(short)(n-m*5),k,0,0,state);
      addin((long)-lngi,(short)(n-m*3),(short)(k+2),2,-4,state);
     }
    }
    else if (g==2){ /* 2 component link five-one */
      addin((long)-lngi,(short)(n+m),--k,-2,2,state);
      addin((long)-lngi,(short)(n-m),k,-1,2,state);
      addin(lngi,(short)(n+m*3),(short)(k+2),0,0,state);
      addin((long)-lngi,(short)(n+m),(short)(k+4),0,0,state);
    }
    else { /* trefoil and circle + circle */
     i= *(state->sign)+(state->sign)[1]+(state->sign)[2]+(state->sign)[3]+(state->sign)[4];
     m=1;
     if (i<20){
      m= -1;
      i= 40-i;
     }
     j= k-2;
     if (i!=30){
       addin(lngi,(short)(n-m*4),j,-1,2,state);
       addin((long)(4*lngi),(short)(n-m*2),j,-1,2,state);
       addin((long)(5*lngi),n,j,0,0,state);
       addin((long)-lngi,(short)(n+m*2),k,-2,-2,state);
      k+=2;
      addin(lngi,n,k,-4,-2,state);
      addin(lngi,(short)(n-m*2),k,0,0,state);
     }
     else {
       addin(lngi,(short)(n-m*8),j,0,0,state);
       addin((long)(5*lngi),(short)(n-m*4),j,0,0,state);
       addin((long)(2*lngi),(short)(n-m*2),j,0,0,state);
       addin((long)(-2*lngi),(short)(n-m*6),k,-2,-2,state);
      k+=2;
      addin(lngi,(short)(n-m*4),k,-5,-2,state);
      addin(lngi,(short)(n-m*2),k,-3,-2,state);
     }
    }
   }
  }
  else if ((state->numcrs)==4){
   if (g==0){ /* knot 4-1 */
     addin((long)-lngi,(short)(n-2),k,0,0,state);
     addin((long)-lngi,n,k,-1,2,state);
     addin((long)-lngi,(short)(n+2),k,0,0,state);
   }
   else if (g>0){
    i=1;
    if (*(state->crsbuf)[*p]==0){
     if (*(state->sign)==6) i= -1; /* 2 link 4^2-1 */
     addin((long)-lngi,(short)(n+i*3),--k,-1,2,state);
     addin((long)-lngi,(short)(n+i*5),k,0,0,state);
     addin((long)-lngi,(short)(n+i),(short)(k+2),0,0,state);
    }
    else {
     if (*(state->sign)==2) i= -1;
     addin((long)-lngi,(short)(n-i*5),--k,-1,2,state);
     addin((long)-lngi,(short)(n-i*3),k,-3,2,state);
     addin((long)-lngi,(short)(n-i*3),(short)(k+4),0,0,state);
    }
   }
   else if (g== -3){ /* 3 link */
    j= *(state->sign)+(state->sign)[1]+(state->sign)[2]+(state->sign)[3];
    i= 8- j/2;
    if (i==0){
      addin((long)-lngi,(short)(n-2),k,-1,-2,state);
      addin((long)(-2*lngi),n,k,-1,-2,state);
      addin((long)-lngi,(short)(n+2),k,-1,-2,state);
      addin(lngi,n,(short)(k+2),0,0,state);
    }
    else {
     j= k-2;
     addin(lngi,(short)(n-2+i),j,0,0,state);
     addin((long)(2*lngi),(short)(n+i),j,-1,2,state);
     addin(lngi,(short)(n+2+i),j,0,0,state);
     addin(lngi,(short)(n+i/2),(short)(k+2),-2,-2,state);
    }
   }
   else { /* circle + circle  and  circle + circle */
    j= *(state->sign)+(state->sign)[1]+(state->sign)[2]+(state->sign)[3];
    i= 8- j/2;
    if (i==0){
      addin(lngi,(short)(n-3),--k,-1,-2,state);
      addin((long)(3*lngi),(short)(n-1),k,-1,-2,state);
      addin((long)(3*lngi),(short)(n+1),k,-1,-2,state);
      addin(lngi,(short)(n+3),k,-1,-2,state);
     k+=2;
     addin((long)-lngi,(short)(n-1),k,0,0,state);
     addin((long)-lngi,(short)(n+1),k,0,0,state);
    }
    else {
     k-=3;
     addin((long)-lngi,(short)(n-3+i),k,0,0,state);
     addin((long)(-3*lngi),(short)(n-1+i),k,0,0,state);
     addin((long)(-3*lngi),(short)(n+1+i),k,0,0,state);
     addin((long)-lngi,(short)(n+3+i),k,0,0,state);
     k+=2;
     addin((long)(2*lngi),(short)(n-2+(i*3)/4),k,0,0,state);
     addin((long)(4*lngi),(short)(n+(i*3)/4),k,0,0,state);
     addin((long)(2*lngi),(short)(n+2+(i*3)/4),k,0,0,state);
     k+=2;
     addin((long)-lngi,(short)(n-1+i/2),k,0,0,state);
     addin((long)-lngi,(short)(n+1+i/2),k,0,0,state);
    }
   }
  }
  else if ((state->numcrs)==3){
    addin((long)-lngi,(short)(n-i*4),k,0,0,state);
    addin(lngi,(short)(n-i*2),(short)(k+2),-2,-2,state);
  }
  else if ((state->numcrs)!=0){
    addin(lngi,(short)(n-i),--k,-1,2,state);
    addin(lngi,(short)(n-i*3),k,0,0,state);
  }
  else addin(lngi,n,k,0,0,state);
  if ((state->notbeg)!=0){
   (state->numcrs)= (state->buf)[(state->notbeg)-1];
   (state->numlps)= (state->buf)[(state->notbeg)-4];
   (state->poslnk)= (state->buf)[(state->notbeg)-5];
   (state->neglnk)= (state->buf)[(state->notbeg)-6];
   xpow= (state->buf)[(state->notbeg)-3];
   ypow= (state->buf)[(state->notbeg)-2];
   if ((state->buf)[(state->notbeg)-7]!=0) ypow= -ypow;
   if ((state->buf)[(state->notbeg)-8]!=0) xpow|= 512;
   (state->notbeg)-= 8+ (state->numcrs)*8;
   i=0;
   p= *(state->crsbuf);
   c= (state->buf)+(state->notbeg);
   while (i!=(state->numcrs)){
    j=4;
    while (j--!=0){
     *(p++)= *(c++);
     *(p++)= *(c++)&6;
    }
    (state->sign)[i++]= (*(c-7))>>4;
   }
   if ((state->numcrs)>chksiz){
    lp1= *(state->poly);
    sp= *(state->bilion);
    lngi= XCNT*XCNT;
    while (lngi--!=0){
     if (*lp1> cmpval){
      if ((*sp)++ == MAXBIL) (state->count)[1]=1;
      *lp1-= cmpval;
     }
     else if (*lp1< -cmpval){
      if ((*sp)-- == -MAXBIL -1) (state->count)[1]=1;
      *lp1+= cmpval;
     }
     ++lp1;
     ++sp;
    }
    /* if (stats>0){ */
/*      lseek (stats,(long)(strlen(nbuf)+1),0); */
/*      lngi=0; */
/*      i= XCNT-1; */
/*      while (count[i]==0) --i; */
/*      while (i!=0){ */
/*       write (stats,t,ntc((long)(i+1),t),state); */
/*       write (stats,"     ",5); */
/*       write (stats,t,ntc(count[i],t)+1,state); */
/*       lngi+= count[i--]; */
/*      } */
/*      write (stats,"\ntotal: ",8); */
/*      write (stats,t,ntc(lngi,t,state)+1); */
/*      write (stats,"\nbiggest stacked knots:\n",24); */
/*      i=12; */
/*      sp= (state->bstlst); */
/*      while (i--!=0) sp[i]=0; */
/*      lngi=notbeg-1; */
/*      while (lngi>0){ */
/*       k= (state->buf)[lngi]; */
/*       sp[++i]= k; */
/*       lngi-= k*8 + 8; */
/*       if (i==10) i= -1; */
/*      } */
/*      k= 0; */
/*      while (*sp!=0){ */
/*       j= ntc((long) *(sp++),t,state); */
/*       k+=j+1; */
/*       write (stats,t,j); */
/*       write (stats," ",1); */
/*      } */
/*      write (stats,"                                            ",44-k); */
/*      write (stats,"\n",1); */
/*     } */
   }
  }
  else (state->numcrs)=0;
 }
 else {
/* locate fake bigons, 0 branch bigons, and 3 cross circuits for (state->bstlst) */
  k= (state->numcrs);
  p= (state->crsbuf)[k]; /* only set k, other crossings will be set in their loops   */
  c= (state->sign)+k;    /* If a tangle is added to the knot with twists, I can find */
  while (k--!=0){ /* it sometimes and remove the twists! */
   n= (state->bstlst)[k];
   p-=8;
   j= *--c;
   i= *p;
   g= p[4];
   if (g==p[j]){   /* good fake bigon behind me */
    if (p[5]!=j){
     c= (state->crsbuf)[g];  /* wrong! twisted tangle ahead! remove twists! */
     mrecon (p+(j^4),c+ (p[j+1]^4),state);
     mrecon (p,c+ (p[5]^4),state);
     sqush2 (k,g,state);
     goto STEP1;
    }
    if (n>0) n= -4-n;
    else n+= -8;
   }
   if (p[j]==(state->crsbuf)[g][(state->sign)[g]^j^4] && (p[5]&2)!=0 && nopro!=(state->numcrs)){
    if ((p[j+1]&2)!=0){
     if (n<0) n-=12;
     else if (n>0) n= -8-n;  /* 0 branch has remov. bigon */
    }
    else { /* 0 branch has fake bigon */
     g=1;
     if (p[5]==j) g=2;  /* bad fake bigon 1 point, good is 2 points */
     if (n<0) n-=g;
     else if (n>0) n+=g;
    }
   }
   j^=4;
   if (p[j]==i){  /* good fake bigon ahead of me */
    if (p[1]!=j){
     q= (state->crsbuf)[i];  /* wrong! twisted tangle behind! remove twists! */
     mrecon (p+4,q+ (p[1]^4),state);
     mrecon (p+ *c,q+ (p[j+1]^4),state);
     sqush2 (k,i,state);
     goto STEP1;
    }
    if (n>0) n= -4-n;
    else n+= -8;
   }
   else if (p[j^4]==i){  /* bad fake bigon ahead */
    if (n<0) n-= 8;
    else if (n<8) n= 8;
    if ((state->bstlst)[i]<0) (state->bstlst)[i]-=8;     /* mutualize */
    else if ((state->bstlst)[i]<8) (state->bstlst)[i]=8;
   }
   if (p[j]==(state->crsbuf)[i][(state->sign)[i]^j^4] && (p[1]&2)!=0 && nopro!=(state->numcrs)){
    if ((p[j+1]&2)!=0){
     if (n<0) n-=12;
     else if (n>0) n= -8-n;  /* 0 branch has remov. bigon */
    }
    else { /* 0 branch has fake bigon (1 point) */
     g=1;
     if (p[1]==j) g=2;  /* bad fake bigon 1 point, good is 2 points */
     if (n<0) n-=g;
     else if (n>0) n+=g;
    }
   }
   (state->bstlst)[k]=n;
  }
  c= (state->stc);
  i= -8;
  if (nopro==(state->numcrs)) m=0;
  while (m>0){
   m-=3;
   n= *(c++);
   g= *(c++);
   p= (state->crsbuf)[n];  /* busting vertex of 3 cross circuit makes a removeable */
   sp= (state->bstlst)+n;  /* + or - link on the 0 branch -- free 2 reduction */
   k= p[*(c++)];
   j= p[g];
   if ((p[g|1]&2) != 0) p= (state->crsbuf)[j];
   else p= (state->crsbuf)[j]+(state->sign)[j];
   if (*p==k) --i;
   if (*sp<0) *sp+= i-4;
   else if (*sp==0) *sp= 4-i;
   else *sp= i- *sp;
  }
  i= nopro= (state->numcrs);
  k=n= 0;
  sp= (state->bstlst)+(state->numcrs);
  while (i--!=0){
   if (*--sp<n){
    n= *sp;
    j=i;
   }
   else if (*sp>k){
    k= *sp;
    m=i;
   }
  }
  if (n==0){
   if (k==0){
    j= dspair[2];
    nopro=0;
   }
   else j=m;
  }
  else if (n>-14 && k>12) j=m; /* since exponential growth is 1.41^n */
  ypow += bust(j,xpow,ypow,state);
  xpow^=512;
  ++(state->count)[(state->numcrs)-1];
 }
 // This point must be reached fairly often, so we're going to put in the timeout code here.
 end = clock();
 cpu_time_used = ((double)(end - start))/CLOCKS_PER_SEC;
 if (timeout > 0 && cpu_time_used > timeout) {
   free(outpoly);
   free(state);
   return NULL;
 }
 if ((state->numcrs)!=0) goto STEP1;
 goto NEWNOT;
}

void rmcir(short top,short vnum,short vbrnch,unsigned char *ovundr,int flag,lmpoly_state_t *state)
/* short top, vnum, vbrnch;
   unsigned char *ovundr;
   int flag; */
{
 register short i, j, k, h, m, n;
 register unsigned char *s, *c, *p;
 register long *lp1, *lp2;
 short g;
 if ((m=vbrnch) >=0){
  i= vnum;
  s= (state->crsbuf)[i];
  k= s[m];
  g= s[m|1];
  if (m!=0) c= s+4;
  else c= s+((state->sign)[i]^4);
  p= (state->crsbuf)[k]+g;
  h= *p= *c;
  n= p[1]= c[1];
  p= (state->crsbuf)[h]+n;
  *p= k;
  p[1]= g;
  s[6]= s[4]= s[2]= *s= i;
 }
 else ++(state->numlps);
 s= (state->clist);
 j= -1;
 while (top-1 > (i=k= ++j)){
  while (++i<top){
   if (s[i] < s[k]) k=i;
  }
  if (k!=j){    /* order crossings in (state->clist) */
   i= s[j];
   s[j]= s[k];
   s[k]=i;
   m= ovundr[j];
   ovundr[j]= ovundr[k];
   ovundr[k]= m;
  }
 }
      /* relink */
 m= top;
 while (m--!=0){
  i= s[m];
  g= ovundr[m];
  k= (state->crsbuf)[i][g];
  j= (state->crsbuf)[i][g|1];
  h= (state->crsbuf)[k][j]= (state->crsbuf)[i][g|4];
  n= (state->crsbuf)[k][j|1]= (state->crsbuf)[i][g|5];
  (state->crsbuf)[h][n]= k;
  (state->crsbuf)[h][n|1]= j;
 }
 if (flag!=0) return;
      /* squish crossings out of the list */
 i= *s;
 lp1= (long *) (state->crsbuf)[i];
 lp2=lp1+(state->suplng);
 p= (state->sign)+i;
 c=p+1;
 j= 0;
 k=1;
 while (++i!=(state->numcrs)){
  if (s[k]== i){
   lp2+= (state->suplng);
   ++c;
   if (++k==top) k=0;
  }
  else {
   *(lp1++)= *(lp2++);
   if ((state->suplng)==2) *(lp1++)= *(lp2++);
   *(p++)= *(c++);
  }
 }
     /* renumber branch pointers (subtract 1 if the crossing was after) */
 (state->numcrs)-= top;
 k= (state->numcrs);
 c= (state->crsbuf)[k];
 while (k--!=0){
  j=5;
  while (--j!=0){
   c-=2;
   i=top;
   s= (state->clist)+i;
   while (i--!=0){
    if (*c> *--s) --*c;
   }
  }
 }
 return;
}

int ntc(long i,char *buf,lmpoly_state_t *state)
{
 long j;
 int r;
 char *p;
 p=buf;
 r=0;
 if (i<0){
  i= -i;
  r=1;
  *(p++)= '-';
 }
 if (i>999999999) j=1000000000;
 else if (i>99999999) j=100000000;
 else if (i>9999999) j=10000000;
 else if (i>999999) j=1000000;
 else if (i>99999) j=100000;
 else if (i>9999) j=10000;
 else if (i>999) j=1000;
 else if (i>99) j=100;
 else if (i>9) j=10;
 else j=1;
 while (j!=0){
  ++r;
  *(p++)= (i/j)%10+0x30;
  j/=10;
 }
 *p='\n';
 return (r);
}

void delpow(lmpoly_state_t *state)
{
 register long *lp1, *lp2, *lsp, *osp;
 register short i, j, k, m, s, mpos;
 short numlnk, mpow, u, v, x, y, z, order[2], prslen[2], pstlen[2];
 lp1=(state->b);
 if ((state->numlps)<4){
  if ((state->numlps)==1){
   *(lp1++)= -1;
   *lp1= -1;
  }
  else if ((state->numlps)==3){
   *(lp1++)= -1;
   *(lp1++)= -3;
   *(lp1++)= -3;
   *lp1= -1;
  }
  else {
   *(lp1++)= 1;
   *(lp1++)= 2;
   *lp1= 1;
  }
 }
 else {
  *lp1= 1;
  (state->b)[1]= 4;
  (state->b)[2]= 6;
  k= (state->b)[3]= 4;
  (state->b)[4]= 1;
  while (k++!=(state->numlps)){
   j=k;
   (state->b)[j]=1;
   while (--j!=0) (state->b)[j]+= (state->b)[j-1];
  }
  if (((state->numlps)&1)!=0) while (--k>=0) (state->b)[k]= -(state->b)[k];
 }
 *prslen= *pstlen= prslen[1]= pstlen[1]= mpos= mpow= 0;
 osp= (state->b)+ (state->numlps)+1;
 if ((state->poslnk)>(state->neglnk)){
  order[1]= (state->poslnk);
  *order= (state->neglnk);
 }
 else {
  order[1]= (state->neglnk);
  *order= (state->poslnk);
  mpos=2;
 }
 z=2;
 while (z--!=0){
  numlnk=order[z];
  mpos^=2;
  while (numlnk--!=0){
   m= ++mpow;
   s=j=k= (state->numlps)+1;
   u= *pstlen;
   v= pstlen[1];
   ++prslen[z];
   *pstlen= x= *prslen;
   pstlen[1]= y= prslen[1];
   lp2= osp;
   lp1= lp2+ j+ m;
   if (z==0) lp1+= order[1];
   osp= lp1;
   while (j--!=0) *--lp1=0;
   while (m--!=0){
    i=j= s;
    lsp= lp1;
    if (mpos== (z|2) && y>0) *--lsp=0;
    if (mpos!=0) lp1+=k;
    else lp1+=i;
    while (j--!=0) *--lp1-= *--lp2;
    lp2+=i;
    lp1= lsp;
    *--lp1=0;
    while (i--!=0) *--lp1= *--lp2;
    lp1= lsp;
    lsp= lp2;
    lp2= lp1-1;
    j=s;
    while (j--!=0) *--lp1+= *--lp2;
    --lp1;
    if (x-->0) ++k;
    if (y-->0){
     ++k;
     if (mpos==z) *--lp1=0;
    }
    if (u-->0) ++s;
    if (v-->0) ++s;
    lp2=lsp;
   }
  }
 }
}

void addin(long kcoeff,short ypow,short xpow,int mult,int reach,lmpoly_state_t *state)
/*short ypow, xpow;
long kcoeff;
int mult, reach;*/
{
  register long *lp1, *lp2 = NULL, addon, *p; // Added init for lp2, JC.
  register short i, j, k, numlnk, loops;
  i= (state->poslnk);
  k= (state->neglnk);
  numlnk= i+k;
  loops= (state->numlps)+ numlnk;
  ypow+= XCNT-xpow-1+(state->lowx)+ 2*(k-i);
  xpow-=(state->lowx)+loops;
  p=(state->b);
  while (numlnk-->=0){
    if (reach!=0) lp2= (state->poly)[xpow+reach]+ypow-reach;
    lp1= (state->poly)[xpow]+ ypow;
    j= loops;
    while (j-->=0){
      addon= kcoeff* *(p++);
      *lp1+= addon;
      lp1+=2;
      if (reach!=0){
	*lp2+= mult*addon;
	lp2+=2;
      }
    }
    if (--i<0) --loops;
    if (--k<0) --loops;
    else ypow-=2;
    xpow+=2;
  }
}

void untwst(short l,short n,short i,short ex,lmpoly_state_t *state)
/* short l, n, i, ex; */
{
 register unsigned char *pl, *pn;
 register short g, h, m, k, d, e;
 short j;
 j= (state->sign)[i]; /* flip 3 crossing circuits */
 k= j^ex;
 ex^=4;
 m= j^ex;
 if ((state->sign)[l]-ex==2){
  g= j^2;
  h= g^4;
  j= 6;
 }
 else {
  h= j^2;
  g= h^4;
  j= 2;
 }
 pl= (state->crsbuf)[l];
 pn= (state->crsbuf)[n];
 d= pl[m];
 e= pl[m|1];
 pl[m]= pn[g];
 pl[m|1]= pn[g|1];
 if (d!=n){
  pn[g]= d;
  pn[g|1]= e;
 }
 else {
  pn[g]= l;
  pn[g|1]= k;
 }
 (state->crsbuf)[pl[m]][pl[m|1]]= l;
 (state->crsbuf)[pl[m]][pl[m|1]|1]= m;
 d= pn[h];
 e= pn[h|1];
 pn[h]= pl[k];
 pn[h|1]= m= pl[k|1];
 if (d!=l){
  pl[k]= d;
  pl[k|1]= e;
 }
 else {
  pl[k]= n;
  pl[k|1]= g;
 }
 (state->crsbuf)[pn[h]][m]= n;
 (state->crsbuf)[pn[h]][m|1]= h;
 h= pn[g|1];
 (state->crsbuf)[pn[g]][h]= n;
 (state->crsbuf)[pn[g]][h|1]= g;
 d= pl[k];
 e= pl[k|1];
 (state->crsbuf)[d][e]= l;
 (state->crsbuf)[d][e|1]= k;
 pl[ex]= g= (state->crsbuf)[i][ex];
 pl[ex|1]= h= (state->crsbuf)[i][ex|1];
 (state->crsbuf)[g][h]= l;
 (state->crsbuf)[g][h|1]= ex;
 pn[j]= g= (state->crsbuf)[i][k];
 pn[j|1]= h= (state->crsbuf)[i][k|1];
 (state->crsbuf)[g][h]= n;
 (state->crsbuf)[g][h|1]= j;
}

int conchk(lmpoly_state_t *state)
{
 register short i, j;
 int vnum, vbrnch, top;
 char str[20], local_cbuf[XCNT*2];
 vnum=vbrnch= 0;
 (state->lowx)= 1;
 i=top= (state->numcrs)*2;
 while (i>=0) local_cbuf[i--]=0;
 i=0;
 while (i!=top){
  --(state->lowx);
  while (local_cbuf[i]==0){
   local_cbuf[i]=1;
   i=vnum;
   j=vbrnch;
   vnum= (state->crsbuf)[i][j];
   vbrnch= (state->crsbuf)[i][j|1];
   if (vnum>=(state->numcrs)){
    write (1,"error at crossing #",19);
    write (1,str,ntc((long)(i+1),str,state));
    write (1,", knot can't have a crossing #",30);
    write (1,str,ntc((long)(vnum+1),str,state)+1);
    return (1);
   }
   if ((state->crsbuf)[vnum][vbrnch]!=i || ((state->crsbuf)[vnum][vbrnch|1])!=j){
    write (1,"the possible connection from branch ",36);
    write (1,str,ntc((long)(i+1),str,state));
    *str= 'a'+(j>>1);
    write (1,str,1);
    write (1," to branch ",11);
    write (1,str,ntc((long)(vnum+1),str,state));
    *str= 'a'+(vbrnch>>1);
    write (1,str,1);
    write (1," is inconsistent\n",17);
    return (1);
   }
   if (vbrnch==0){
    write (1,"a-c part of crossing #",22);
    write (1,str,ntc((long)(vnum+1),str,state));
    write (1," is oriented backwards\n",23);
    write (1,"with respect to crossing #",26);
    write (1,str,ntc((long)(i+1),str,state)+1);
    return (1);
   }
   else if (vbrnch== ((state->sign)[vnum] & 6)){
    write (1,"b-d part of crossing #",22);
    write (1,str,ntc((long)(vnum+1),str,state));
    write (1," is oriented backwards\n",23);
    write (1,"with respect to crossing #",26);
    write (1,str,ntc((long)(i+1),str,state)+1);
    return (1);
   }
   i= (vnum<<1) + ((vbrnch>>1)&1);
   vbrnch^=4;
  }
  i=0;
  while (local_cbuf[i]!=0) ++i;
  vnum= i/2;
  if ((i&1)!=0) vbrnch= (state->sign)[vnum] & 6;
  else vbrnch=0;
 }
 return(0);
}

void mypause(int sig)
{
  //(state->restrt)=1;
 //signal(SIGTERM,pause);  /* pause is now a system function */
 signal(sig,mypause);
}

int bust(short tobust,short xpow,short ypow,lmpoly_state_t *state)
/* short tobust, xpow, ypow; */
{
 register unsigned char *fastp, *c;
 register long *lp1, *lp2;
 register short i, j, m, k, l, n, *sp;
 unsigned char *inbuf[XCNT], tsign[XCNT], *bstcrs;
 int dir;
 //long lngi;
 bstcrs= (state->crsbuf)[tobust];
 fastp=(state->sign);
 c=tsign;
 lp2= (long *) *(state->crsbuf);
 lp1= (long *) ((state->buf)+(state->notbeg));
 j=0;
 while (j!=(state->numcrs)){
  if (j==tobust) ++fastp;
  *(c++)= *(fastp++);
  inbuf[j++]= (unsigned char *) lp1;
  if ((state->suplng)==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
 }
 fastp= bstcrs;
 l= (state->sign)[tobust]/2;
 dir= 2-l;
 m= fastp[7-l];
 n= fastp[8-l];
 j= fastp[3*l-3];
 k= fastp[3*l-2];
 inbuf[j][k]= m;
 inbuf[j][k|1]= n;
 inbuf[m][n]= j;
 inbuf[m][n|1]= k;
 m= fastp[5-l];
 n= fastp[6-l];
 j= fastp[3-l];
 k= fastp[4-l];
 inbuf[j][k]= m;
 inbuf[j][k|1]= n;
 inbuf[m][n]= j;
 inbuf[m][n|1]= k;
 k= (state->numcrs)-1;
     /* squish crossing out of the list (tsign has already been squished) */
 i= k-tobust;
 lp1= (long *) inbuf[tobust];
 lp2=lp1+(state->suplng);
 while (i--!=0){
  if ((state->suplng)==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
 }
     /* renumber branch pointers (subtract 1 if the crossing was after) */
 (state->notbeg)+= (state->numcrs)*8;
 fastp= *inbuf;
 while (++i!=k){
  fastp[1]|= tsign[i]<<4;
  j=5;
  while (--j!=0){
   if (*fastp>=tobust) --*fastp;
   fastp+=2;
  }
 }
 fastp= (state->buf)+(state->notbeg)-7;
 *fastp= *(fastp-1)=0;
 if (ypow+dir<0){
  *fastp=1;
  (state->buf)[(state->notbeg)-2]= -ypow-dir;
 }
 else (state->buf)[(state->notbeg)-2]= ypow+dir;
 if ((xpow&512) ==0) (state->buf)[(state->notbeg)-8]=1;
 *++fastp= (state->neglnk);
 *++fastp= (state->poslnk);
 *++fastp= (state->numlps);
 *++fastp= xpow+1;
 (state->buf)[(state->notbeg)-1]= i;
 fastp= bstcrs;
 sp= (short *) fastp;
 j= *sp;
 *sp= sp[l];
 sp[l]= sp[2];
 sp[2]= sp[l^2];
 sp[l^2]= j;
 (state->sign)[tobust]^=4;
 (state->crsbuf)[*fastp][fastp[1]|1]=0;
 (state->crsbuf)[fastp[2]][fastp[3]|1]=2;
 (state->crsbuf)[fastp[4]][fastp[5]|1]=4;
 (state->crsbuf)[fastp[6]][fastp[7]|1]=6;
 return (dir*2);
}

void mrecon(unsigned char *c,unsigned char *p,lmpoly_state_t *state)
/*unsigned char *c, *p;*/
{
 register short j, k, m, n;
 register unsigned char *ccb;
 k= *c;
 j= c[1];
 ccb= (state->crsbuf)[k]+j;
 n= *ccb= *p;
 m= ccb[1]= p[1];
 ccb= (state->crsbuf)[n]+m;
 *ccb= k;
 ccb[1]= j;
}

void squish(short sqush,lmpoly_state_t *state)
/* short sqush; */
{
 register short i, j, tosqsh;
 register unsigned char *p, *c;
 register long *lp1, *lp2;
       /* squish crossing out of the list */
 tosqsh= sqush;
 j= --(state->numcrs)-tosqsh;
 lp1= (long *) (state->crsbuf)[tosqsh];
 lp2=lp1+(state->suplng);
 p= (state->sign)+tosqsh;
 c=p+1;
 while (j--!=0){
  if ((state->suplng)==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
  *(p++)= *(c++);
 }
    /* renumber branch pointers (subtract 1 if the crossing was after) */
 i= (state->numcrs);
 p= (state->crsbuf)[i];
 while (i--!=0){
  j=5;
  while (--j!=0){
   p-=2;
   if (*p>=tosqsh) --*p;
  }
 }
}

void sqush2(short n,short g,lmpoly_state_t *state)
/* short n, g; */
{
 long *lp1, *lp2;
 register unsigned char *p, *c;
 register short m, i, j;
     /* squish crossings out of the list */
 i= n;
 m= g;
 if (m>i){
  j=m;
  m=i;
  i=j;
 }
 lp1= (long *) (state->crsbuf)[m];
 lp2=lp1+(state->suplng);
 p=(state->sign)+m;
 c=p+1;
 j= i-m;
 while (--j!=0){
  if ((state->suplng)==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
  *(p++)= *(c++);
 }
 lp2+=(state->suplng);
 ++c;
 j= (state->numcrs)-i;
 (state->numcrs)-=2;
 while (--j!=0){
  if ((state->suplng)==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
  *(p++)= *(c++);
 }
    /* renumber branch pointers (subtract 1 if the crossing was after) */
 p= (state->crsbuf)[(state->numcrs)];
 n= (state->numcrs);
 while (n--!=0){
  j=5;
  while (--j!=0){
   p-=2;
   if (*p>=i) --*p;
   if (*p>=m) --*p;
  }
 }
}

void triple(short m,short i,unsigned char *p,lmpoly_state_t *state)
/*short m, i;
  unsigned char *p; */
{
 register unsigned char *crsi, *crsk;
 register short j, d, e, n, h, k;
 crsi= p;
 n= m^4;
 k= crsi[m];
 crsk= (state->crsbuf)[k];
 j= crsi[m|1];
 d= crsk[j]= crsi[n];
 e= crsk[j|1]= crsi[n|1];
 (state->crsbuf)[d][e]= k;
 (state->crsbuf)[d][e|1]= j;
 j^=4;
 d= crsi[m]= crsk[j];
 e= crsi[m|1]= crsk[j|1];
 (state->crsbuf)[d][e]= i;
 (state->crsbuf)[d][e|1]= m;
 crsi[n]= k;
 crsi[n|1]= j;
 crsk[j]= i;
 crsk[j|1]= n;
 j^=2;
 n= crsk[j];
 h= crsk[j|1]^4;
 d= crsk[j]= (state->crsbuf)[n][h];
 e= crsk[j|1]= (state->crsbuf)[n][h|1];
 (state->crsbuf)[d][e]= k;
 (state->crsbuf)[d][e|1]= j;
 j^=4;
 n= crsk[j];
 h= crsk[j|1]^4;
 d= crsk[j]= (state->crsbuf)[n][h];
 e= crsk[j|1]= (state->crsbuf)[n][h|1];
 (state->crsbuf)[d][e]= k;
 (state->crsbuf)[d][e|1]= j;
 j=m^2;
 n= crsi[j];
 h= crsi[j|1]^4;
 d= crsi[j]= (state->crsbuf)[n][h];
 e= crsi[j|1]= (state->crsbuf)[n][h|1];
 (state->crsbuf)[d][e]= i;
 (state->crsbuf)[d][e|1]= j;
 j^=4;
 m= crsi[j];
 h= crsi[j|1]^4;
 d= crsi[j]= (state->crsbuf)[m][h];
 e= crsi[j|1]= (state->crsbuf)[m][h|1];
 (state->crsbuf)[d][e]= i;
 (state->crsbuf)[d][e|1]= j;
 sqush2(n,m,state);
}

int tstcir(short h,short *bstlst,short *dspair,short *skflag,lmpoly_state_t *state)
/* short h, *bstlst, *dspair, *skflag; */
{
 register short d, g, i, j, m, length, *q, *sp;
 register unsigned char *p, *tp = NULL, *vp; // Added init for tp. JC
 short k = 0, n, lngsum, badcrs = 0, vertex, lngpos = 0, i1, i2, i3, j1, j2, j3;
 sp=q= (state->tt);  /* # of overs and unders were passed in through *(state->tt) and [2] */
 lngsum= n= 0;
 if (*(state->tt) > (state->tt)[2]) n=2;
 else if (*(state->tt) == (state->tt)[2]) lngsum=1;   /* need to test BOTH overs & unders */
 vp= (state->clist);
 vertex=h;
 while (vertex-->0) *(sp++)= *(vp++);   /* make a working copy of (state->clist) */
 i= h-1;
 if ((h&1)!=(length=0)) vertex= (state->clist)[i];
 else {
  if ((state->donlnk)[XCNT]==0){   /* if my check-off list is now invalid */
   d= (state->numcrs);             /* re-set it all to 0 */
   while (d-->0) (state->donlnk)[d]= 0;
   (state->donlnk)[XCNT]=1;
  }
  if (((state->donlnk)[(state->clist)[i]]&((state->t)[i]+2))!=0) return(0);
 }
 m=j= h/2;
 *sp= -1;
 p=state->t;
 while (j-->length){
  i=2;
  while (*q==1024){
   ++q;
   ++p;
  }
  d= *q;
  if (*(p++)==0) i+= (state->sign)[d];
  g= (state->crsbuf)[d][i&6];
  sp= ++q;
  while (*sp>=0 && *sp!=g) ++sp;
  if (*sp<0) j= -10;
  else if (*sp==vertex) length= -1;
  else *sp= 1024;
 }
 if (j== -11){
  d=h;
  sp=q= (state->tt);
  vp= (state->clist);
  while (d-->0) *(sp++)= *(vp++);
  length=0;
  p=state->t;
  while (m-->length){
   i=6;
   while (*q==1024){
    ++q;
    ++p;
   }
   d= *q;
   if (*(p++)==0) i+= (state->sign)[d];
   g= (state->crsbuf)[d][i&6];
   sp= ++q;
   while (*sp>=0 && *sp!=g) ++sp;
   if (*sp<0) m= -10;
   else if (*sp==vertex) length= -1;
   else *sp= 1024;
  }
  j-=m;
 }
 if (vertex< (length=lngsum=d=m= 0)){  /* if a link */
  tp=p= state->t;
  while ((n^*(p++))==0) ++d;  /* find first bad crossing and move it to the */
  vp= (state->clist)+d;                /* front so that (state->clist) cannot break a string */
  --p;                        /* of good crossings in a link from */
  sp= (state->tt);                     /* (state->clist)[h-1] to (state->clist)[0] */
  i=m= h-d;
  while (m-->0){
   (state->donlnk)[*vp]|= *(p++)+2;  /* also check off all list elements to assure */
   *(sp++)= *(vp++);        /* this link never comes to tstcir again with */
  }                         /* a different beginning crossing */
  m=d;
  vp= (state->clist);
  while (m-->0){
   (state->donlnk)[*vp]|= *tp+2;
   *(sp++)= *vp;
   *(vp++)= *(tp++);
  }
  p=state->t;
  while (i-->0) *(p++)= *(tp++);
  vp= (state->clist);
  while (d-->0) *(p++)= *(vp++);
  d=m= 0;
 }
 else {
  sp= (state->tt);
  vp= (state->clist);
  i= h;
  while (i-->0) *(sp++)= *(vp++);
 }
 p=state->t;
 g= h|1;         /* find the longest good string & put in "length" */
 (state->t)[g]= (state->t)[g-1];   /* also, the longest with only 1 bad crossing in "lngsum" */
 (state->t)[g-1]= n^2;    /* badcrs is the bad crossing in lngsum */
 while (++i<g){
  if ((n^*(p++)) ==0) ++d;
  else {
   if (d>=length){
    length=d;
    k= i-d;
   }
   if (m+d>=lngsum){
    lngsum= m+d+1;
    if (d==0) tp= p-3;  /* tp points to a generic good crossing in lngsum */
    else tp= p-2;
    if (d!= (lngpos=i)) badcrs= (state->tt)[i-d-1];
    else badcrs= (state->tt)[i];
   }
   m=d;
   d=0;
  }
 }
 i= h/2;
 if (vertex>=0 && length>=i){ /* can this circuit untwist its vertex? */
   if (k==0 && (*(state->t)^*p)!=0) return (skinny(1,h,(short) 0,state));
   else if (k+length+1==h && ((state->t)[k]^*p)==0) return (skinny(1,h,i,state));
 }
 if (j==0 && length>=i){     /* circuit is not skinny, and can be skinnied */
  if (vertex>= (i1=0) && length==i){
   m= k+i;
   *(p-1)= *p;     /* infinite loops happen when the loose strand is shared */
   d= (state->tt)[m-1];        /* between 2 intersecting circuits -- skflag prevents */
   while (m>=k){           /* looping, allows me to skinny once per bust!   */
    vp= (state->crsbuf)[(state->tt)[m]];
    if ((state->t)[m]!=0) j=2;
    if (vp[j]==d || vp[j^4]==d) m=i1= -1;  /* possible inf. loop circuit! */
    else j=0;
    m-=i;
    if (k==0) d= (state->tt)[h-1];
    else d=(state->tt)[k-1];
   }
  }
  if (i1== (j=0)) return (skinny(0,h,k,state)); /* OK to skinny it! */
  if (*skflag==0){
   *skflag=1;
   return (skinny(0,h,k,state)); /* OK to do ONE skinny */
  }
 }
 if (lngsum>=i){
  if (lngsum == h-2) lngsum=h;   /* a link w/only 2 crossings vanishes! */
  d= (lngsum-i)*8;  /* and find all the extra benefits of badcrs */
  if (j==0) ++d;    /* give one point for doing a skinny */
  if (vertex>=0){
   if (lngpos==lngsum && (*tp^*p)!=0) d+=4;
   else if (lngpos+1==h && (*tp^*p)==0) d+=4;
   else if (lngsum+1==h) d+=4;
  }
  if (d==8 && length==lngsum-1) d=0; /* this is really just a fake bigon */
  if (j==0 && d!=0 && (state->numcrs)+ 2*(i-lngsum)>11) d+=2;
  if ((k=lngpos-lngsum)!= (m=0) && d!=0){
   if ((state->t)[k]!=0) m= (state->sign)[(state->tt)[k]];
   if ((state->crsbuf)[(state->tt)[k]][m] == (state->tt)[k-1]) ++d;
   if (k!=1 && ((state->t)[k-2]^(state->t)[k-1])!=0) ++d;
  }
  if (lngpos+1 != g && d!=0){
   if ((state->t)[lngpos]!= (m=0)) m= (state->sign)[(state->tt)[lngpos]];
   if ((state->crsbuf)[(state->tt)[lngpos]][m] == (state->tt)[lngpos-1]) ++d;
   if (lngpos+2 !=g && ((state->t)[lngpos+1]^(state->t)[lngpos])!=0) ++d;
  }
  if (d<6 && j==0){    /* min. one point for making skinny circuits */
   if (bstlst[badcrs]==0) bstlst[badcrs]= *dspair= 1;
  }
  else {
   *dspair=1;
   (state->donlnk)[badcrs]=0;  /* from badcrs's perspective, link may be better */
   if (bstlst[badcrs]<0){
    if (bstlst[badcrs]> -8-d) bstlst[badcrs]= -8-d;
   }
   else if (bstlst[badcrs]<d) bstlst[badcrs]=d;
  }
 }
 else {
  k= h/4;
  if (*dspair<k) k= *dspair;
  m=1;
  while (m++<k){
   tp=state->t;
   *(state->gapsto)=j= 0;
   d= m+1;
   while (--d>0){
    while ((*(tp++)^n) ==0) ++j;
    (state->gapsto)[d]= j++;
   }
   while (j<g){
    while ((*(tp++)^n) ==0) ++j;
    length= (j-(state->gapsto)[d]-i)*4 +2;
    if (length>0 && (m<*dspair || (m== *dspair && length>=dspair[1]))){
     if (m<*dspair || length>dspair[1]) dspair[3]=0;
     *dspair=m;
     dspair[1]= length;
     length=lngpos= 0;
     if (d!=m) lngpos= d+1;
     j1= j;
     j2= (state->gapsto)[lngpos];
     h=m;
     while (h-->0){
      if (++lngpos>m) lngpos=0;
      vertex= (state->tt)[j2];
      j3= (state->gapsto)[lngpos];
      lngsum= j1-j3;
      p= (state->crsbuf)[vertex];
      i1= (state->sign)[vertex];
      i2= p[i1];
      i3= p[i1|1];
      i1^=4;
      if ((state->crsbuf)[i2][(i3+i1)&6]== p[4]){
	//if ((p[5]^i3&2) !=0) lngsum= -1 -lngsum;
	/* According to wikipedia, the & should come first, so this is */
	if ((p[5]^(i3&2)) !=0) lngsum= -1 -lngsum;
       else lngsum= -10 -lngsum;
      }
      i2= p[i1];
      i3= p[i1|1];
      if ((state->crsbuf)[i2][(i3+i1)&6]== *p){
	if ((p[1]^(i3&2)) !=0){  // The same is true here
        if (lngsum<0) lngsum-=5;
        else lngsum= -1 -lngsum;
       }
       else if (lngsum<0) lngsum-= 10;
       else lngsum= -10 -lngsum;
      }
      if (lngsum<0){
       if (length>lngsum){
        length= lngsum;
        badcrs= vertex;
       }
      }
      else if (length>=0 && length<lngsum){
       length= lngsum;
       badcrs= vertex;
      }
      j1=j2;
      j2=j3;
     }
     i1= dspair[3];
     if (length<0){
      if (i1>length){
       dspair[2]= badcrs;
       dspair[3]= length;
      }
     }
     else if (i1>=0 && i1<length){
      dspair[2]= badcrs;
      dspair[3]= length;
     }
     k=0;
    }
    (state->gapsto)[d--]= j++;
    if (d<0) d=m;
   }
  }
 }
 return(0);
}

int skinny(int twist,short lencir,short lngpos,lmpoly_state_t *state)
/* short lencir, lngpos;
   int twist; */
{
 register short a, b, g, n, vnum, dir, *c;
 register unsigned char *p, *q, *cp;
 short i, j, k, m, last, f, v2; //*sp;
 vnum= -1;     /* assume circuit is a link */
 dir= 6;       /* and I will be placing loose strand LEFT of fixed strand */
 if ((lencir&1)!= (b=0)){
  vnum= (state->tt)[lencir-1];                /* if not a link */
  if (((state->t)[lencir-1]=(state->t)[lencir]) ==0) b= (state->sign)[vnum];
    /* tstcir wrecks (state->t)[lencir-1] on non-links but saves a copy! */
 }
 if (vnum>=0){
  g=vnum;      /* SET DIR PROPERLY!! */
  while (b>=0){        /* if VERTEX is left of fixed strand, set dir RIGHT */
   p=(state->crsbuf)[g]+b;  /* there is also a special case if I am untwisting */
   g= *p;
   b= *(p+1)^4;  /* walk out from the vertex -- do I run into the circuit? */
   a=lencir;
   c= (state->tt);
   while (a>0 && *(c++)!=g) --a;
   if (a>0){
    if (*--c!=vnum){
     if (b==0) dir= (state->sign)[*c]^4;     /* yes */
     else dir=b;
    }
    else {
     g= *(state->tt);     /* no -- I ran back into the vertex again */
     b=6;              /* trace the first left region of the circuit */
     if (*(state->t)==0) b=(state->sign)[g]-2;
     n= (state->tt)[lencir-2];  /* does it connect to the other side of the circuit? */
     while (g!=vnum && g!=n){
      p= (state->crsbuf)[g]+b;
      g= *p;
      b= (*(p+1)+2)&6;
     }
     if (g==vnum) dir=2;  /* no -- vertex is left of fixed strand */
    }
    b= -1;
   }
  }
  if (twist!=0) dir^=4;
 }
 a=j= lencir/2;   /* to move strand, remove it - then reinsert crossings by */
 p= (state->clist);       /* hand.  Put moveable crossings in (state->clist) & call rmcir */
 c= (state->tt)+lngpos;
 while (a-->0) *(p++)= *(c++);
 --(state->numlps);                   /* rmcir adds a (state->numlps) for fun */
 rmcir (j,(short) 0,a,(state->t)+lngpos,1,state);
 i=k= lngpos+j;
 if (vnum>=0){
  --lencir;
  if (lngpos==0 && twist==0) twist= -1;
 }
 else if (i==lencir) i=0;
 if (twist> (m=b= 0) && lngpos!=0){
  if ((state->t)[lencir]==0) b= (state->sign)[vnum];
  last= (state->crsbuf)[vnum][b];   /* set last properly if I am untwisting vertex */
  f= (state->crsbuf)[vnum][b|1];
 }
 else {
  last= (state->tt)[i];  /* set last to automatically hook up one of the bigon ends */
  f= 4;
  if ((state->t)[i]==0) f= (state->sign)[last]^4;
 }
 if ((state->t)[--k]==0) m=4;
 while (j-->0){
  if (i==lencir) i=0;
  g=dir;
  vnum= (state->tt)[i];   /* vnum is the crossing on the fixed strand */
  b= (state->sign)[vnum];
  if ((state->t)[i++]==0) g+=b;
  else b^=4;
  p= (state->crsbuf)[last]+f;  /* reconnect the loose end from the last loop to the */
  *(p++)= n= (state->tt)[k--];  /* present n */
  q= (state->crsbuf)[n];   /* n is the crossing on the loose strand - q points to it */
  (state->sign)[n]= a= b^m;   /* this calculation for (state->sign)[n] is a 2 step trick */
  if (m==0){
   a=0;
   b=dir;
  }
  else b^=dir;
  *p=a;
  q[a]= last;   /* reconnect last going FOREWARDS along loose strand from n */
  q[a+1]=f;
  f= a^4;
  g&=6;
  cp= (state->crsbuf)[vnum]+g;
  v2= *cp;             /* v2 is one step out the dir branch of fixed strand */
  *(cp++)= last= n;
  a= *cp;    /* reconnect vnum to n instead (insert n) */
  *cp=b;
  p= q+b;
  *(p++)= vnum;   /* mutual from n to vnum */
  *p=g;
  b^=4;
  p= (state->crsbuf)[v2]+a;  /* insert n from v2 side */
  *(p++)= n;
  *p=b;
  p= q+b;     /* mutual connect from n */
  *(p++)= v2;
  *p=a;
 }
 if (twist!=0){       /* reconnect the loose end properly */
  vnum=(state->tt)[lencir];
  if ((state->t)[lencir]!= (dir=b=0)) dir=b=(state->sign)[vnum];
  if (twist>0){
   p= (state->crsbuf)[vnum];
   if (lngpos==0){
    n= p[b^4];
    b= p[b^5];
    if ((state->t)[lencir]== (dir=0)) dir=(state->sign)[vnum];
   }
   else {
    n= (state->tt)[--i];
    if ((state->t)[i]== (b=0)) b=(state->sign)[n];
   }
   a= p[dir];
   g= p[dir|1];
   k= p[dir^5];
   dir= p[dir^4];
   (state->crsbuf)[a][g]= dir;
   (state->crsbuf)[a][g|1]=k;
   (state->crsbuf)[dir][k]=a;
   (state->crsbuf)[dir][k|1]=g;
  }
  else twist=0;
 }
 else if ((state->t)[--i]== (b=0)) b= (state->sign)[(state->tt)[i]];
 if (twist==0) n= (state->tt)[i];
 (state->crsbuf)[last][f]= n;
 (state->crsbuf)[last][f|1]= b;
 (state->crsbuf)[n][b]= last;
 (state->crsbuf)[n][b|1]= f;
 if (twist!=0) squish(vnum,state); /* if I untwisted vertex, now eliminate it */
 return(1);
}
