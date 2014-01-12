/*          FAST VERSION       */
/* A load balanced calculation of the oriented polynomial  Ver 1.1f */
/* Bruce Ewing and Kenneth C. Millett  (author: Ewing)    Feb 5, 01 */
/* Department of Mathematics */
/* University of California */
/* Santa Barbara, CA 93106 */
/* there is a hard limit of 255 crossings in this program */

#include <signal.h>
#include <sys/types.h>
#include <sys/times.h>

  /* values the USER needs to define!! */
#define XCNT       250     /* maximum crossings in knot (at most 255) */
#define XCNTSQ   62500     /* YOU MUST CALCULATE XCNT SQUARED! */
#define XCSQTR  187500     /* YOU MUST CALCULATE XCNTSQ times THREE! */


/* system values -- it's best if these are not touched */
#define EOFCHR       0     /* integer value of my "end of file" character */
#define MAXBIL   32767     /* positive value where overflow warning triggers */

long plybuf[XCNTSQ], *poly[XCNT], notbeg, count[XCNT], b[XCNT+1];
unsigned char sign[XCNT], donlnk[XCNT+1], t[XCNT+2], crsbuf[XCNT+1][8];
unsigned char buf[XCSQTR], cbuf[10242], clist[XCNT+2], stc[XCNT*2];
short numcrs, numlps, poslnk, neglnk, lowx, restrt, gapsto[65];
short tt[XCNT+2], bstlst[XCNT], bilbuf[XCNTSQ], *bilion[XCNT], suplng;

main(argc,argv)
int argc;
char *argv[];
{
 register short i, j, k, h, m, n, *sp;
 register unsigned char *p, *c;
 int in, out, stats, len, len2, pause();
 long lngi, *lp1, *lp2, kstrt, cmpval;
 short g, xpow, ypow, dspair[4], maxcrs, chksiz, nopro, skflag;
 unsigned char nbuf[82], *q, *s;
 struct tms hi;
 maxcrs=XCNT;    /* maximum crossings in a knot, including all limits */
 chksiz= 50;     /* knot size trigger where bilion is checked, stats printed */
 *sign= 255;      /* do 3 tests on the machine to verify trick operation */
 i= (short) *sign;
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
 *t=1;
 t[1]=2;
 t[2]=3;
 t[3]=4;
 lp1= (long *) t;
 lp2= (long *) nbuf;
 *lp2= *lp1;
 if (*nbuf!=1 || nbuf[1]!=2 || nbuf[2]!=3 || nbuf[3]!=4){
  write (1,"4 byte copy problem -- use PORTABLE PROGRAM version!\n",53);
  exit (1);
 }
 if (sizeof(lngi)== (suplng=2)){
  write(1,"halt/restart won't work, and watch for coefficient overflows\n",61);
  cmpval= 10000;
  chksiz= 20;
 }
 else cmpval= 1000000000;     /* how big a value fits in a poly memory slot? */
 if (sizeof(lngi)== 8) suplng= 1;
 if (argc<2){
  write (1,"usage:  lmpoly -rs knotfile [knotfile ...]\n",43);
  exit (1);
 }
 if ((out=open("lmknot.out",1))== -1){
  if ((out=creat("lmknot.out",0644))== -1){
   write (1,"fatal: cannot create output file 'lmknot.out'\n",46);
   exit (1);
  }
 }
 else lseek (out,(long) 0,2);
 close (out);
 *count= restrt= i= j= kstrt= nopro= 0;
 stats= -1;
 if (*argv[1]=='-'){
  if (argc==2 && argv[1][1]=='r'){
   restrt=1;
   if (argv[1][2]=='s') stats= 0;
  }
  else if (argv[1][1]=='s'){
   stats= 0;
   ++argv;
   --argc;
  }
 }
 signal(SIGTERM,pause);
 lp1= plybuf;
 sp= bilbuf;
 while (j<XCNT){
  poly[j]= lp1;
  bilion[j]= sp;
  lp1+= 1+ 2* (XCNT- ++j);
  sp+= 1+ 2* (XCNT-j);
 }
NEWFIL:
 if (stats>0){     /* remove stats info for previous knot file */
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
 }
 if (--argc==0) exit(0);
 if (restrt!=0){           /* if I am restarting an old calculation */
  if ((in=open("lmpoly.restrt",0))== -1){
   write (1,"couldn't find restart file\n",27);
   exit (0);
  }
  read (in,nbuf,82);
  read (in,cbuf,12);
  *argv= (char *) cbuf;
  read (in,buf,7);
  numcrs= *buf;
  xpow= (short) buf[2];
  ypow= (short) buf[1];
  numlps= (short) buf[3];
  poslnk= (short) buf[4];
  neglnk= (short) buf[5];
  if ((buf[6]&1) !=0) ypow= -ypow;
  if ((buf[6]&2) !=0) xpow|= 512;
  read (in,*crsbuf,numcrs*8);
  i= 0;
  while (i!=numcrs){
   sign[i]= crsbuf[i][1]>>4;
   crsbuf[i++][1]&=6;
  }
  p= (unsigned char *) &lowx;
  read (in,p,2);
  p= (unsigned char *) &kstrt;
  read (in,p,4);
  p= (unsigned char *) &notbeg;
  read (in,p,4);
  read (in,buf, (int) notbeg);
  read (in,(char *) count,XCNT*4);
  read (in,(char *) plybuf,4*XCNT*XCNT);
  read (in,(char *) bilbuf,2*XCNT*XCNT);
  close (in);
  in=0;
  s= cbuf+14;
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
 if (restrt!=0){
  if (stats>0){
   write (stats,nbuf,strlen(nbuf));
   write (stats,"\n",1);
  }
  restrt=0;
  goto STEP1; /* restart file has set up all NEWNOT data already */
 }
 c=cbuf;
 c[read(in,c,10240)]= EOFCHR;
 while (*c=='\n') ++c;
 s=c;
NEWNOT:
 i= XCNT-1;
 while (i!=0 && count[i]==0) --i;
 if (*count!= (numcrs=0)){
  if (stats>=0){
   if ((out=open("lmknot.stats",1))== -1) out=creat("lmknot.stats",0644);
   else lseek (out,(long) 0,2);
  }
  if (stats>=0 && out>0){
   write (out,nbuf,strlen(nbuf));
   write (out,"\n",1);     /* write out statistics for knot just completed */
   lngi=0;
   while (i!=0){
    write (out,t,ntc((long)(i+1),t));
    write (out,"     ",5);
    write (out,t,ntc(count[i],t)+1);
    lngi+= count[i--];
   }
   write (out,"\ntotal: ",8);
   write (out,t,ntc(lngi,t)+1);
   write (out,"run: ",5);
   times(&hi);
   lngi= ((hi.tms_utime-kstrt)*5)/3;
   write (out,t,ntc(lngi/100,t));
   write (out,".",1);
   len= ntc(lngi%100,t);
   write (out,"0",2-len);
   write (out,t,len);
   write (out," s\n\n",4);
   kstrt= hi.tms_utime;
   if (count[2]<0) write (out,"program error: knot became inconsistent\n",40);
   close (out);
  }
  k= XCNT;
  i=0;
  while (k!=0 && i==0){
   i= 2+ 2*(XCNT-k);
   lp1= poly[--k];
   sp= bilion[k];
   while (--i!=0 && *lp1==0 && *(sp++)==0) ++lp1;
  }
  if (i!=0){
   if ((out= open("lmknot.out",1))<0) out=open("temp.out",1);
   lseek (out,(long) 0,2);
   write (out,nbuf,strlen(nbuf));
   write (out,"\n",1);     /* write out polynomial for knot just completed */
   if (count[2]<0) write (out,"program error: knot became inconsistent\n",40);
   if (count[1]!=0) write (out,"coefficient overflow error: output BAD\n",39);
   len=m= -1;
   while (m++!=k){
    if (lowx== (i=0)) write (out,"[",1);
    n= XCNT-m-1;
    j= n*2;
    while (j!=n && poly[m][j]==0 && bilion[m][j]==0) --j;
    while (i!=n && poly[m][i]==0 && bilion[m][i]==0) ++i;
    if (len==0 || lowx>=0 || i!=j || poly[m][i]!=0 || bilion[m][i]!=0){
     while (i<=j){
      if (i==n) write (out,"[",1);
      h=0;
      lngi= poly[m][i];
      if (lngi>=cmpval || lngi<= -cmpval){
       h= lngi/cmpval;
       lngi-= h* cmpval;
      }
      h+= bilion[m][i];
      if (h*lngi <0){   /* bilion and poly are different signs */
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
       write (out,t,ntc((long) h,t));
       len= ntc(lngi,t);
       if (cmpval==10000) len2= 4-len;
       else len2= 9-len;
       write (out,"00000000",len2);
       write (out,t,len);
      }
      else write (out,t,ntc(lngi,t));
      if (i++ ==n) write (out,"]",1);
      if (i<=j) write (out," ",1);
     }
     if (lowx== (len=0)) write (out,"]",1);
     write (out,"\n",1);
    }
    ++lowx;
   }
  }
  write (out,"\n",1);
  close (out);
 }
 i= XCNT;
 while (i!=0) count[--i]=0;
 c=s;
 *nbuf=0;
 if (*c!='1' || (*(c+1)!='+' && *(c+1)!='-')){
  p= nbuf;
  while (*c!='\n' && *c!= EOFCHR) *(p++)= *(c++);
  *p= 0;
  if (*c != EOFCHR && stats>0){
   lseek (stats,(long) 0,0);
   write (stats,nbuf,strlen(nbuf));
   write (stats,"\n",1);
  }
 }
 if (*c=='\n') ++c;
 p= *crsbuf;
 while (*c!='\n' && *c!=EOFCHR){
  if (numcrs==maxcrs){
   write (1,"too many crossings in knot\n",27);
   goto NEWFIL;
  }
  while (*c>='0' && *c<='9') ++c;
  if (*c=='+') sign[numcrs]=6;      /* sign[] says what crossings are + or - */
  else if (*c=='-') sign[numcrs]=2;
  else if (*c==EOFCHR) --c;
  else {
   write (1,"the format of this knot is unreadable, skipping file\n",53);
   close (in);
   goto NEWFIL;
  }
  ++c;
  j=4;
  while (j--!= (i=0) && *c!=EOFCHR){
   while (*c>='0' && *c<='9'){
    i*=10;
    i+= *(c++)-0x30;
   }
   if (i==0){
    write (1,"the format of this knot is unreadable, skipping file\n",53);
    close (in);
    goto NEWFIL;
   }
   *(p++)= i-1;
   if (*c!=EOFCHR){
    if (*c<'a' || *c>'d'){
     write (1,"the format of this knot is unreadable, skipping file\n",53);
     close (in);
     goto NEWFIL;
    }
    *(p++)= (*(c++)-'a')*2;
   }
  }
  ++numcrs;
  if (*c=='\n') ++c;
 }
 if (*c==EOFCHR){
  c=cbuf;
  len=10240;
  while (*s!=EOFCHR){
   *(c++)= *(s++);
   --len;
  }
  if (in!=0 && (i=read(in,c,len))!=0){
   c[i]=EOFCHR;
   c=cbuf;
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
  if (numcrs<2) goto NEWFIL;
 }
 while (*c=='\n') ++c;
 s=c;
 lp1= *poly;
 sp= *bilion;
 lngi= XCNT*XCNT;
 while (--lngi!=0) *(lp1++)= *(sp++)= 0;
 notbeg=numlps=poslnk=neglnk= xpow= ypow= 0;
 if (conchk()!=0) exit(0);
 *count=1;
STEP1:
 skflag=0;
RESTRT:
 if (restrt!=0){
  if ((out=creat("lmpoly.restrt",0644))!= -1){
   write (out,nbuf,82);
   i= strlen(*argv)-11;
   if (i<0) i=0;
   write (out,*argv+i,12);
   *cbuf= numcrs;
   if (ypow<0){
    ypow= -ypow;
    cbuf[6]=1;
   }
   else cbuf[6]=0;
   if ((xpow&512) !=0) cbuf[6]|= 2;
   cbuf[1]= ypow;
   cbuf[2]= xpow;
   cbuf[3]= numlps;
   cbuf[4]= poslnk;
   cbuf[5]= neglnk;
   write (out,cbuf,7);
   i= -1;
   while (++i!=numcrs) crsbuf[i][1]|= sign[i]<<4;
   write (out,*crsbuf,numcrs*8);
   p= (unsigned char *) &lowx;
   write (out,p,2);
   times(&hi);
   kstrt-= hi.tms_utime;
   p= (unsigned char *) &kstrt;
   write (out,p,4);
   p= (unsigned char *) &notbeg;
   write (out,p,4);
   write (out,buf, (int) notbeg);
   write (out,(char *) count,XCNT*4);
   write (out,(char *) plybuf,4*XCNT*XCNT);
   write (out,(char *) bilbuf,2*XCNT*XCNT);
   close (out);
   close (stats);
   exit(0);
  }
  write (1,"couldn't create restart file 'lmpoly.restrt'\n",45);
  restrt=0;
 }
/* FIRST -- remove figure 8 loops, 2 cases  "#-#d#c#b#a"  "#+#b#a#d#c" */
/* SECOND -- remove monogons (must recheck previous crossings) */
/* 4 cases  "#+#b#a????"  "#+????#d#c"  "#-#d????#a"  "#-??#c#b??" */
/* THIRD -- remove bigons -- make sure loops don't get forgotten. */
/* 2 cases "k?mcmb???? m???kbka??"  "k?mc????md  m?????kakd" */
/* bigon yanker MUST be completely thorough or loops can get lost later */
 *dspair=g= i= numcrs;
 c=p= crsbuf[i];
 k=2;
 while (i--!=0){  /* loop through all crossings */
  p-=8;
  if (*p==i){
   if (p[4]==i) ++numlps;
   else mrecon (p+4,p+sign[i]);
   squish (i);
   goto STEP1;
  }
  else if (p[4]==i){
   mrecon (p,p+(sign[i]^4));
   squish (i);
   goto STEP1;
  }
   /* this bigon test checks that ka points to some mc, */
    /* and kd to md, or kb to mb */
  if (p[1]==4 && ((*p==p[6] && p[7]==6) || (*p==p[2] && p[3]==2))){
   if ((j= *p) ==p[2]) k=6;
   c= crsbuf[j];
   if (p[4]==j) ++numlps;
   if (p[2]==p[6]) ++numlps;
   mrecon (p+4,c);
   mrecon (p+k,c+k);
   sqush2 (i,j);
   goto STEP1;
  }
  bstlst[i]= 0;
 }
 if (numcrs<6) g=0;
/* find very good and ok triples */
 while (g--!= (lngi=0)){
  c-=8;
  p= crsbuf[*c]+ ((c[1]^2)&2);
  q= crsbuf[c[4]]+ ((c[5]^2)&2);
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
    ++numlps; /* knot's a distant union with a 2-link */
    if (sign[g]==2) ++neglnk;
    else ++poslnk;
    sqush2 (g,i);
    goto STEP1;
   }
   j=h;
   c+=h;
   if (sign[g]!=sign[i]) ++numlps; /* loose ring */
   else {
    if (sign[g]==2) ++neglnk;
    else ++poslnk;
    j= h^2;
   }
   p= crsbuf[i]+j;  /* only do one mrecon unless I HAVE to do 2 */
   q= c+4;
   if (*c==i){
    if ((c[1]&4) ==0) c= p+4;
    else c=p;
   }
   else if (*q==i){
    if ((q[1]&4) ==0) q= p+4;
    else q=p;
   }
   else mrecon (p,p+4);
   mrecon (c,q);
   sqush2 (g,i);
   goto STEP1;
  }
  if (lngi!=0){
   k= lngi;
   *dspair=1;
   if ((k&3)!=0){
    if ((c[3]&c[7]&2)!=0){
     triple((short) ((k&2)<<1),g,c);
     goto STEP1;
    }
    if (((c[3]^c[7])&2)!=0){
     if ((c[((k&2)<<1)|1]&2)==0){
      if (sign[c[(k&2)<<1]]!=sign[g]) bstlst[i]= bstlst[j]= -19;
      else bstlst[i]= bstlst[j]= 8;
     }
     else {
      n=8;
      if (sign[c[(k&2)<<1]]==sign[g]) n= -8;
      sp= bstlst+i;
      if ((c[3]&2)!=0) sp= bstlst+j;
      if (*sp==0 || (n<0 && n<*sp)) *sp= n;
     }
    }
    else {
     if (bstlst[g]>=0) bstlst[g]+= 8;
     else bstlst[g]-= 8;
    }
   }
   if (k>3){
    if (((c[1]|c[5])&2)==0){
     triple((short) ((k&4)|2),g,c);
     goto STEP1;
    }
    if (((c[1]^c[5])&2)!=0){
     if ((c[(k&4)|3]&2)!=0){
      if (sign[c[(k&4)|2]]!=sign[g]) bstlst[*c]= bstlst[c[4]]= -19;
      else bstlst[*c]= bstlst[c[4]]= 8;
     }
     else {
      n=8;
      if (sign[c[(k&4)|2]]==sign[g]) n= -8;
      sp= bstlst + c[4];
      if ((c[1]&2)!=0) sp= bstlst+ *c;
      if (*sp==0 || (n<0 && n<*sp)) *sp= n;
     }
    }
    else {
     if (bstlst[g]>=0) bstlst[g]+= 8;
     else bstlst[g]-= 8;
    }
   }
  }
 }
 /* circuit remover */
 j= numcrs-1;
 if (j<3) j= -1;
 else if (j<11){    /* <12 will blow up */
  p= *crsbuf+1;
  i=numcrs;
  while (i--!=0 && *p!=4) p+=8; /* only non-alternating (a to c) knots */
  if (i< (n=0)) j= -1;
  else p= crsbuf[j]+5;   /* remove twists */
  while (j>=0){
   n^=4;
   k= *(p--)^4;
   if ((k&2) == 0){
    c= crsbuf[*p]+k;
    i= *c;
    k= c[1]^4;
    if ((k&2) != 0 && crsbuf[i][k]==j){
     untwst ((short) *p,i,j,n);
     squish (j);
     goto STEP1;
    }
   }
   if (n==0) --j;
   p-=3;
  }
 }
 g= donlnk[XCNT]= m= 0;
 while ((n=j) >= (*tt=tt[2]=h= 0)){
  k=g;
  do {
   p= crsbuf[j]+k;
   j= *p;
   k= p[1]^4;
   i= h;
   p=clist;
   while (i--!=0 && *(p++) != j);
   if (i>=0){
    j= n;              /* circuit crosses itself */
    goto NXTCIR; /* skip circuit, wait for a sub-piece that doesn't */
   }
   else {
    ++tt[k&2];
    *p= j;
    t[h++]= (k&2)^2;
   }
  } while (n!=j);
  if (k != g) --tt[k&2];
  else k= -1;
  if (*tt==0 || tt[2]==0){
   rmcir(h,j,k,t,0);
   goto STEP1;
  }
  if (tstcir(h,bstlst,dspair,&skflag) !=0) goto RESTRT;
  if (h==3){
   stc[m++]=j; /* store vertex of 3 crossing circuit for bstlst later */
   stc[m++]=g;
   stc[m++]=k^4;
  }
NXTCIR:
  if (g==0) g= sign[j];
  else {
   --j;
   g=0;
  }
 }
/* if knot is tiny, store polynomial coefficients and restart a stored knot */
 if (numcrs<6){
  lngi= 1;
  if (numcrs==0) --numlps;
  delpow();
  if ((xpow&512)!=0){
   lngi= -1;
   xpow&=511;
  }
  if (numcrs==4){
   c= *crsbuf;
   p= c+4;
   if (*p==*c || c[2]==p[2]){
    g= -3;
    if (*p==p[2] && *c==c[2]) g= -4;
   }
   else {
    g= 0;  /* 0 is flag for knot */
    i= *sign+sign[1]+sign[2]+sign[3];
    if (i!=16) g=2;     /* flag for 2 links   */
   }
  }
  else if (numcrs==5){
   n=10;
   i=j= 0; /* step through knot, if don't get to end, it's a link */
   do {
    p= crsbuf[i]+j;
    i= *p;
    j= p[1]^4;
    --n;
   } while (i!=0 || j!=0);
   if (n>0){
    i=h= 0; /* 5 crossing 2 link  or twisted 3 link  or tref and 2 links */
    g=m= 1;      /* or tref + link */
    k= *sign+sign[1]+sign[2]+sign[3]+sign[4];
    if (k>10 && k<30){
     if (k>20){
      while (sign[i]==6) ++i;
     }
     else {
      while (sign[i]==2) ++i;
     }
     if (k==14 || k==26) goto FOURX;
     c= crsbuf[i];
     p= c+4;
     if (*p==*c || p[2]==c[2]){
      if (*p==p[2] && *c==c[2]) g=3;
      else {
       h=4;
       if (sign[i]==6) m= -1;
      }
     }
     else {
      g=2;
      j=0;
      if (sign[*c]!=sign[i]) j=4;
      k= c[j];
      if (crsbuf[k][j]!=i){
       if (sign[i]==6) m= -1;
      }
      else if (sign[i]==2) m= -1;
     }
    }
    if (g==1){
     if (h==0 && sign[i]==2) m= -1;
     i=j=k= 5;
     c= t;
     while (j!=0) c[--j]=0;
     while (k--!=0){
      j= sign[k]; /* trefs + circles have fake bigons at all crossings */
      p= crsbuf[k]; /* so find all fake bigons */
      if (p[1]==(j^4) && (*p==p[2] || *p==p[6])) c[k]=c[*p]= 1;
      if (p[5]==j && p[4]==p[j]) c[k]=c[p[4]]= 1;
     }
     while (--i>=0 && c[i]!=0);
     if (i>=0){
FOURX: /* found twisted 3 link -- "squish" twist before making polynomial */
      while (++i!=6) sign[i-1]= sign[i];
      numcrs=4;
      g= -3;  /* flag for 4 crossing 3 links */
     }
     else {
      p= *crsbuf;
      if (n==4 || (n==8 && *p==p[4] && p[2]==p[6])) g=3;
     }
    }
   }
   else g= -1; /* it's a knot! */
  }
  i= (*sign)/2 -2;
  k= xpow;
  n=ypow;
  p= *crsbuf;
  if (numcrs==5){
   if (g<0){ /* knots */
    if (i>0){
     if (*p==p[6] && p[2]==p[4] && p[10]==p[12]) g= -2;
    }
    else if (*p==p[2] && p[4]==p[6] && p[8]==p[10]) g= -2;
    if (g== -1){ /* knot five-two */
     addin(lngi,(short)(n-i*6),k,0,0);      /* the flags at the end of addin  */
     addin(lngi,(short)(n-i*4),k,-1,2);     /* let you do two addins at once, */
     addin((long)-lngi,(short)(n-i*2),k,-1,2);    /* for addins with equal ypows.   */
    }
    else { /* knot five-one */
     k+=2;
     addin((long)(-4*lngi),(short)(n-i*4),k,0,0);
     addin((long)-lngi,(short)(n-i*6),k,-2,-2);
     addin(lngi,(short)(n-i*4),(short)(k+2),3,-4);
    }
   }
   else { /* links */
    if (g==1){ /* trefoil + circle */
     if (h!=0){
      addin(lngi,(short)(n-m*3),++k,-1,-2);
      addin((long)(3*lngi),(short)(n-m),k,-1,-2);
      addin(lngi,(short)(n+m),k,-2,-2);
      addin((long)-lngi,(short)(n-m),(short)(k+2),0,0);
     }
     else {
      addin((long)-lngi,(short)(n-m*7),--k,0,0);
      addin((long)(-3*lngi),(short)(n-m*5),k,0,0);
      k+=2;
      addin((long)(3*lngi),(short)(n-m*3),k,0,0);
      addin((long)(2*lngi),(short)(n-m*5),k,0,0);
      addin((long)-lngi,(short)(n-m*3),(short)(k+2),2,-4);
     }
    }
    else if (g==2){ /* 2 component link five-one */
     addin((long)-lngi,(short)(n+m),--k,-2,2);
     addin((long)-lngi,(short)(n-m),k,-1,2);
     addin(lngi,(short)(n+m*3),(short)(k+2),0,0);
     addin((long)-lngi,(short)(n+m),(short)(k+4),0,0);
    }
    else { /* trefoil and circle + circle */
     i= *sign+sign[1]+sign[2]+sign[3]+sign[4];
     m=1;
     if (i<20){
      m= -1;
      i= 40-i;
     }
     j= k-2;
     if (i!=30){
      addin(lngi,(short)(n-m*4),j,-1,2);
      addin((long)(4*lngi),(short)(n-m*2),j,-1,2);
      addin((long)(5*lngi),n,j,0,0);
      addin((long)-lngi,(short)(n+m*2),k,-2,-2);
      k+=2;
      addin(lngi,n,k,-4,-2);
      addin(lngi,(short)(n-m*2),k,0,0);
     }
     else {
      addin(lngi,(short)(n-m*8),j,0,0);
      addin((long)(5*lngi),(short)(n-m*4),j,0,0);
      addin((long)(2*lngi),(short)(n-m*2),j,0,0);
      addin((long)(-2*lngi),(short)(n-m*6),k,-2,-2);
      k+=2;
      addin(lngi,(short)(n-m*4),k,-5,-2);
      addin(lngi,(short)(n-m*2),k,-3,-2);
     }
    }
   }
  }
  else if (numcrs==4){
   if (g==0){ /* knot 4-1 */
    addin((long)-lngi,(short)(n-2),k,0,0);
    addin((long)-lngi,n,k,-1,2);
    addin((long)-lngi,(short)(n+2),k,0,0);
   }
   else if (g>0){
    i=1;
    if (*crsbuf[*p]==0){
     if (*sign==6) i= -1; /* 2 link 4^2-1 */
     addin((long)-lngi,(short)(n+i*3),--k,-1,2);
     addin((long)-lngi,(short)(n+i*5),k,0,0);
     addin((long)-lngi,(short)(n+i),(short)(k+2),0,0);
    }
    else {
     if (*sign==2) i= -1;
     addin((long)-lngi,(short)(n-i*5),--k,-1,2);
     addin((long)-lngi,(short)(n-i*3),k,-3,2);
     addin((long)-lngi,(short)(n-i*3),(short)(k+4),0,0);
    }
   }
   else if (g== -3){ /* 3 link */
    j= *sign+sign[1]+sign[2]+sign[3];
    i= 8- j/2;
    if (i==0){
     addin((long)-lngi,(short)(n-2),k,-1,-2);
     addin((long)(-2*lngi),n,k,-1,-2);
     addin((long)-lngi,(short)(n+2),k,-1,-2);
     addin(lngi,n,(short)(k+2),0,0);
    }
    else {
     j= k-2;
     addin(lngi,(short)(n-2+i),j,0,0);
     addin((long)(2*lngi),(short)(n+i),j,-1,2);
     addin(lngi,(short)(n+2+i),j,0,0);
     addin(lngi,(short)(n+i/2),(short)(k+2),-2,-2);
    }
   }
   else { /* circle + circle  and  circle + circle */
    j= *sign+sign[1]+sign[2]+sign[3];
    i= 8- j/2;
    if (i==0){
     addin(lngi,(short)(n-3),--k,-1,-2);
     addin((long)(3*lngi),(short)(n-1),k,-1,-2);
     addin((long)(3*lngi),(short)(n+1),k,-1,-2);
     addin(lngi,(short)(n+3),k,-1,-2);
     k+=2;
     addin((long)-lngi,(short)(n-1),k,0,0);
     addin((long)-lngi,(short)(n+1),k,0,0);
    }
    else {
     k-=3;
     addin((long)-lngi,(short)(n-3+i),k,0,0);
     addin((long)(-3*lngi),(short)(n-1+i),k,0,0);
     addin((long)(-3*lngi),(short)(n+1+i),k,0,0);
     addin((long)-lngi,(short)(n+3+i),k,0,0);
     k+=2;
     addin((long)(2*lngi),(short)(n-2+(i*3)/4),k,0,0);
     addin((long)(4*lngi),(short)(n+(i*3)/4),k,0,0);
     addin((long)(2*lngi),(short)(n+2+(i*3)/4),k,0,0);
     k+=2;
     addin((long)-lngi,(short)(n-1+i/2),k,0,0);
     addin((long)-lngi,(short)(n+1+i/2),k,0,0);
    }
   }
  }
  else if (numcrs==3){
   addin((long)-lngi,(short)(n-i*4),k,0,0);
   addin(lngi,(short)(n-i*2),(short)(k+2),-2,-2);
  }
  else if (numcrs!=0){
   addin(lngi,(short)(n-i),--k,-1,2);
   addin(lngi,(short)(n-i*3),k,0,0);
  }
  else addin(lngi,n,k,0,0);
  if (notbeg!=0){
   numcrs= buf[notbeg-1];
   numlps= buf[notbeg-4];
   poslnk= buf[notbeg-5];
   neglnk= buf[notbeg-6];
   xpow= buf[notbeg-3];
   ypow= buf[notbeg-2];
   if (buf[notbeg-7]!=0) ypow= -ypow;
   if (buf[notbeg-8]!=0) xpow|= 512;
   notbeg-= 8+ numcrs*8;
   i=0;
   p= *crsbuf;
   c= buf+notbeg;
   while (i!=numcrs){
    j=4;
    while (j--!=0){
     *(p++)= *(c++);
     *(p++)= *(c++)&6;
    }
    sign[i++]= (*(c-7))>>4;
   }
   if (numcrs>chksiz){
    lp1= *poly;
    sp= *bilion;
    lngi= XCNT*XCNT;
    while (lngi--!=0){
     if (*lp1> cmpval){
      if ((*sp)++ == MAXBIL) count[1]=1;
      *lp1-= cmpval;
     }
     else if (*lp1< -cmpval){
      if ((*sp)-- == -MAXBIL -1) count[1]=1;
      *lp1+= cmpval;
     }
     ++lp1;
     ++sp;
    }
    if (stats>0){
     lseek (stats,(long)(strlen(nbuf)+1),0);
     lngi=0;
     i= XCNT-1;
     while (count[i]==0) --i;
     while (i!=0){
      write (stats,t,ntc((long)(i+1),t));
      write (stats,"     ",5);
      write (stats,t,ntc(count[i],t)+1);
      lngi+= count[i--];
     }
     write (stats,"\ntotal: ",8);
     write (stats,t,ntc(lngi,t)+1);
     write (stats,"\nbiggest stacked knots:\n",24);
     i=12;
     sp= bstlst;
     while (i--!=0) sp[i]=0;
     lngi=notbeg-1;
     while (lngi>0){
      k= buf[lngi];
      sp[++i]= k;
      lngi-= k*8 + 8;
      if (i==10) i= -1;
     }
     k= 0;
     while (*sp!=0){
      j= ntc((long) *(sp++),t);
      k+=j+1;
      write (stats,t,j);
      write (stats," ",1);
     }
     write (stats,"                                            ",44-k);
     write (stats,"\n",1);
    }
   }
  }
  else numcrs=0;
 }
 else {
/* locate fake bigons, 0 branch bigons, and 3 cross circuits for bstlst */
  k= numcrs;
  p= crsbuf[k]; /* only set k, other crossings will be set in their loops   */
  c= sign+k;    /* If a tangle is added to the knot with twists, I can find */
  while (k--!=0){ /* it sometimes and remove the twists! */
   n= bstlst[k];
   p-=8;
   j= *--c;
   i= *p;
   g= p[4];
   if (g==p[j]){   /* good fake bigon behind me */
    if (p[5]!=j){
     c= crsbuf[g];  /* wrong! twisted tangle ahead! remove twists! */
     mrecon (p+(j^4),c+ (p[j+1]^4));
     mrecon (p,c+ (p[5]^4));
     sqush2 (k,g);
     goto STEP1;
    }
    if (n>0) n= -4-n;
    else n+= -8;
   }
   if (p[j]==crsbuf[g][sign[g]^j^4] && (p[5]&2)!=0 && nopro!=numcrs){
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
     q= crsbuf[i];  /* wrong! twisted tangle behind! remove twists! */
     mrecon (p+4,q+ (p[1]^4));
     mrecon (p+ *c,q+ (p[j+1]^4));
     sqush2 (k,i);
     goto STEP1;
    }
    if (n>0) n= -4-n;
    else n+= -8;
   }
   else if (p[j^4]==i){  /* bad fake bigon ahead */
    if (n<0) n-= 8;
    else if (n<8) n= 8;
    if (bstlst[i]<0) bstlst[i]-=8;     /* mutualize */
    else if (bstlst[i]<8) bstlst[i]=8;
   }
   if (p[j]==crsbuf[i][sign[i]^j^4] && (p[1]&2)!=0 && nopro!=numcrs){
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
   bstlst[k]=n;
  }
  c= stc;
  i= -8;
  if (nopro==numcrs) m=0;
  while (m>0){
   m-=3;
   n= *(c++);
   g= *(c++);
   p= crsbuf[n];  /* busting vertex of 3 cross circuit makes a removeable */
   sp= bstlst+n;  /* + or - link on the 0 branch -- free 2 reduction */
   k= p[*(c++)];
   j= p[g];
   if ((p[g|1]&2) != 0) p= crsbuf[j];
   else p= crsbuf[j]+sign[j];
   if (*p==k) --i;
   if (*sp<0) *sp+= i-4;
   else if (*sp==0) *sp= 4-i;
   else *sp= i- *sp;
  }
  i= nopro= numcrs;
  k=n= 0;
  sp= bstlst+numcrs;
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
  ypow += bust(j,xpow,ypow);
  xpow^=512;
  ++count[numcrs-1];
 }
 if (numcrs!=0) goto STEP1;
 goto NEWNOT;
}

rmcir(top,vnum,vbrnch,ovundr,flag)
short top, vnum, vbrnch;
unsigned char *ovundr;
int flag;
{
 register short i, j, k, h, m, n;
 register unsigned char *s, *c, *p;
 register long *lp1, *lp2;
 short g;
 if ((m=vbrnch) >=0){
  i= vnum;
  s= crsbuf[i];
  k= s[m];
  g= s[m|1];
  if (m!=0) c= s+4;
  else c= s+(sign[i]^4);
  p= crsbuf[k]+g;
  h= *p= *c;
  n= p[1]= c[1];
  p= crsbuf[h]+n;
  *p= k;
  p[1]= g;
  s[6]= s[4]= s[2]= *s= i;
 }
 else ++numlps;
 s= clist;
 j= -1;
 while (top-1 > (i=k= ++j)){
  while (++i<top){
   if (s[i] < s[k]) k=i;
  }
  if (k!=j){    /* order crossings in clist */
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
  k= crsbuf[i][g];
  j= crsbuf[i][g|1];
  h= crsbuf[k][j]= crsbuf[i][g|4];
  n= crsbuf[k][j|1]= crsbuf[i][g|5];
  crsbuf[h][n]= k;
  crsbuf[h][n|1]= j;
 }
 if (flag!=0) return(0);
      /* squish crossings out of the list */
 i= *s;
 lp1= (long *) crsbuf[i];
 lp2=lp1+suplng;
 p= sign+i;
 c=p+1;
 j= 0;
 k=1;
 while (++i!=numcrs){
  if (s[k]== i){
   lp2+= suplng;
   ++c;
   if (++k==top) k=0;
  }
  else {
   *(lp1++)= *(lp2++);
   if (suplng==2) *(lp1++)= *(lp2++);
   *(p++)= *(c++);
  }
 }
     /* renumber branch pointers (subtract 1 if the crossing was after) */
 numcrs-= top;
 k= numcrs;
 c= crsbuf[k];
 while (k--!=0){
  j=5;
  while (--j!=0){
   c-=2;
   i=top;
   s= clist+i;
   while (i--!=0){
    if (*c> *--s) --*c;
   }
  }
 }
 return(0);
}

ntc(i,buf)
long i;
char *buf;
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

delpow()
{
 register long *lp1, *lp2, *lsp, *osp;
 register short i, j, k, m, s, mpos;
 short numlnk, mpow, u, v, x, y, z, order[2], prslen[2], pstlen[2];
 lp1=b;
 if (numlps<4){
  if (numlps==1){
   *(lp1++)= -1;
   *lp1= -1;
  }
  else if (numlps==3){
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
  b[1]= 4;
  b[2]= 6;
  k= b[3]= 4;
  b[4]= 1;
  while (k++!=numlps){
   j=k;
   b[j]=1;
   while (--j!=0) b[j]+= b[j-1];
  }
  if ((numlps&1)!=0) while (--k>=0) b[k]= -b[k];
 }
 *prslen= *pstlen= prslen[1]= pstlen[1]= mpos= mpow= 0;
 osp= b+ numlps+1;
 if (poslnk>neglnk){
  order[1]= poslnk;
  *order= neglnk;
 }
 else {
  order[1]= neglnk;
  *order= poslnk;
  mpos=2;
 }
 z=2;
 while (z--!=0){
  numlnk=order[z];
  mpos^=2;
  while (numlnk--!=0){
   m= ++mpow;
   s=j=k= numlps+1;
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

addin(kcoeff,ypow,xpow,mult,reach)
short ypow, xpow;
long kcoeff;
int mult, reach;
{
 register long *lp1, *lp2, addon, *p;
 register short i, j, k, numlnk, loops;
 i= poslnk;
 k= neglnk;
 numlnk= i+k;
 loops= numlps+ numlnk;
 ypow+= XCNT-xpow-1+lowx+ 2*(k-i);
 xpow-=lowx+loops;
 p=b;
 while (numlnk-->=0){
  if (reach!=0) lp2= poly[xpow+reach]+ypow-reach;
  lp1= poly[xpow]+ ypow;
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

untwst(l,n,i,ex)
short l, n, i, ex;
{
 register unsigned char *pl, *pn;
 register short g, h, m, k, d, e;
 short j;
 j= sign[i]; /* flip 3 crossing circuits */
 k= j^ex;
 ex^=4;
 m= j^ex;
 if (sign[l]-ex==2){
  g= j^2;
  h= g^4;
  j= 6;
 }
 else {
  h= j^2;
  g= h^4;
  j= 2;
 }
 pl= crsbuf[l];
 pn= crsbuf[n];
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
 crsbuf[pl[m]][pl[m|1]]= l;
 crsbuf[pl[m]][pl[m|1]|1]= m;
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
 crsbuf[pn[h]][m]= n;
 crsbuf[pn[h]][m|1]= h;
 h= pn[g|1];
 crsbuf[pn[g]][h]= n;
 crsbuf[pn[g]][h|1]= g;
 d= pl[k];
 e= pl[k|1];
 crsbuf[d][e]= l;
 crsbuf[d][e|1]= k;
 pl[ex]= g= crsbuf[i][ex];
 pl[ex|1]= h= crsbuf[i][ex|1];
 crsbuf[g][h]= l;
 crsbuf[g][h|1]= ex;
 pn[j]= g= crsbuf[i][k];
 pn[j|1]= h= crsbuf[i][k|1];
 crsbuf[g][h]= n;
 crsbuf[g][h|1]= j;
}

conchk()
{
 register short i, j;
 int vnum, vbrnch, top;
 char str[20], cbuf[XCNT*2];
 vnum=vbrnch= 0;
 lowx= 1;
 i=top= numcrs*2;
 while (i>=0) cbuf[i--]=0;
 i=0;
 while (i!=top){
  --lowx;
  while (cbuf[i]==0){
   cbuf[i]=1;
   i=vnum;
   j=vbrnch;
   vnum= crsbuf[i][j];
   vbrnch= crsbuf[i][j|1];
   if (vnum>=numcrs){
    write (1,"error at crossing #",19);
    write (1,str,ntc((long)(i+1),str));
    write (1,", knot can't have a crossing #",30);
    write (1,str,ntc((long)(vnum+1),str)+1);
    return (1);
   }
   if (crsbuf[vnum][vbrnch]!=i || (crsbuf[vnum][vbrnch|1])!=j){
    write (1,"the possible connection from branch ",36);
    write (1,str,ntc((long)(i+1),str));
    *str= 'a'+(j>>1);
    write (1,str,1);
    write (1," to branch ",11);
    write (1,str,ntc((long)(vnum+1),str));
    *str= 'a'+(vbrnch>>1);
    write (1,str,1);
    write (1," is inconsistent\n",17);
    return (1);
   }
   if (vbrnch==0){
    write (1,"a-c part of crossing #",22);
    write (1,str,ntc((long)(vnum+1),str));
    write (1," is oriented backwards\n",23);
    write (1,"with respect to crossing #",26);
    write (1,str,ntc((long)(i+1),str)+1);
    return (1);
   }
   else if (vbrnch== (sign[vnum] & 6)){
    write (1,"b-d part of crossing #",22);
    write (1,str,ntc((long)(vnum+1),str));
    write (1," is oriented backwards\n",23);
    write (1,"with respect to crossing #",26);
    write (1,str,ntc((long)(i+1),str)+1);
    return (1);
   }
   i= (vnum<<1) + ((vbrnch>>1)&1);
   vbrnch^=4;
  }
  i=0;
  while (cbuf[i]!=0) ++i;
  vnum= i/2;
  if ((i&1)!=0) vbrnch= sign[vnum] & 6;
  else vbrnch=0;
 }
 return(0);
}

pause()
{
 restrt=1;
 signal(SIGTERM,pause);
}

bust(tobust,xpow,ypow)
short tobust, xpow, ypow;
{
 register unsigned char *fastp, *c;
 register long *lp1, *lp2;
 register short i, j, m, k, l, n, *sp;
 unsigned char *inbuf[XCNT], tsign[XCNT], *bstcrs;
 int dir;
 long lngi;
 bstcrs= crsbuf[tobust];
 fastp=sign;
 c=tsign;
 lp2= (long *) *crsbuf;
 lp1= (long *) (buf+notbeg);
 j=0;
 while (j!=numcrs){
  if (j==tobust) ++fastp;
  *(c++)= *(fastp++);
  inbuf[j++]= (unsigned char *) lp1;
  if (suplng==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
 }
 fastp= bstcrs;
 l= sign[tobust]/2;
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
 k= numcrs-1;
     /* squish crossing out of the list (tsign has already been squished) */
 i= k-tobust;
 lp1= (long *) inbuf[tobust];
 lp2=lp1+suplng;
 while (i--!=0){
  if (suplng==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
 }
     /* renumber branch pointers (subtract 1 if the crossing was after) */
 notbeg+= numcrs*8;
 fastp= *inbuf;
 while (++i!=k){
  fastp[1]|= tsign[i]<<4;
  j=5;
  while (--j!=0){
   if (*fastp>=tobust) --*fastp;
   fastp+=2;
  }
 }
 fastp= buf+notbeg-7;
 *fastp= *(fastp-1)=0;
 if (ypow+dir<0){
  *fastp=1;
  buf[notbeg-2]= -ypow-dir;
 }
 else buf[notbeg-2]= ypow+dir;
 if ((xpow&512) ==0) buf[notbeg-8]=1;
 *++fastp= neglnk;
 *++fastp= poslnk;
 *++fastp= numlps;
 *++fastp= xpow+1;
 buf[notbeg-1]= i;
 fastp= bstcrs;
 sp= (short *) fastp;
 j= *sp;
 *sp= sp[l];
 sp[l]= sp[2];
 sp[2]= sp[l^2];
 sp[l^2]= j;
 sign[tobust]^=4;
 crsbuf[*fastp][fastp[1]|1]=0;
 crsbuf[fastp[2]][fastp[3]|1]=2;
 crsbuf[fastp[4]][fastp[5]|1]=4;
 crsbuf[fastp[6]][fastp[7]|1]=6;
 return (dir*2);
}

mrecon(c,p)
unsigned char *c, *p;
{
 register short j, k, m, n;
 register unsigned char *ccb;
 k= *c;
 j= c[1];
 ccb= crsbuf[k]+j;
 n= *ccb= *p;
 m= ccb[1]= p[1];
 ccb= crsbuf[n]+m;
 *ccb= k;
 ccb[1]= j;
}

squish(sqush)
short sqush;
{
 register short i, j, tosqsh;
 register unsigned char *p, *c;
 register long *lp1, *lp2;
       /* squish crossing out of the list */
 tosqsh= sqush;
 j= --numcrs-tosqsh;
 lp1= (long *) crsbuf[tosqsh];
 lp2=lp1+suplng;
 p= sign+tosqsh;
 c=p+1;
 while (j--!=0){
  if (suplng==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
  *(p++)= *(c++);
 }
    /* renumber branch pointers (subtract 1 if the crossing was after) */
 i= numcrs;
 p= crsbuf[i];
 while (i--!=0){
  j=5;
  while (--j!=0){
   p-=2;
   if (*p>=tosqsh) --*p;
  }
 }
}

sqush2(n,g)
short n, g;
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
 lp1= (long *) crsbuf[m];
 lp2=lp1+suplng;
 p=sign+m;
 c=p+1;
 j= i-m;
 while (--j!=0){
  if (suplng==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
  *(p++)= *(c++);
 }
 lp2+=suplng;
 ++c;
 j= numcrs-i;
 numcrs-=2;
 while (--j!=0){
  if (suplng==2) *(lp1++)= *(lp2++);
  *(lp1++)= *(lp2++);
  *(p++)= *(c++);
 }
    /* renumber branch pointers (subtract 1 if the crossing was after) */
 p= crsbuf[numcrs];
 n= numcrs;
 while (n--!=0){
  j=5;
  while (--j!=0){
   p-=2;
   if (*p>=i) --*p;
   if (*p>=m) --*p;
  }
 }
}

triple(m,i,p)
short m, i;
unsigned char *p;
{
 register unsigned char *crsi, *crsk;
 register short j, d, e, n, h, k;
 crsi= p;
 n= m^4;
 k= crsi[m];
 crsk= crsbuf[k];
 j= crsi[m|1];
 d= crsk[j]= crsi[n];
 e= crsk[j|1]= crsi[n|1];
 crsbuf[d][e]= k;
 crsbuf[d][e|1]= j;
 j^=4;
 d= crsi[m]= crsk[j];
 e= crsi[m|1]= crsk[j|1];
 crsbuf[d][e]= i;
 crsbuf[d][e|1]= m;
 crsi[n]= k;
 crsi[n|1]= j;
 crsk[j]= i;
 crsk[j|1]= n;
 j^=2;
 n= crsk[j];
 h= crsk[j|1]^4;
 d= crsk[j]= crsbuf[n][h];
 e= crsk[j|1]= crsbuf[n][h|1];
 crsbuf[d][e]= k;
 crsbuf[d][e|1]= j;
 j^=4;
 n= crsk[j];
 h= crsk[j|1]^4;
 d= crsk[j]= crsbuf[n][h];
 e= crsk[j|1]= crsbuf[n][h|1];
 crsbuf[d][e]= k;
 crsbuf[d][e|1]= j;
 j=m^2;
 n= crsi[j];
 h= crsi[j|1]^4;
 d= crsi[j]= crsbuf[n][h];
 e= crsi[j|1]= crsbuf[n][h|1];
 crsbuf[d][e]= i;
 crsbuf[d][e|1]= j;
 j^=4;
 m= crsi[j];
 h= crsi[j|1]^4;
 d= crsi[j]= crsbuf[m][h];
 e= crsi[j|1]= crsbuf[m][h|1];
 crsbuf[d][e]= i;
 crsbuf[d][e|1]= j;
 sqush2(n,m);
}

tstcir(h,bstlst,dspair,skflag)
short h, *bstlst, *dspair, *skflag;
{
 register short d, g, i, j, m, length, *q, *sp;
 register unsigned char *p, *tp, *vp;
 short k, n, lngsum, badcrs, vertex, lngpos, i1, i2, i3, j1, j2, j3;
 sp=q= tt;  /* # of overs and unders were passed in through *tt and [2] */
 lngsum= n= 0;
 if (*tt > tt[2]) n=2;
 else if (*tt == tt[2]) lngsum=1;   /* need to test BOTH overs & unders */
 vp= clist;
 vertex=h;
 while (vertex-->0) *(sp++)= *(vp++);   /* make a working copy of clist */
 i= h-1;
 if ((h&1)!=(length=0)) vertex= clist[i];
 else {
  if (donlnk[XCNT]==0){   /* if my check-off list is now invalid */
   d= numcrs;             /* re-set it all to 0 */
   while (d-->0) donlnk[d]= 0;
   donlnk[XCNT]=1;
  }
  if ((donlnk[clist[i]]&(t[i]+2))!=0) return(0);
 }
 m=j= h/2;
 *sp= -1;
 p=t;
 while (j-->length){
  i=2;
  while (*q==1024){
   ++q;
   ++p;
  }
  d= *q;
  if (*(p++)==0) i+= sign[d];
  g= crsbuf[d][i&6];
  sp= ++q;
  while (*sp>=0 && *sp!=g) ++sp;
  if (*sp<0) j= -10;
  else if (*sp==vertex) length= -1;
  else *sp= 1024;
 }
 if (j== -11){
  d=h;
  sp=q= tt;
  vp= clist;
  while (d-->0) *(sp++)= *(vp++);
  length=0;
  p=t;
  while (m-->length){
   i=6;
   while (*q==1024){
    ++q;
    ++p;
   }
   d= *q;
   if (*(p++)==0) i+= sign[d];
   g= crsbuf[d][i&6];
   sp= ++q;
   while (*sp>=0 && *sp!=g) ++sp;
   if (*sp<0) m= -10;
   else if (*sp==vertex) length= -1;
   else *sp= 1024;
  }
  j-=m;
 }
 if (vertex< (length=lngsum=d=m= 0)){  /* if a link */
  tp=p= t;
  while ((n^*(p++))==0) ++d;  /* find first bad crossing and move it to the */
  vp= clist+d;                /* front so that clist cannot break a string */
  --p;                        /* of good crossings in a link from */
  sp= tt;                     /* clist[h-1] to clist[0] */
  i=m= h-d;
  while (m-->0){
   donlnk[*vp]|= *(p++)+2;  /* also check off all list elements to assure */
   *(sp++)= *(vp++);        /* this link never comes to tstcir again with */
  }                         /* a different beginning crossing */
  m=d;
  vp= clist;
  while (m-->0){
   donlnk[*vp]|= *tp+2;
   *(sp++)= *vp;
   *(vp++)= *(tp++);
  }
  p=t;
  while (i-->0) *(p++)= *(tp++);
  vp= clist;
  while (d-->0) *(p++)= *(vp++);
  d=m= 0;
 }
 else {
  sp= tt;
  vp= clist;
  i= h;
  while (i-->0) *(sp++)= *(vp++);
 }
 p=t;
 g= h|1;         /* find the longest good string & put in "length" */
 t[g]= t[g-1];   /* also, the longest with only 1 bad crossing in "lngsum" */
 t[g-1]= n^2;    /* badcrs is the bad crossing in lngsum */
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
    if (d!= (lngpos=i)) badcrs= tt[i-d-1];
    else badcrs= tt[i];
   }
   m=d;
   d=0;
  }
 }
 i= h/2;
 if (vertex>=0 && length>=i){ /* can this circuit untwist its vertex? */
  if (k==0 && (*t^*p)!=0) return (skinny(1,h,(short) 0));
  else if (k+length+1==h && (t[k]^*p)==0) return (skinny(1,h,i));
 }
 if (j==0 && length>=i){     /* circuit is not skinny, and can be skinnied */
  if (vertex>= (i1=0) && length==i){
   m= k+i;
   *(p-1)= *p;     /* infinite loops happen when the loose strand is shared */
   d= tt[m-1];        /* between 2 intersecting circuits -- skflag prevents */
   while (m>=k){           /* looping, allows me to skinny once per bust!   */
    vp= crsbuf[tt[m]];
    if (t[m]!=0) j=2;
    if (vp[j]==d || vp[j^4]==d) m=i1= -1;  /* possible inf. loop circuit! */
    else j=0;
    m-=i;
    if (k==0) d= tt[h-1];
    else d=tt[k-1];
   }
  }
  if (i1== (j=0)) return (skinny(0,h,k)); /* OK to skinny it! */
  if (*skflag==0){
   *skflag=1;
   return (skinny(0,h,k)); /* OK to do ONE skinny */
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
  if (j==0 && d!=0 && numcrs+ 2*(i-lngsum)>11) d+=2;
  if ((k=lngpos-lngsum)!= (m=0) && d!=0){
   if (t[k]!=0) m= sign[tt[k]];
   if (crsbuf[tt[k]][m] == tt[k-1]) ++d;
   if (k!=1 && (t[k-2]^t[k-1])!=0) ++d;
  }
  if (lngpos+1 != g && d!=0){
   if (t[lngpos]!= (m=0)) m= sign[tt[lngpos]];
   if (crsbuf[tt[lngpos]][m] == tt[lngpos-1]) ++d;
   if (lngpos+2 !=g && (t[lngpos+1]^t[lngpos])!=0) ++d;
  }
  if (d<6 && j==0){    /* min. one point for making skinny circuits */
   if (bstlst[badcrs]==0) bstlst[badcrs]= *dspair= 1;
  }
  else {
   *dspair=1;
   donlnk[badcrs]=0;  /* from badcrs's perspective, link may be better */
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
   tp=t;
   *gapsto=j= 0;
   d= m+1;
   while (--d>0){
    while ((*(tp++)^n) ==0) ++j;
    gapsto[d]= j++;
   }
   while (j<g){
    while ((*(tp++)^n) ==0) ++j;
    length= (j-gapsto[d]-i)*4 +2;
    if (length>0 && (m<*dspair || (m== *dspair && length>=dspair[1]))){
     if (m<*dspair || length>dspair[1]) dspair[3]=0;
     *dspair=m;
     dspair[1]= length;
     length=lngpos= 0;
     if (d!=m) lngpos= d+1;
     j1= j;
     j2= gapsto[lngpos];
     h=m;
     while (h-->0){
      if (++lngpos>m) lngpos=0;
      vertex= tt[j2];
      j3= gapsto[lngpos];
      lngsum= j1-j3;
      p= crsbuf[vertex];
      i1= sign[vertex];
      i2= p[i1];
      i3= p[i1|1];
      i1^=4;
      if (crsbuf[i2][(i3+i1)&6]== p[4]){
       if ((p[5]^i3&2) !=0) lngsum= -1 -lngsum;
       else lngsum= -10 -lngsum;
      }
      i2= p[i1];
      i3= p[i1|1];
      if (crsbuf[i2][(i3+i1)&6]== *p){
       if ((p[1]^i3&2) !=0){
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
    gapsto[d--]= j++;
    if (d<0) d=m;
   }
  }
 }
 return(0);
}

skinny(twist,lencir,lngpos)
short lencir, lngpos;
int twist;
{
 register short a, b, g, n, vnum, dir, *c;
 register unsigned char *p, *q, *cp;
 short i, j, k, m, last, f, v2, *sp;
 vnum= -1;     /* assume circuit is a link */
 dir= 6;       /* and I will be placing loose strand LEFT of fixed strand */
 if ((lencir&1)!= (b=0)){
  vnum= tt[lencir-1];                /* if not a link */
  if ((t[lencir-1]=t[lencir]) ==0) b= sign[vnum];
    /* tstcir wrecks t[lencir-1] on non-links but saves a copy! */
 }
 if (vnum>=0){
  g=vnum;      /* SET DIR PROPERLY!! */
  while (b>=0){        /* if VERTEX is left of fixed strand, set dir RIGHT */
   p=crsbuf[g]+b;  /* there is also a special case if I am untwisting */
   g= *p;
   b= *(p+1)^4;  /* walk out from the vertex -- do I run into the circuit? */
   a=lencir;
   c= tt;
   while (a>0 && *(c++)!=g) --a;
   if (a>0){
    if (*--c!=vnum){
     if (b==0) dir= sign[*c]^4;     /* yes */
     else dir=b;
    }
    else {
     g= *tt;     /* no -- I ran back into the vertex again */
     b=6;              /* trace the first left region of the circuit */
     if (*t==0) b=sign[g]-2;
     n= tt[lencir-2];  /* does it connect to the other side of the circuit? */
     while (g!=vnum && g!=n){
      p= crsbuf[g]+b;
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
 p= clist;       /* hand.  Put moveable crossings in clist & call rmcir */
 c= tt+lngpos;
 while (a-->0) *(p++)= *(c++);
 --numlps;                   /* rmcir adds a numlps for fun */
 rmcir (j,(short) 0,a,t+lngpos,1);
 i=k= lngpos+j;
 if (vnum>=0){
  --lencir;
  if (lngpos==0 && twist==0) twist= -1;
 }
 else if (i==lencir) i=0;
 if (twist> (m=b= 0) && lngpos!=0){
  if (t[lencir]==0) b= sign[vnum];
  last= crsbuf[vnum][b];   /* set last properly if I am untwisting vertex */
  f= crsbuf[vnum][b|1];
 }
 else {
  last= tt[i];  /* set last to automatically hook up one of the bigon ends */
  f= 4;
  if (t[i]==0) f= sign[last]^4;
 }
 if (t[--k]==0) m=4;
 while (j-->0){
  if (i==lencir) i=0;
  g=dir;
  vnum= tt[i];   /* vnum is the crossing on the fixed strand */
  b= sign[vnum];
  if (t[i++]==0) g+=b;
  else b^=4;
  p= crsbuf[last]+f;  /* reconnect the loose end from the last loop to the */
  *(p++)= n= tt[k--];  /* present n */
  q= crsbuf[n];   /* n is the crossing on the loose strand - q points to it */
  sign[n]= a= b^m;   /* this calculation for sign[n] is a 2 step trick */
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
  cp= crsbuf[vnum]+g;
  v2= *cp;             /* v2 is one step out the dir branch of fixed strand */
  *(cp++)= last= n;
  a= *cp;    /* reconnect vnum to n instead (insert n) */
  *cp=b;
  p= q+b;
  *(p++)= vnum;   /* mutual from n to vnum */
  *p=g;
  b^=4;
  p= crsbuf[v2]+a;  /* insert n from v2 side */
  *(p++)= n;
  *p=b;
  p= q+b;     /* mutual connect from n */
  *(p++)= v2;
  *p=a;
 }
 if (twist!=0){       /* reconnect the loose end properly */
  vnum=tt[lencir];
  if (t[lencir]!= (dir=b=0)) dir=b=sign[vnum];
  if (twist>0){
   p= crsbuf[vnum];
   if (lngpos==0){
    n= p[b^4];
    b= p[b^5];
    if (t[lencir]== (dir=0)) dir=sign[vnum];
   }
   else {
    n= tt[--i];
    if (t[i]== (b=0)) b=sign[n];
   }
   a= p[dir];
   g= p[dir|1];
   k= p[dir^5];
   dir= p[dir^4];
   crsbuf[a][g]= dir;
   crsbuf[a][g|1]=k;
   crsbuf[dir][k]=a;
   crsbuf[dir][k|1]=g;
  }
  else twist=0;
 }
 else if (t[--i]== (b=0)) b= sign[tt[i]];
 if (twist==0) n= tt[i];
 crsbuf[last][f]= n;
 crsbuf[last][f|1]= b;
 crsbuf[n][b]= last;
 crsbuf[n][b|1]= f;
 if (twist!=0) squish(vnum); /* if I untwisted vertex, now eliminate it */
 return(1);
}

