      program companion516
ccc 
c     for nontrivial  knot companions
c
c     parameter (repete = 10000)
      parameter (sts = 4)
      parameter (pi = 3.14159265358979)
ccc
c
      real*8
     :       X,             Y,              Z,
     1       rdomang,       diff,           cdiam,
     2       newcurv,       nlcurv,         x1,
     3       y1,            z1,             zj,
     4       t1,            tj,             sig1,
     5       sgg1,          qq1,            qq2,
     6       pp1,           pp2,            raport,
     7       subangle,      xxc,            sdiam,
     8       diam,          elength,        avcross,
     9       totalcross,    tdiam,          avdiam,
     :       crgm,          trgm,           avrgm,
     :       xxy,           xxz,            cx,
     :       cy,            cz,             xxx,
     :       ndiam,         xj,             yj,
     :       radius,        sradius,        rope,
     :       ssss,          rrss,           length
c
      real*8
     :       points(1500,3),  ronecords(3),  movepts(1500,3),
     2       setpts(1500,3),  rrrr,          matricks(3,3),
     3       rotpts(1500,3),  newpts(1500,3),  newvect(1500,3),
     4       newc(1500),                     pts(1500,3),
     5       d_pts(1500,3),   r_pts(1500,3),   r_delta(1500,3),
     6       t(1500),         rotp(1500,3),    rotd(1500,3),
     7       ss(6),         p(1500,6),       sols(2),
     8       solutions(6001,2), coords(3),   coordinates(6001,3),
     :       orso(6001,2),    ordsols(6001,2,2001), sutions(60001,2)
c
      real*8
     :       tvalues(6001),    svalues(6001),  point(3),    
     3       xx(3001),         yy(3001),         zz(3001), 
     4       reverse(30001,2), p_app(30001,6),   tintix(30001),
     5       p_end(30001,6),   inter(30001),     th(30001,6) 
c
      character*1
     :       ilet(30001),      olet(30001),   
     :       cdwr(30001,9,4),  lett(30001,24)    
c
      character*3
     : ch1, ch2, ch4, ch5, ch6 
c  
      character*4
     :       cd(30001,9), minus     
c
      integer  rone,         rtwo,           r,
     :         counter,      tag(4500),       iseed,
     :         r_tag(4500),   stp,            failure,
     :         rotag(4500),   sn,             snn,
     :         ith(30001,4),  inlet(30001),    onlet(30001),
     :         num(30001),    walk,           number,
     :         irstrt,       istef1,         maxcross
c
c     character*80
c    : header
c      Set Maximum Number of Crossings
c
 533   format(3(2x,f21.17))
 534   format(a80)
 535   format(3(2x,f21.17))

      maxcross = 999
c
      open(7,file='zmatrix', status='unknown',form ='formatted')
      open(8,file='zpoints', status='unknown',form ='formatted')
      open(9,file='zkdata', status='unknown',form ='formatted')
c
      write(6,*) 'number of edges [max 101] - integer number'
      read(5,*) number
      write(6,*) 'number = ', number
      ndiam = number 
c
      write(6,*) 'how many cases?'
      read(5,*) repete
      write(6,*) 'number of cases = ', repete
c
      write(6,*) 'restart? - please type yes=1 / no=0 '
      read(5,*) irstrt
      write(6,*) irstrt
 1993 format(i1)
      if(irstrt.eq.1) then
       write(6,*) 'now reading restart file res.pts'
       open(11, file = 'res.pts', status ='old', form = 'formatted')
       do 1831 i = 1,number
       read(11,533) points(i,1), points(i,2), points(i,3)
       write(6,533) points(i,1), points(i,2), points(i,3)
 1831  continue
      else
c 
c - generation of planar points -
c
       subangle = (-2.0 *  pi)/number
c
       do 1995 i = 1, number
       points(i,1) = dcos(i*subangle)
       points(i,2) = dsin(i*subangle)
       points(i,3) = 0.0
 1995  continue
       write(6,*) 'new run - no restart'
      endif
c
c
      write(6,*) 'starting points'     
      do 1997 i = 1, number
      write(6,533)  points(i,1), points(i,2), points(i,3)
 1997 continue
c   
      write(6,*) 'record points in zpoints? - please type yes=1 / no=0 '
      read(5,*) zpts
      write(6,*) ' record points? ', zpts
c
      write(6,*)'enter value for iseed call - integer number  -'
      read(5,*) iseed
      write(6,*)'seed = ',iseed
      istef1 = irand(iseed)
c
      write(6,*)'istef1 = ',istef1
c
      read(5,*) rope
      write(6,*) 'ropelength is ',rope
c
      read(5,*) sradius
      write(6,*) 'perturbation radius is ',sradius
      write(6,*) 'so far so good'
      elength = 0
      do 1999 j = 1, number -1
      elength = 0
      do 1998 i = 1, 3
      xxc = points(j+1,i) - points(j,i)
      elength = elength + xxc*xxc
 1998 continue
      length = length + dsqrt(elength)
 1999 continue
c
      do 2000 i = 1, 3
      xxc = points(1,i) - points(number,i)
      elength = elength + xxc*xxc
 2000 continue
      length = length + dsqrt(elength)
      write(6,*) ' knot length = ',length
c
      radius = sradius
      write(6,*) ' perturbation radius = ',radius
c
c 
      counter = 0
      tdiam = 0
      avdiam = 0
      totalcross = 0
      avcross = 0
      ircont = 1
      intix = 0 
      jrcont = 1
      walk = 1
      trgm = 0
      avrgm = 0
      newcurv = 100
c
      do 1001 kcm = 1, walk
c
      do 1011 ish = 1, repete
      intix = intix + 1
c
 192  continue
c 
c -   set to zero .........
c   
      lp = 1
      lrev = 1
      do 561 i = 1, 3001
      sutions(i,1) = 0.0
      reverse(i,1) = 0.0
      sutions(i,2) = 0.0
      reverse(i,2) = 0.0
c	
      p_app(i,1) = 0.0
      p_app(i,2) = 0.0
      p_app(i,3) = 0.0
      p_app(i,4) = 0.0
      p_app(i,5) = 0.0
      p_app(i,6) = 0.0
c
	
      p_end(i,1) = 0.0
      p_end(i,2) = 0.0
      p_end(i,3) = 0.0
      p_end(i,4) = 0.0
      p_end(i,5) = 0.0
      p_end(i,6) = 0.0
c
      inter(i) = 0.0
      num(i) = 0
      tintix(i) = 0.0
      inlet(i)  = 0
      onlet(i) = 0
      ilet(i) = ' '
      olet(i) = ' '
      ith(i,1) = 0
      th(i,1) = 0.0
      ith(i,2) = 0
      th(i,2) = 0.0
      ith(i,3) = 0
      th(i,3) = 0.0
      ith(i,4) = 0
      th(i,4) = 0.0
  561 continue
      do 567 j = 1, 9 
      do 569 i = 1, 3001 
      cd(i,j) = ' '
      cdwr(i,j,1) = ' '
      cdwr(i,j,2) = ' '
      cdwr(i,j,3) = ' '
      cdwr(i,j,4) = ' '
  569 continue
  567 continue
      do 572 j = 1, 24 
      do 574 i = 1, 3001 
      lett(i,j) = ' '
  574 continue
  572 continue
c
      X = 0.0
      Y = 0.0
      Z = 0.0
      do 1501 i = 1, number
      do 1203 j = 1,3
      movepts(i,j) = 0.0
 1203 continue
 1501 continue
      do 1207 i = 1,3
      ronecords(i) = 0.0
 1207 continue
      do 1209 i = 1, number
      xx(i) = 0.0
      yy(i) = 0.0
      zz(i) = 0.0
      r_tag(i) = 0
      tag(i) =  0
      rotag(i) = 0
      setpts(i,1) = 0.0
      setpts(i,2) = 0.0
      setpts(i,3) = 0.0
 1209 continue
c
c
c
c
      seed = iseed
c
c
c 	Random Vertex Perturbation of Radius = radius
c
      do 1211 i = 1,number
	rrrr = rand(0)
	ssss = rand(0)*2*pi
c     write(6,*) 'rrrr = ',rrrr,'  ssss = ',ssss
        rrss = rand(0)
c     write(6,*) ' rrss = ',rrss
c       rrss = rrss**(1.0/3.0)
c
c     write(6,*) ' rrss = ',rrss
      if(rrrr.lt.0.5) then
         rdomang = (-1)*pi*rrrr
      else
         rdomang = (rrrr - 0.5)*pi
      endif
      point(1) = rrss*radius*dcos(rdomang)
      point(2) = rrss*radius*dsin(rdomang)*dcos(ssss)
      point(3) = rrss*radius*dsin(rdomang)*dsin(ssss)
      do 1213 j = 1,3
      newpts(i,j) = points(i,j) + point(j)
 1213 continue
 1211 continue
c
c	initialization
c
      do 1215 i = 1,3
      do 1217 j = 1, number
      newvect(j,i) = 0.0
      pts(j,i) = 0.0
      d_pts(j,i) = 0.0
      r_pts(j,i) = 0.0
      r_delta(j,i) = 0.0
      rotp(j,i) = 0.0
      rotd(j,i) = 0.0
 1217 continue
 1215 continue
      do 1219 i = 1,number 
      newc(i) = 0.0
 1219 continue
c
c
c
C	we now have the new knot
C
c -   newvect -----
c
      do 31 i = 1,(number-1)
      newvect(i,1) = newpts(i+1,1) - newpts(i,1)
      newvect(i,2) = newpts(i+1,2) - newpts(i,2)
      newvect(i,3) = newpts(i+1,3) - newpts(i,3)
  31  continue
      newvect(number,1) = newpts(1,1) - newpts(number,1)
      newvect(number,2) = newpts(1,2) - newpts(number,2)
      newvect(number,3) = newpts(1,3) - newpts(number,3)
c
      newcurv = 0.0
c
      do 33 i = 1, (number-1)
      newc(i) = dacos( (newvect(i,1)*newvect(i+1,1)
     :          +       newvect(i,2)*newvect(i+1,2)
     :          +       newvect(i,3)*newvect(i+1,3) )/
     :         (dsqrt(  newvect(i,1)*newvect(i,1)
     :          +       newvect(i,2)*newvect(i,2)
     :          +       newvect(i,3)*newvect(i,3)   )*             
     :          dsqrt(  newvect(i+1,1)*newvect(i+1,1)
     :          +       newvect(i+1,2)*newvect(i+1,2)
     :          +       newvect(i+1,3)*newvect(i+1,3) ) ) )
      newcurv = newcurv + newc(i)
 33   continue
      nlcurv = dacos(  (newvect(number,1)*newvect(1,1)
     :          +       newvect(number,2)*newvect(1,2)
     :          +       newvect(number,3)*newvect(1,3) )/
     :         (dsqrt(  newvect(number,1)*newvect(number,1)
     :          +       newvect(number,2)*newvect(number,2)
     :          +       newvect(number,3)*newvect(number,3) )*             
     :          dsqrt(  newvect(1,1)*newvect(1,1) 
     :          +       newvect(1,2)*newvect(1,2)
     :          +       newvect(1,3)*newvect(1,3) ) ) ) 
      newcurv = newcurv + nlcurv
c
c       Calculate  diameter
c     write(6,*) 'calculate diameter'
      cdiam = 0
       do 560 idm = 1, number
        do 562 jdm = idm + 1, number
        diam = 0
         do 563 kdm = 1, 3
         xxc = newpts(idm,kdm) - newpts(jdm,kdm)
         diam = diam + xxc*xxc
  563    continue
       diam = dsqrt(diam)
c      write(6,*) 'diam = ',diam
        if(diam.gt.cdiam) then
        cdiam = diam
        endif
  562   continue
  560  continue
c      write(6,*) ' cdiam = ',cdiam
       tdiam = tdiam + cdiam
c        End of Diameter Calculation
c
c        Calculate Radius of Gryation
       cx = 0
       cy = 0
       cz = 0
        do 570 irm = 1, number
	    cx = cx + newpts(irm,1)
            cy = cy + newpts(irm,2)
            cz = cz + newpts(irm,3)
  570   continue
	 cx = cx/number
         cy = cy/number
         cz = cz/number
	 crgm = 0
       do 580 idm = 1, number
 	   xxx = newpts(idm,1) - cx
 	   xxy = newpts(idm,2) - cy
	   xxz = newpts(idm,3) - cz
	   crgm = crgm + dsqrt(xxx*xxx + xxy*xxy + xxz*xxz)
  580    continue
c      write(6,*) 'diam = ',diam
	 rgm = crgm/number
       trgm = trgm + rgm
       counter = counter + 1
c	   End of Radius of Gyration Calculation
c
c
c533   format(f9.6,x,f9.6,x,f10.6)
c535   format(f9.6,x,f9.6,x,f9.6)
c
c      do 1035 i = 1, number
c      points(i,1) = newpts(i,1)
c      points(i,2) = newpts(i,2)
c      points(i,3) = newpts(i,3)
c1035 continue
c
       if (zpts.eq.1) then
		write(8,3779) counter, rgm, cdiam
 3779 format(' # ',i6,'  rgm =  ',f21.17,'  diam = ',f21.17)
c
         do 373 i = 1, number 
         write(8,533) newpts(i,1), newpts(i,2), newpts(i,3)
  373 continue
      endif
c
c
       if (mod(intix,1000).eq.0) then
	 write(9,*) counter, '"case = "',kcm,'"."',ish,
     :    '" diam = "',cdiam,'" rgm = "',rgm,'"    points = "'
         do 374 i = 1, number
         write(9,533) newpts(i,1), newpts(i,2), newpts(i,3)
  374 continue
      endif
 3374 continue
      lrev = 1
      lp = 1
c
c --- first while -------------------------------------------------------|
c                                                                        |
      if(newcurv.le.pi/10000) go to 999 
c
      do 35 i = 1, number
      pts(i,1) = newpts(i,1)
      pts(i,2) = newpts(i,2)
      pts(i,3) = newpts(i,3)
  35  continue
c
      nn = number
c
c     write(6,*) counter, '"next case = "',kcm,'"."',ish,
c    :              '"pts = "'
c     write(6,*) 'number = ',number
c      do 973 i = 1, number
c      write(6,533) pts(i,1), pts(i,2), pts(i,3)
c973  continue
c     write(6,*) ' d_pts = '
c
      do 37 j=1,(nn-1)
      d_pts(j,1) = pts(j+1,1) - pts(j,1)
      d_pts(j,2) = pts(j+1,2) - pts(j,2)
      d_pts(j,3) = pts(j+1,3) - pts(j,3)
c     write(6,*) d_pts(j,1)
c     write(6,*) d_pts(j,2)
c     write(6,*) d_pts(j,3)
  37  continue
c
      d_pts(nn,1) = pts(1,1) - pts(nn,1)
      d_pts(nn,2) = pts(1,2) - pts(nn,2)
      d_pts(nn,3) = pts(1,3) - pts(nn,3)
c     write(6,*) d_pts(nn,1)
c     write(6,*) d_pts(nn,2)
c     write(6,*) d_pts(nn,3)
 
c
c -tag------------------
      do 36 j=1,nn
      tag(j) = j
  36  continue
c - end tag
c
c - rotate ---------------------
c
      do 739 j=1,(nn-1)
      r_pts(j,1) = pts(j+1,1)
      r_pts(j,2) = pts(j+1,2)
      r_pts(j,3) = pts(j+1,3)
      r_delta(j,1) = d_pts(j+1,1)
      r_delta(j,2) = d_pts(j+1,2)
      r_delta(j,3) = d_pts(j+1,3)
 739  continue
      r_pts(nn,1) = pts(1,1)
      r_pts(nn,2) = pts(1,2)
      r_pts(nn,3) = pts(1,3)
      r_delta(nn,1) = d_pts(1,1)
      r_delta(nn,2) = d_pts(1,2)
      r_delta(nn,3) = d_pts(1,3)
c
c - rotate tag
      do 738 i = 1, (nn-1)
      r_tag(i) = tag(i+1)
c      write(6,*) r_tag(i)
 738  continue
      r_tag(nn) = tag(1)
c
c - end rotate ------------------
c
c ------ second while ------------------------------------------|
c                                                               |
      ind = 0
      stp = nn
c
c - set to zero ---------------------
      do 517 i = 1, number 
      do 519 k = 1, number-1 
      ordsols(i,1,k) = 0.0
      ordsols(i,2,k) = 0.0
  519 continue
  517 continue
      snn  = 0
      sn = 0
c
c - set to zero ---------------------
c
  704 continue
c
      if( ind.ge.stp) go to 703
      ind = ind + 1
c - rotate ---------------------
      do 39 j=1,nn
      rotp(j,1) = r_pts(j,1)
      rotp(j,2) = r_pts(j,2)
      rotp(j,3) = r_pts(j,3)
      rotd(j,1) = r_delta(j,1)
      rotd(j,2) = r_delta(j,2)
      rotd(j,3) = r_delta(j,3)
  39  continue
c
c - rotate tag
      do 38 i = 1, nn
      rotag(i) = r_tag(i)
  38  continue
c
c - end rotate ------------------
c
c - initialize p, solutions, coordinates
c
      npts = 0
      do 7501 i = 1,number
      p(i,1) = 0.0
      p(i,2) = 0.0
      p(i,3) = 0.0
      p(i,4) = 0.0
      p(i,5) = 0.0
      p(i,6) = 0.0
 7501 continue
      do 503 i = 1,number
      solutions(i,1) = 0.0
      solutions(i,2) = 0.0 
      orso(i,1) = 0.0
      orso(i,2) = 0.0
 503  continue
      do 505 i = 1,number
      coordinates(i,1) = 0.0
      coordinates(i,2) = 0.0
      coordinates(i,3) = 0.0
 505  continue
      do 507 i =1,number
      tvalues(i) = 0.0
      svalues(i) = 0.0
      t(i) = 0.0
 507  continue
      t1 = 0.0
      tj = 0.0
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
      zj = 0.0
      qq1 = 0.0
      pp1 = 0.0
      pp2 = 0.0
      do 511 i = 1,6
      ss(i) = 0.0
  511 continue
      do 513 i = 1,2
      sols(i) = 0.0
  513 continue
      do 515 i = 1,3
      coords(i) = 0.0
  515 continue
c
c - end of initialize: p, solutions, coordinates
c
      do 41 k = 3, (nn-1)
      diff = r_delta(k,1)*r_delta(1,2) - r_delta(1,1)*  
     1       r_delta(k,2)
c     write(6,*) 'diff = ',diff
c    
      if(diff.eq.0.0) go to 40 
      t(r_tag(1)) = ((r_pts(k,2) - r_pts(1,2))*r_delta(k,1) +
     :             (r_pts(1,1) - r_pts(k,1))*r_delta(k,2) )/
     :             diff
      t(r_tag(k)) = ((r_pts(k,2) - r_pts(1,2))*r_delta(1,1) +
     :             (r_pts(1,1) - r_pts(k,1))*r_delta(1,2) )/
     :             diff
c
      t1 = t(r_tag(1))
      tj = t(r_tag(k))
      z1 = r_pts(1,3) + r_delta(1,3) * t1
      x1 = r_pts(1,1) + r_delta(1,1) * t1
      y1 = r_pts(1,2) + r_delta(1,2) * t1
      zj = r_pts(k,3) + r_delta(k,3) * tj
      xj = r_pts(k,1) + r_delta(k,1) * tj
      yj = r_pts(k,2) + r_delta(k,2) * tj
c     write(6,*) 't1 = ',t1,' tj = ',tj
c     write(6,*) 'x1 = ',x1,' xj = ',xj
c     write(6,*) 'y1 = ',y1,' yj = ',yj
      sig1 = z1 - zj
c     write(6,*) ' sig1 = ', sig1
c     write(6,*) ' at 801 lp = ',lp
c
      if( sig1.lt.0) then 
        sn = -1
        go to 801
        endif
      if( sig1.eq.0) then 
        sn = 0
        go to 801
        endif
      if( sig1.gt.0) then 
        sn = 1 
        go to 801
        endif
c
 801  continue
c
      qq1 = r_delta(1,1)
      qq2 = r_delta(1,2)
      pp1 = r_delta(k,1)
      pp2 = r_delta(k,2)
c
      sgg1 = qq1*pp2 - qq2*pp1
      if( sgg1.lt.0) then 
        snn = -1
        go to 802
        endif
      if( sgg1.eq.0) then 
        snn = 0
        go to 802
        endif
      if( sgg1.gt.0) then  
        snn = 1
        go to 802
        endif
c
 802  continue
c
       if (0.0.lt.t1.and.1.0.gt.t1) then
c
          if(0.0.lt.tj.and.1.0.gt.tj) then
c
              if(sn.eq.0) go to 40
              failure = 0
c
              npts = npts + 1
c
              ss(1) = 0.0
              ss(2) = snn
              ss(3) = sn
c
              ss(4) = x1
              ss(5) = y1
              ss(6) = t1
c
              do 51 i = 1,6
              p(npts,i) = ss(i)
  51          continue
c
c     write(6,*) 'ind =',ind,' npts =',npts,' p =',(p(npts,i), i =1,6) 
c
c
              sols(1) = t1
              sols(2) = tj
c
              do 53 i =1,2
              solutions(npts,i) = sols(i)
  53          continue
c
           do 63 i = 1, npts
           orso(i,1) = solutions(i,1)
           orso(i,2) = solutions(i,2)
  63       continue
c
              coords(1) = x1
              coords(2) = y1
              coords(3) = z1
c
              do 55 i=1,3
              coordinates(npts,i) = coords(i)
  55          continue
c
          else
              junk = 0
          endif
        else
          junk = 0
        endif
  40  continue
c
      sn = 1
      failure = 2
c
      if(t1.eq.1.0) go to 42 
         failure = 0
c
      if(t1.eq.0.0) go to 42 
         failure = 0
c
      if(tj.eq.1.0) go to 42 
         failure = 0
c
      if(tj.eq.0.0) go to 42 
         failure = 0
c
      if(sn.eq.0) go to 42 
         failure = 0
c
 41   continue
c
 42   continue
c
       if(npts.ge.199) then
         write(6,*) ' ........ attention .....'
         write(6,*) 'npts = ', npts
         write(6,*) ' ........ attention .....'
       endif
c
c      write(6,*) '551 npts = ',npts
c      write(6,*) '551 lp = ',lp
c
              do 553 j = 1, npts
              do 551 i = 1,6
              p_app(lp,i) = p(j,i)
 551    continue
c       write(6,*) 'lp =',lp,' npts=',npts,' p_app =',(p(j,i), i =1,6)
              lp = lp + 1
 553    continue
c      
c       write(6,*) 'control of lp =',lp
         if(npts.gt.0) then
c             write(6,*) 'reduction of lp'
              lp = lp - 1
         endif
        if(npts.gt.1) then      
              do 1553 j = 1, npts      
              inter(j) = p_app(lp-npts+j,6)
 1553         continue
c       write(6,*) ' lp = ',lp
c
c       write(6,*) ' npts = ',npts,'  inter = ',(inter(j), j=1,npts)
c
              call ssort(npts,inter)
c
              do 1555  j = 1, npts 
                 do 1557 k = 1, npts 
                    if(inter(j).eq.p_app(lp-npts+k,6)) then
                       do 1559 m = 1, 6
                          p_end(lp-npts+j,m) = p_app(lp-npts+k,m)
 1559                  continue
                    endif
 1557             continue
 1555          continue
c
        elseif (npts.eq.1) then
              do 1563 i = 1,6
              p_end(lp,i) = p_app(lp,i)
 1563         continue
c       write(6,*) ' lp = ',lp
c
c       write(6,*) ' npts = ',npts,'p_end  = ',(p_end(lp,j), j=1,6)

        else
        endif
      continue
c       
        if(npts.gt.0) then
         lp = lp + 1
        endif
c 
      do 61 i=1,npts
      tvalues(i) = p(i,6)
  61  continue
c
      do 62 i=1,npts
      svalues(i) = tvalues(i)
  62  continue 
      if(npts.gt.1) then
c
      call ssort(npts,svalues)
c
      call ssort(npts, orso)
c
      endif
      do 65 i = 1, npts
         do 67 j = 1,npts
           if (orso(i,1).eq.solutions(j,1))  then
               orso(i,2) = solutions(j,2)
           endif
  67     continue
  65  continue
c
      do 69 k = 1, npts
         ordsols(ind,1,k) = orso(k,1)
         ordsols(ind,2,k) = orso(k,2)
 69   continue
c
      if(ind.lt.stp) then
c
c - rotate ---------------------
c
      do 839 j=1,(nn-1)
      r_pts(j,1) = rotp(j+1,1)
      r_pts(j,2) = rotp(j+1,2)
      r_pts(j,3) = rotp(j+1,3)
      r_delta(j,1) = rotd(j+1,1)
      r_delta(j,2) = rotd(j+1,2)
      r_delta(j,3) = rotd(j+1,3)
 839  continue
      r_pts(nn,1) = rotp(1,1)
      r_pts(nn,2) = rotp(1,2)
      r_pts(nn,3) = rotp(1,3)
      r_delta(nn,1) = rotd(1,1)
      r_delta(nn,2) = rotd(1,2)
      r_delta(nn,3) = rotd(1,3)
c - rotate tag
      do 838 i = 1, (nn-1)
      r_tag(i) = rotag(i+1)
c      write(6,*) r_tag(i)
 838  continue
      r_tag(nn) = rotag(1)
c - end rotate ------------------
      go to 704 
      endif
c
 703  continue
c
      if(failure.ne.0) go to 999
      junk = 0
c
c     means not in Gen Posn
c
      l=1
      do 71 i = 1, number 
         do 73 k = 1, number-1 
            if( ordsols(i,1,k).eq.0.0) then
                 go to 71
             else 
                sutions(l,1) = ordsols(i,1,k)
                sutions(l,2) = ordsols(i,2,k)
                l = l + 1
            endif
 73      continue
 71   continue
      lmax = l - 1
c     write(6,*) ' lmax = ', lmax
c
      if (lmax.ge.2*maxcross) then
       totalcross = totalcross + lmax/2
       go to 999
c        write(6,*) ' ATTENTION '
c        write(6,*) ' lmax greater that 1001 - increase memory - ', lmax
c        write(6,*) ' ATTENTION '
      endif
c
      failure = 1
      if(lmax.eq.0) then
c      cdiam = 0
c      do 1560 idm = 1, number
c       do 1562 jdm = idm + 1, number
c       diam = 0
c        do 3563 kdm = 1, 3
c        xxc = newpts(idm,kdm) - newpts(jdm,kdm)
c        diam = diam + xxc*xxc
c3563    continue
c      diam = dsqrt(diam)/elength
c      write(6,*) 'diam = ',diam
c       if(diam.gt.cdiam) then
c       cdiam = diam
c       endif
c1562   continue
c1560  continue
c      write(6,*) ' cdiam = ',cdiam
c        End of Diameter Calculation
       write(7,3779) counter, rgm, cdiam
       write(7,*) '1+1b1a1d1c'
       write(7,1795)
c      if(cdiam.lt.sdiam) then
c       write(6,*) ' cdiam = ',cdiam
c       do 8036 i = 1, number
c       points(i,1) = newpts(i,1)
c       points(i,2) = newpts(i,2)
c       points(i,3) = newpts(i,3)
c8036   continue
c      else
        if(cdiam.lt.ndiam) then
c        do 8037 i = 1, number
c        points(i,1) = newpts(i,1)
c        points(i,2) = newpts(i,2)
c        points(i,3) = newpts(i,3)
c8037   continue
         ndiam = cdiam
c        write(6,*) ' ndiam = ',ndiam
        else
        go to 999
        endif
c      endif
c
       go to 999
      endif
      failure = 0
      if(lmax.le.4) then
c       cdiam = 0
c       do 2560 idm = 1, number
c       do 2562 jdm = idm + 1, number
c       diam = 0
c        do 2563 kdm = 1, 3
c        xxc = newpts(idm,kdm) - newpts(jdm,kdm)
c        diam = diam + xxc*xxc
c2563    continue
c      diam = dsqrt(diam)/elength
c      write(6,*) 'diam = ',diam
c       if(diam.gt.cdiam) then
c       cdiam = diam
c       endif
c2562   continue
c2560  continue
c      write(6,*) ' cdiam = ',cdiam
c        End of Diameter Calculation
       write(7,3779) counter, rgm, cdiam 
       write(7,*) '1+2c2b1d1c'
       write(7,*) '2-2d1b1a2a'
       write(7,1795)
c      if(cdiam.lt.sdiam) then
c       write(6,*) ' cdiam = ',cdiam
c       do 9036 i = 1, number
c       points(i,1) = newpts(i,1)
c       points(i,2) = newpts(i,2)
c       points(i,3) = newpts(i,3)
c9036   continue
c      else
        if(cdiam.lt.ndiam) then
c        do 9037 i = 1, number
c        points(i,1) = newpts(i,1)
c        points(i,2) = newpts(i,2)
c        points(i,3) = newpts(i,3)
c9037   continue
         ndiam = cdiam
c        write(6,*) ' ndiam = ',ndiam
        else
        go to 999
        endif
c      endif
       go to 999
      endif
c
c
c     means has 0 crossings ------------
c      write(6,*) ' .......... sutions ............. '
c      do 75 l = 1, lmax
c      write(6,711)  sutions(l,1), sutions(l,2) 
c 75   continue 
  711  format(2x,2(f10.6,2x))
c
c
c
      do 77 k = 1, lmax
      reverse(k,1) = sutions(k,2)
      reverse(k,2) = sutions(k,1)
 77   continue
c      write(6,*) '.......... reverse ...............'
c      do 79 l = 1, lmax
c      write(6,711) reverse(l,1), reverse(l,2)
c 79   continue
c
c      write(6,*) ' .......p_end ................'
c      do 81 l  = 1,  lmax
c      write(6,811) ( p_end(l,i), i = 1,6 )
c 81   continue
  811  format(6(2x,f10.6))
c
      do 83 k = 1, lmax
         do 85 j = 1, lmax
            if(sutions(j,1).eq.reverse(k,1)) then
              tintix(lrev) = j 
              lrev =  lrev +  1
c             write(6,*) 'j= ',j,' k= ',k,' lrev= ',lrev
            endif
  85     continue
  83  continue
c
c      write(6,*) ' .........  t(lrev) .....'
c      write(6,*) (tintix(i), i = 1, lmax)
c
      do 91 i = 1, lmax
      th(i,1) = p_end(i,2)
      th(i,2) = p_end(i,3)
      th(i,3) = i
      th(i,4) = tintix(i)
  91  continue
      do 95 i = 1, lmax
      ith(i,1) = idint(th(i,1))
      ith(i,2) = idint(th(i,2))
      ith(i,3) = idint(th(i,3))
      ith(i,4) = idint(th(i,4))
  95  continue
c     write(6,*) ' ... Thistlethwaite code ...'
c     do 97 i = 1, lmax
c     write(6,923) (ith(i,j), j = 1, 4)
c 97  continue
 923  format(4(2x,i6))
c     write(6,*) ' '
c
c -  calculate the standard (Millett/Ewing convention)
c    matrix of the link, starting with the Thistlethwaite-
c    Representation in ith() -
c
      do 101 k =1, lmax
      if(mod(ith(k,3),2).eq.0) then
         num(k) = ith(k,3)/2
      else
         num(k) = ith(k,4)/2
      endif
 101  continue
      do 103 k =1, lmax
      if(ith(k,2).eq.1) then
         ilet(k) = 'c'
         inlet(k) = 5
         olet(k) = 'a'
         onlet(k) = 3
      else
         if(ith(k,1).eq.1) then
            ilet(k) = 'd'
            inlet(k) = 6
            olet(k) = 'b'
            onlet(k) = 4
         else
            ilet(k) = 'b'
            inlet(k) = 4
            olet(k) = 'd'
            onlet(k) = 6
         endif
      endif
 103  continue
c - start of the while5 loop --------------|
c
         write(ch2, '(i3)' ) num(1)
         write(ch4, '(i3)' ) num(lmax)
c
      do 501 i = 1, lmax
c
      if(ith(i,2).eq.1) then
         j = num(i)
c - stef
         write(ch1, '(i3)' ) j
c
       if(i.lt.lmax) then
         write(ch5, '(i3)' ) num(i+1)
       else
         write(ch5, '(i3)' ) num(1)
       endif
c
       if(i.gt.1) then
         write(ch6, '(i3)' ) num(i-1)
       else
         write(ch6, '(i3)' ) num(lmax)
       endif 
c - stef
c
         cd(j,1) = ch1 
         iminus = ith(i,1)*ith(i,2)
         write(minus, '(i3)') iminus
         cd(j,2) = minus 
c
           if(cd(j,2)(1:1).eq.'-') then
              go to 1071
              elseif (cd(j,2)(2:2).eq.'-') then
              go to 1071
              elseif (cd(j,2)(3:3).eq.'-') then
              go to 1071
              elseif (cd(j,2)(4:4).eq.'-') then
              go to 1071
           else
              cd(j,2)='+'
              go to 1072
           endif
c
 1071           cd(j,2)='-'
c
 1072     continue
c
         if(i.eq.lmax) then
            cd(j,3) = ch2 // ilet(1)
         else 
            cd(j,3) = ch5 // ilet(i+1)
c
         endif
         if(i.eq.1) then
            cd(j,5) = ch4 // olet(lmax)
         else
            cd(j,5) = ch6 // olet(i-1)
         endif
         if(i.eq.lmax) then
            cd(num(1),inlet(1)) = ch1 // 'a'
         else
            cd(num(i+1),inlet(i+1)) = ch1 // 'a'
         endif
         if(i.eq.1) then
            cd(num(lmax),onlet(lmax)) = ch1 // 'c'
         else
            cd(num(i-1),onlet(i-1)) = ch1 // 'c'
         endif
       else 
         if((ith(i,1)*ith(i,2)).eq.1) then
            j = num(i)
c - stef
         write(ch1, '(i3)' ) j
c
       if(i.lt.lmax) then
         write(ch5, '(i3)' ) num(i+1)
       else
         write(ch5, '(i3)' ) num(1)
       endif
c
       if(i.gt.1) then
         write(ch6, '(i3)' ) num(i-1)
       else
         write(ch6, '(i3)' ) num(lmax)
       endif
c - stef
            if(i.eq.lmax) then
              cd(j,6) = ch2 // ilet(1)
            else
              cd(j,6) = ch5 // ilet(i+1)
            endif
            if(i.eq.1) then
              cd(j,4) = ch4 // olet(lmax)
            else
              cd(j,4) = ch6 // olet(i-1)
            endif
            if(i.eq.lmax) then
              cd(num(1), inlet(1)) = ch1 // 'd'
            else
              cd(num(i+1), inlet(i+1)) = ch1 // 'd'
            endif
            if(i.eq.1) then
              cd(num(lmax),onlet(lmax)) = ch1 // 'b'
            else
              cd(num(i-1),onlet(i-1)) = ch1 // 'b'
            endif
         else 
            j = num(i)
c - stef
         write(ch1, '(i3)' ) j
c
       if(i.lt.lmax) then
         write(ch5, '(i3)' ) num(i+1)
       else
         write(ch5, '(i3)' ) num(1)
       endif
c
       if(i.gt.1) then
         write(ch6, '(i3)' ) num(i-1)
       else
         write(ch6, '(i3)' ) num(lmax)
       endif
c - stef
            if(i.eq.lmax) then
              cd(j,4) = ch2 // ilet(1)
            else
              cd(j,4) = ch5 // ilet(i+1)
            endif
            if(i.eq.1) then
              cd(j,6) = ch4 // olet(lmax)
            else
              cd(j,6) = ch6 // olet(i-1)
            endif
            if(i.eq.lmax) then
              cd(num(1), inlet(1)) = ch1 // 'b'
            else
              cd(num(i+1), inlet(i+1)) = ch1 // 'b'
            endif
            if(i.eq.1) then
              cd(num(lmax),onlet(lmax)) = ch1 // 'd'
            else
              cd(num(i-1),onlet(i-1)) = ch1 // 'd'
            endif
         endif 
       endif
 501   continue
c  
c - end of while#5 loop
c
      if(failure.ne.0) then
        go to 999
      else
        junk = 0
      endif
c
c - the next section prints out the matrix and calculates
c - the writhe
c - start of while#6 loop
c
       do 773 ist = 1,lmax/2
       do 777 jst = 1,6
       cdwr(ist,jst,1) = cd(ist,jst)(1:1)
       cdwr(ist,jst,2) = cd(ist,jst)(2:2) 
       cdwr(ist,jst,3) = cd(ist,jst)(3:3)
       cdwr(ist,jst,4) = cd(ist,jst)(4:4)
 777   continue
 773   continue
c
c      write(6,*) 'after diam calculation, cdiam = ',cdiam
c      if(cdiam.lt.sdiam) then
c       write(6,*) ' cdiam = ',cdiam
c       do 1036 i = 1, number
c       points(i,1) = newpts(i,1)
c       points(i,2) = newpts(i,2)
c       points(i,3) = newpts(i,3)
c1036   continue
c      else
        if(cdiam.lt.ndiam) then
c        do 1037 i = 1, number
c        points(i,1) = newpts(i,1)
c        points(i,2) = newpts(i,2)
c        points(i,3) = newpts(i,3)
c1037   continue
         ndiam = cdiam
c        write(6,*) ' ndiam = ',ndiam
        else
        go to 575
        endif
c      endif
c
 3778 format('"#"',i6,'"case = "',i6,'"-"',i6,'"diam =  "',f21.17)
c     write(8,3778) counter,kcm,ish,sdiam
c     do 573 i = 1, number
c     write(8,533) newpts(i,1), newpts(i,2), newpts(i,3)
c573  continue
c
 575  write(7,3779) counter, rgm, cdiam 
c
      totalcross = totalcross + lmax/2
       do 775 ist = 1, lmax/2
          ia = 1
       do 779 jst = 1, 6
       do 781 kst = 1, 4 
       if (cdwr(ist,jst,kst).ne.' ') then
           lett(ist,ia) = cdwr(ist,jst,kst)
           ia = ia + 1
       endif
 781   continue
 779   continue
       ib = ia - 1
       write(7,1793) (lett(ist,lst), lst=1, ib)
 1793  format(24a1)
 775   continue	
c
       write(7,1795)
 1795  format(a1)
c
       if(failure.eq.0) then
          failure = 5 
       else
          junk = 0
       endif
c
c - next section is for when the knot has no crossings,
c - or is not in general position
c
      if(failure.lt.5) then
c
         if(failure.eq.1) then
            write(6,*) 'NO CROSSINGS'
         else
            junk = 0
         endif
c
         if(failure.eq.2) then
            write(6,*) 'Not in General Position'
          else
            junk = 0
         endif
c
         if(failure.eq.3) then
            write(6,*) 'The same point was chosen twice.
     : Execution stopped'
         else
            junk = 0
         endif
c
      endif
c
 999  continue
 1011 continue
 1001 continue
      avcross = totalcross/counter
      avdiam = tdiam/counter
      avrgm = trgm/counter
      write(6,*) 'average number of crossings is  ',avcross
      write(6,*) 'average diameter is ',avdiam
      write(6,*) 'average radius of gyration is ',avrgm
      write(6,*) 'closing minimal diameter is ',ndiam
      close(7)
      close(8)
      close(9)
      stop
      end
      subroutine ssort(n,ra)
      real*8 ra(n), rra 
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end
