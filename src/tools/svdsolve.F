 subroutine svdSolve (A,B,IDIAG,NEQ,FACT,BACK)

 use debugmodule, only: lustat, itest, idebug
!
! Routine to perform Singular Value Decomposition and/or backsubstitution
! of symmetric system of equations.
!     A(NA)        matrix in column compacted form
!    B(NEQ)        right side vector (after backsubstitution column
!            solution)
!    IDIAG(NEQ)    diagonal adresses of A(NA)
!    NEQ        number of equations
!    FACT        .TRUE.  = Factor A
!            .FALSE. = Do not factor A
!    BACK        .TRUE.  = Backsubstitute B(NEQ)
!            .FALSE. = Do not backsubstitute B(NEQ)

! SNGL    parameter (IFLT=4)
 parameter (IFLT=8)
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
!-pass
 dimension A(*),B(NEQ),IDIAG(NEQ)
 logical FACT,BACK
!-locl
 pointer (Plarge,ALARGE), (PW,W), (PV,V), (PX,X), (PTMP,TMP)
 dimension ALARGE(1), W(1), V(1), X(1), TMP(1)
 logical FIRST
 save FIRST,Plarge,PW,PV,PX
 external FLOAT
!-init
 data FIRST/.TRUE./
 N  = NEQ
 M  = NEQ
 NP = NEQ
 MP = NEQ
 luw= 95

 if (FIRST) then
Plarge = MALLOC(NEQ*NEQ * IFLT)
PW     = MALLOC(NEQ     * IFLT)
PV     = MALLOC(NEQ*NEQ * IFLT)
PX     = MALLOC(NEQ     * IFLT)
PTMP   = MALLOC(NEQ     * IFLT)
call openf(luw,'EigenVal.dat','unknown')
FIRST  = .FALSE.
if (iecho.eq.1) then
write(stderr,5) FLOAT((2*NEQ*NEQ+3*NEQ)*IFLT)/1000.
 5        format(3X,'svdsolve: allocated ',F7.1,' kbytes'/ &
       3X,'eigenvalues will be written to "EigenVal.dat"'/)
endif
 endif

 if (FACT) then
call unpack (A,ALARGE,IDIAG,NEQ)
if (idebug.ne.0) then
write(stderr,6)
 6        format(/3X,'unpacked full stiffness matrix:')
call putmat (ALARGE,NEQ,NEQ,stderr)
write(stderr,7)
 7        format(/3X,'load vector:')
call putvec (B,NEQ,stderr)
endif
if (iecho.eq.1) then
write(stderr,10,advance='no')
10        format(/3X,'Starting SVD...')
endif
call svdcmp (ALARGE,N,M,NP,MP,W,V,TMP)
if (iecho.eq.1) write(stderr,15)
15        format('done')
call condit (W,NEQ,luw)
 endif
 if (BACK) then
if (iecho.eq.1) then
write(stderr,20,advance='no')
20        format(/3X,'Backsubstituting SVD solution...')
endif
call svbksb (ALARGE,W,V,N,M,NP,MP,B,X,TMP)
call cpXtoB (X,B,NEQ)
if (iecho.eq.1) write(stderr,15)
 endif
 return
 end
!-------------------------------------------------------------------------------
 subroutine unpack (A,ALARGE,IDIAG,NEQ)
!
! Routine to unpack column compacted matrix A into ALARGE(NEQ,NEQ).
! A contains upper triangle columns, running from top to diagonal
! ("skyline storage" model in Akin [1982]).
!
 implicit double precision (a-h,o-z)
 dimension A(*),ALARGE(NEQ,NEQ),IDIAG(NEQ)

 call CLEAR(ALARGE,NEQ*NEQ)

 ALARGE(1,1) = A(IDIAG(1))
 do 200 j=2,NEQ
jh = IDIAG(j) - IDIAG(j-1)
do 100 i=j-jh+1,j
ij = IDIAG(j) - (j-i)
ALARGE(i,j) = A(ij)
if (i.ne.j) ALARGE(j,i) = A(ij)
100        continue
200    continue
 return
 end
!-------------------------------------------------------------------------------
 subroutine condit (W,NEQ,luw)

 use constants, only: dfmin,dlogmin,dfmax,dlogmax,eps

!
! routine to check conditioning of matrix by checking Eigenvalues range
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
!-pass
 dimension W(NEQ)
 integer NEQ
!-locl
 double precision Wmax,Wmin
 integer i,N
 external FLOAT
 i=1

 write(luw,11) NEQ
11    format('2 ',I6)
 write(luw,12) (FLOAT(i),W(i),i=1,NEQ)
12    format(2(1PG15.6,1X))

 N = 0
 Wmax = -dfmax
 Wmin =  dfmax
 do 100 i=1,NEQ
if (Wmin.gt.W(i)) Wmin = W(i)
100        if (Wmax.lt.W(i)) Wmax = W(i)
 if (Wmax.lt.dfmin) then
write(stderr,1) dfmin
 1        format('condit FATAL error: maximum eigenvalue less than ', &
     1PG12.2)
call stoper()
 endif
 if (iecho.eq.1) then
write(stderr,2) Wmin,Wmax
 2        format(//3X,'singular value decomposition: '/3X, &
     '                 smallest eigenvalue = ',1PE14.6/3X, &
     '                  largest eigenvalue = ',1PE14.6/  )
 endif
 if (OUTFIL(1)) write(luout(1),2) Wmin,Wmax

 Wmin = EPS * Wmax
 do 200 i=1,NEQ
if (Wmin.gt.W(i)) then
N = N + 1
W(i) = 0.
endif
200    continue
 if (iecho.eq.1) then
if (N.gt.0) then
write(stderr,3) N
 3        format(3X,'number of eigenvalues forced to zero = ', &
         I5//)
endif
 endif
 if (OUTFIL(1)) write(luout(1),3) N
 return
 end
!-------------------------------------------------------------------------------
   subroutine svdcmp(a,m,n,mp,np,w,v,rv1)
   USE PAUSEMODULE

! SNGL    implicit real (a-h,o-z,A-H,O-Z)
   implicit double precision (a-h,o-z)
   dimension a(mp,np),w(np),v(np,np),rv1(mp)
   g=0.0
   scale=0.0
   anorm=0.0
   do 25 i=1,n
l=i+1
rv1(i)=scale*g
g=0.0
s=0.0
scale=0.0
if (i.le.m) then
  do 11 k=i,m
    scale=scale+abs(a(k,i))
11        continue
  if (scale.ne.0.0) then
    do 12 k=i,m
      a(k,i)=a(k,i)/scale
      s=s+a(k,i)*a(k,i)
12          continue
    f=a(i,i)
    g=-sign(sqrt(s),f)
    h=f*g-s
    a(i,i)=f-g
    if (i.ne.n) then
      do 15 j=l,n
        s=0.0
        do 13 k=i,m
          s=s+a(k,i)*a(k,j)
13              continue
        f=s/h
        do 14 k=i,m
          a(k,j)=a(k,j)+f*a(k,i)
14              continue
15            continue
    endif
    do 16 k= i,m
      a(k,i)=scale*a(k,i)
16          continue
  endif
endif
w(i)=scale *g
g=0.0
s=0.0
scale=0.0
if ((i.le.m).and.(i.ne.n)) then
  do 17 k=l,n
    scale=scale+abs(a(i,k))
17        continue
  if (scale.ne.0.0) then
    do 18 k=l,n
      a(i,k)=a(i,k)/scale
      s=s+a(i,k)*a(i,k)
18          continue
    f=a(i,l)
    g=-sign(sqrt(s),f)
    h=f*g-s
    a(i,l)=f-g
    do 19 k=l,n
      rv1(k)=a(i,k)/h
19          continue
    if (i.ne.m) then
      do 23 j=l,m
        s=0.0
        do 21 k=l,n
          s=s+a(j,k)*a(i,k)
21              continue
        do 22 k=l,n
          a(j,k)=a(j,k)+s*rv1(k)
22              continue
23            continue
    endif
    do 24 k=l,n
      a(i,k)=scale*a(i,k)
24          continue
  endif
endif
anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
   do 32 i=n,1,-1
if (i.lt.n) then
  if (g.ne.0.0) then
    do 26 j=l,n
      v(j,i)=(a(i,j)/a(i,l))/g
26          continue
    do 29 j=l,n
      s=0.0
      do 27 k=l,n
        s=s+a(i,k)*v(k,j)
27            continue
      do 28 k=l,n
        v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
  endif
  do 31 j=l,n
    v(i,j)=0.0
    v(j,i)=0.0
31        continue
endif
v(i,i)=1.0
g=rv1(i)
l=i
32    continue
   do 39 i=n,1,-1
l=i+1
g=w(i)
if (i.lt.n) then
  do 33 j=l,n
    a(i,j)=0.0
33        continue
endif
if (g.ne.0.0) then
  g=1.0/g
  if (i.ne.n) then
    do 36 j=l,n
      s=0.0
      do 34 k=l,m
        s=s+a(k,i)*a(k,j)
34            continue
      f=(s/a(i,i))*g
      do 35 k=i,m
        a(k,j)=a(k,j)+f*a(k,i)
35            continue
36          continue
  endif
  do 37 j=i,m
    a(j,i)=a(j,i)*g
37        continue
else
  do 38 j= i,m
    a(j,i)=0.0
38        continue
endif
a(i,i)=a(i,i)+1.0
39    continue
   do 49 k=n,1,-1
do 48 its=1,30
  do 41 l=k,1,-1
    nm=l-1
    if ((abs(rv1(l))+anorm).eq.anorm)  go to 2
    if ((abs(w(nm))+anorm).eq.anorm)  go to 1
41        continue
1         c=0.0
  s=1.0
  do 43 i=l,k
    f=s*rv1(i)
    if ((abs(f)+anorm).ne.anorm) then
      g=w(i)
      h=sqrt(f*f+g*g)
      w(i)=h
      h=1.0/h
      c= (g*h)
      s=-(f*h)
      do 42 j=1,m
        y=a(j,nm)
        z=a(j,i)
        a(j,nm)=(y*c)+(z*s)
        a(j,i)=-(y*s)+(z*c)
42            continue
    endif
43        continue
2         z=w(k)
  if (l.eq.k) then
    if (z.lt.0.0) then
      w(k)=-z
      do 44 j=1,n
        v(j,k)=-v(j,k)
44            continue
    endif
    go to 3
  endif
  if (its.eq.30) go_or_quit('no convergence in 30 iterations')
  x=w(l)
  nm=k-1
  y=w(nm)
  g=rv1(nm)
  h=rv1(k)
  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
  g=sqrt(f*f+1.0)
  f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
  c=1.0
  s=1.0
  do 47 j=l,nm
    i=j+1
    g=rv1(i)
    y=w(i)
    h=s*g
    g=c*g
    z=sqrt(f*f+h*h)
    rv1(j)=z
    c=f/z
    s=h/z
    f= (x*c)+(g*s)
    g=-(x*s)+(g*c)
    h=y*s
    y=y*c
    do 45 nm=1,n
      x=v(nm,j)
      z=v(nm,i)
      v(nm,j)= (x*c)+(z*s)
      v(nm,i)=-(x*s)+(z*c)
45          continue
    z=sqrt(f*f+h*h)
    w(j)=z
    if (z.ne.0.0) then
      z=1.0/z
      c=f*z
      s=h*z
    endif
    f= (c*g)+(s*y)
    x=-(s*g)+(c*y)
    do 46 nm=1,m
      y=a(nm,j)
      z=a(nm,i)
      a(nm,j)= (y*c)+(z*s)
      a(nm,i)=-(y*s)+(z*c)
46          continue
47        continue
  rv1(l)=0.0
  rv1(k)=f
  w(k)=x
48      continue
3       continue
49    continue
   return
   end
!-------------------------------------------------------------------------------
   subroutine svbksb(u,w,v,m,n,mp,np,b,x,tmp)
! SNGL    implicit real (a-h,o-z,A-H,O-Z)
   implicit double precision (a-h,o-z)
   dimension u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(mp)
   do 12 j=1,n
s=0.
if(w(j).ne.0.)then
  do 11 i=1,m
    s=s+u(i,j)*b(i)
11        continue
  s=s/w(j)
endif
tmp(j)=s
12    continue
   do 14 j=1,n
s=0.
do 13 jj=1,n
  s=s+v(j,jj)*tmp(jj)
13      continue
x(j)=s
14    continue
   return
   end
!-------------------------------------------------------------------------------
 subroutine cpXtoB (X,B,NEQ)
! SNGL    implicit real (a-h,o-z,A-H,O-Z)
 implicit double precision (a-h,o-z)
 dimension X(NEQ),B(NEQ)
 do 100 i=1,NEQ
100        B(i) = X(i)
 return
 end
