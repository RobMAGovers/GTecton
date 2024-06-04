 function Tcentr (T,TX,XL,IEN,LMTX,IDUM1)

 USE MODELDEFINITION

!
! Calculates the temperature in the center of the element
!
 implicit none
 integer NEN,NSD
 parameter (NEN=4, NSD=3)
 double precision QUART
 parameter (QUART=0.25d0)
!-pass
 integer IEN,LMTX,IDUM1
 double precision T,TX,XL,Tcentr
 dimension T(*),TX(*),IEN(NEN),LMTX(NEN),XL(NSD,NEN)

 integer ITMODE,NTCALC,ITCALC,MAXTIT,IADVEC,IDIFT
 integer NTYPE,NUMEL,NUMAT,IOPT,IPOINT,NUMPR,NUMSLP, &
   NUMFN,NPRE,LGDEF,IRESDU,NUMSTR,IGRAV,IVLIM, &
   NUMWNK,NSURF,NSED,INCOMP, NOCOMPR,NSLSKEW

!-locl
 double precision TL
 dimension TL(NEN)
!
!    Localize temperatures
 call LTEMP (IEN,T,TL,NEN)
!
!    Modify for discontinuities
 if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) &
   call LTEMPX (LMTX,IEN,TX,TL,NEN)
!
!    Calculate center temperature
 Tcentr = (TL(1)+TL(2)+TL(3)+TL(4))*QUART
!
 return
 end
