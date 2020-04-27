************************************************************************
*                           TRIAL FUNCTION                             *
************************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calbvec( l,theta,phi,plm,bvec,bvecdt,bvecdp )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Evaluating the value of toroidal harmonics (fully normalized)
c at each station whose latitude and longitude are theta and phi.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit none
	real*8 pi
	parameter ( pi=3.1415926535897932d0 )
c
	integer l,m,i
	real*8 theta,phi,x,plm(3,0:3),fact,coef
	complex*16 bvec(3,-2:2),expimp
	complex*16 bvecdt(3,-2:2),bvecdp(3,-2:2)
	real*8 plmdt,xl2
c
	x = dcos( theta )
	xl2 = dble(l) * dble(l+1)
	do 100 m=0,min0(l,3)
	  call calplm( l,m,x,plm(1,m) )
  100	continue
	do 120 m=0,min0(l,2)
	  fact = 1.d0
	  if ( m.ne.0 ) then
	    do 110 i=l-m+1,l+m
	      fact = fact * dble(i)
  110	    continue
	  endif
	  coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
	  expimp = cdexp( dcmplx( 0.d0, dble(m)*phi ) )
	  plmdt = dble(m) * x / dsin( theta ) * plm(1,m) + plm(1,m+1)
c
	  bvec(1,m)  = dcmplx( 0.d0 )
	  bvec(1,-m) = dcmplx( 0.d0 )
	  bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta )
     &	               * coef * plm(1,m) * expimp
	  bvec(2,-m) = dconjg( bvec(2,m) )
	  bvec(3,m) = - coef * plmdt * expimp
c	  bvec(3,m) = - coef * ( dble(m) * x / dsin( theta ) * plm(1,m)
c     &	                         + plm(1,m+1) )
c     &	                * expimp
	  bvec(3,-m) = dconjg( bvec(3,m) )
c calculate derivatives
	  bvecdt(1,m)  = dcmplx( 0.d0 )
	  bvecdt(1,-m) = dcmplx( 0.d0 )
	  bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt / dsin(theta)
     &      - x / ( 1 - x * x ) * plm(1,m) ) * coef * expimp
        bvecdt(2,-m) = dconjg( bvecdt(2,m) )
	  bvecdt(3,m) = ( x / dsin(theta) * plmdt 
     &           - dble(m) * dble(m) / (1 - x * x) * plm(1,m)
     &           + xl2 * plm(1,m) ) * coef * expimp
        bvecdt(3,-m) = dconjg( bvecdt(3,m) )
	  bvecdp(1,m)  = dcmplx( 0.d0 )
	  bvecdp(1,-m) = dcmplx( 0.d0 )
	  bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m)
     &              * coef * expimp
        bvecdp(2,-m) = dconjg( bvecdp(2,m) )
	  bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef * expimp
        bvecdp(3,-m) = dconjg( bvecdp(3,m) )
c
	  if ( mod(m,2).eq.1 ) then
	    bvec(2,-m) = - bvec(2,-m)
	    bvec(3,-m) = - bvec(3,-m)
	    bvecdt(2,-m) = - bvecdt(2,-m)
	    bvecdt(3,-m) = - bvecdt(3,-m)
	    bvecdp(2,-m) = - bvecdp(2,-m)
	    bvecdp(3,-m) = - bvecdp(3,-m)
	  endif
  120	continue
c
	return
	end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calbvecreal( l,theta,plm,bvec,bvecdt,bvecdp )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Evaluating the value of toroidal harmonics (fully normalized)
c at each station whose latitude and longitude are theta and phi.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      integer l,m,i
      real*8 theta,phi,x,plm(3,0:3),fact,coef
      complex*16 bvec(3,0:2)
      complex*16 bvecdt(3,0:2),bvecdp(3,0:2)
      real*8 plmdt,xl2
      real*8 sintm
c
      x = dcos( theta )
      sintm=1.d0/dsin(theta)
      xl2 = dble(l) * dble(l+1)
c
      do 100 m=0,min0(l,3)
      call calplm( l,m,x,plm(1,m) )
  100    continue
      do 120 m=0,min0(l,2)
      fact = 1.d0
      if ( m.ne.0 ) then
      do 110 i=l-m+1,l+m
      fact = fact * dble(i)
  110        continue
      endif
c
      coef = dcmplx(dsqrt( dble(2*l+1)/(4.d0*pi) / fact ))
      plmdt = dble(m) * x * sintm * plm(1,m) + plm(1,m+1)
c
      bvec(1,m)  = dcmplx(0.d0)
      bvec(2,m)  = dcmplx( 0.d0, dble(m) ) * sintm
     &                   * coef * plm(1,m)
      bvec(3,m) = - coef * plmdt
c calculate derivatives
      bvecdt(1,m)  = dcmplx(0.d0)
      bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt * sintm
     &      - x / ( 1 - x * x ) * plm(1,m) ) * coef
      bvecdt(3,m) = (x * sintm * plmdt
     &           - dble(m) * dble(m) / (1 - x * x) * plm(1,m)
     &           + xl2 * plm(1,m) ) * coef
      bvecdp(1,m)  = dcmplx(0.d0)
      bvecdp(2,m)  = - dble(m) * dble(m) * sintm * plm(1,m)
     &              * coef
      bvecdp(3,m)  = -dcmplx(0.d0, dble(m)) * plmdt * coef
c
  120    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calbvecreal2( l,theta,plm,bvec,bvecdt,bvecdp )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Evaluating the value of toroidal harmonics (fully normalized)
c at each station whose latitude and longitude are theta and phi.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      integer l,m,i
      real*8 theta,x,plm(3,0:3),fact,coef
      complex*16 bvec(3,-2:2)
      complex*16 bvecdt(3,-2:2),bvecdp(3,-2:2)
      real*8 plmdt,xl2
      real*8 sintm
c
      x = dcos( theta )
      sintm=1.d0/dsin(theta)
      xl2 = dble(l) * dble(l+1)
c
      do 100 m=0,min0(l,3)
      call calplm( l,m,x,plm(1,m) )
  100    continue
      do 120 m=0,min0(l,2)
      fact = 1.d0
      if ( m.ne.0 ) then
      do 110 i=l-m+1,l+m
      fact = fact * dble(i)
  110        continue
      endif
c
      coef = dcmplx(dsqrt( dble(2*l+1)/(4.d0*pi) / fact ))
      plmdt = dble(m) * x * sintm * plm(1,m) + plm(1,m+1)
c
      bvec(1,m)  = dcmplx(0.d0)
      bvec(1,-m)  = dcmplx(0.d0)
      bvec(2,m)  = dcmplx( 0.d0, dble(m) ) * sintm
     &                   * coef * plm(1,m)
      bvec(2,-m)  = dconjg(bvec(2,m))
      bvec(3,m) = - coef * plmdt
      bvec(3,-m)  = dconjg(bvec(3,m))
c calculate derivatives
      bvecdt(1,m)  = dcmplx(0.d0)
      bvecdt(1,-m) = dcmplx(0.d0)
      bvecdt(2,m) = dcmplx( 0.d0, dble(m) ) * ( plmdt * sintm
     &      - x / ( 1 - x * x ) * plm(1,m) ) * coef
      bvecdt(2,-m) = dconjg(bvecdt(2,m))
      bvecdt(3,m) = (x * sintm * plmdt
     &           - dble(m) * dble(m) / (1 - x * x) * plm(1,m)
     &           + xl2 * plm(1,m) ) * coef
      bvecdt(3,-m) = dconjg(bvecdt(3,m))
      bvecdp(1,m)  = dcmplx(0.d0)
      bvecdp(1,-m)  = dcmplx(0.d0)
      bvecdp(2,m)  = - dble(m) * dble(m) * sintm * plm(1,m)
     &              * coef
      bvecdp(2,-m)  = dconjg(bvecdp(2,m))
      bvecdp(3,m)  = -dcmplx(0.d0, dble(m)) * plmdt * coef
      bvecdp(3,-m)  = dconjg(bvecdp(3,m))
c
      if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
      endif
c
  120    continue
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine calplm( l,m,x,plm )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer l,m,i
	real*8 x,plm(3),pmm,somx2,fact
c
	if ( (m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.d0) )
     &	  pause 'bad arguments'
	if ( l.eq.m ) then
	  pmm = 1.d0
	  if ( m.gt.0 ) then
	    somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
	    fact = 1.d0
	    do 11 i=1,m
	      pmm = -pmm * fact * somx2
	      fact = fact + 2.d0
   11	    continue
	  endif
	  plm(3) = 0.d0
	  plm(2) = 0.d0
	  plm(1) = pmm
	else
	  plm(3) = plm(2)
	  plm(2) = plm(1)
	  if ( l.eq.m+1 ) then
	    plm(1) = x * dble(2*m+1) * plm(2)
	  else
	    plm(1) = ( x * dble(2*l-1) * plm(2)
     &	                 - dble(l+m-1) * plm(3)
     &	                ) / dble(l-m)
	  endif
	endif
c
	return
	end
