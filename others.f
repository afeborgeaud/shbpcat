c others.f for wcalprem.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pinput( maxnlay,maxnzone,maxnr,maxnstack,
     &     re,ratc,ratl,
     &     tlen,np,omegai,imin,imax,
     &     nzone,vrmin,vrmax,rho,vsv,vsh,qmu,
     &     nr,theta,phi,lat,lon,dir,
     &     obslat,obslon,obs,nsta,rsta ) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Parameter Input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer maxnlay,maxnzone,maxnr,maxnstack
      integer np
      integer imin,imax
      integer nzone,nr
      real*8 thetamin,thetamax,dtheta
      real*8 tlen,omegai,re,ratc,ratl
      real*8 vrmin(*),vrmax(*),rho(4,*),vsv(4,*),vsh(4,*)
      real*8 qmu(*)
      real*8 theta(*),phi(*),lat(*),lon(*),stlat,stlon
      real*8 obslat,obslon
      integer i
      integer nsta
      real*8 rsta(*),obslattmp,r0
      character*80 dummy,tmpfile,obs
      character*80 dir
c
      data tmpfile / 'workshbp' /
c
c temporary file open
      open( unit=11, file=tmpfile, status='unknown' )
c writing to the temporary file
 100  continue
      read(5,110) dummy
 110  format(a80)
      if ( dummy(1:1).eq.'c' ) goto 100
      if ( dummy(1:3).eq.'end' ) goto 120
      write(11,110) dummy
      goto 100
 120  continue
c temporary file close
      close(11)
c 
c temporary file open
      open( unit=11, file=tmpfile, status='unknown' )
c reading the parameter
      read(11,*) tlen,np
      read(11,*) re		! relative error (vertical grid)
      read(11,*) ratc		! ampratio (vertical grid cut-off)
      read(11,*) ratl		! ampratio (for l-cutoff)
      read(11,*) omegai
      omegai = - dlog(omegai) / tlen
      read(11,*) imin,imax
      read(11,*) nzone
      if ( nzone.gt.maxnzone )
     &     pause 'nzone is too large. (pinput)'
      do 140 i=1,nzone
         read(11,*) vrmin(i),vrmax(i),
     &        rho(1,i),rho(2,i),rho(3,i),rho(4,i),
     &        vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i),
     &        vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), qmu(i)
 140  continue
c
      read(11,*) r0,obslat,obslon
      obslattmp = obslat
      call translat(obslattmp,obslattmp)
c     
      read(11,110) dir
c
      read(11,110) obs ! name of station
c
      read(11,*) thetamin,thetamax,dtheta
      nr=int((thetamax - thetamin)/dtheta) + 1
      write(*,'(a15 f8.2 f8.2 f8.2 i7)') "Theta sampling ",
     & thetamin,thetamax,dtheta,nr
      if ( nr.gt.maxnr )
     &     pause 'nr is too large. (pinput)'
      do 150 i=1,nr
        theta(i)=(i-1)*dtheta + thetamin
        phi(i)=0.d0
 150  continue
      read(11,*) nsta
      if(nsta.gt.maxnstack)  
     &     pause 'nsdta is toolarge (pinput)'
      do 170 i=1,nsta
         read(11,*) rsta(i)
 170  continue
c     temporary file close
      close(11)
c     
      return
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 pi
      parameter ( pi = 3.1415926535897932d0 )
c       
      real*8 ievla,ievlo,istla,istlo
      real*8 evla,evlo,stla,stlo
      real*8 theta,phi
      real*8 gcarc,az
      real*8 tc,ts
c
c transformation to spherical coordinates
c
      evla = 90.d0 - ievla
      stla = 90.d0 - istla
c
      evla = evla / 1.8d2 * pi
      evlo = ievlo / 1.8d2 * pi
      stla = stla / 1.8d2 * pi
      stlo = istlo / 1.8d2 * pi
c  
      gcarc = dacos( dcos(evla) * dcos(stla)
     &     + dsin(evla) * dsin(stla) * dcos(evlo - stlo) )
c
      tc = ( dcos(stla) * dsin(evla) 
     &     - dsin(stla) * dcos(evla) * dcos(stlo - evlo) )
     &     / dsin(gcarc)
      ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)
c
      az = dacos(tc)
      if( ts .lt. 0.d0 ) az = -1.d0 * az
c
      az = az * 1.8d2 / pi

      gcarc = gcarc * 1.8d2 / pi
c
      theta = gcarc
      phi   = 180.d0 - az
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine translat(geodetic,geocentric)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 flattening,pi
c      parameter ( flattening = 1.d0 / 297.d0 )
c      parameter ( flattening = 1.d0 / 298.25d0)
      parameter(flattening=0.d0)
      parameter ( pi = 3.1415926535897932d0 )
      real*8 geocentric, geodetic

      integer flag
c      read(5,*) geodetic
      flag = 0
      if(geodetic .gt. 90.d0) then
         geodetic = 1.8d2 - geodetic
         flag = 1
      endif
c
      geodetic = geodetic / 1.8d2 * pi
      geocentric = datan( (1.d0 - flattening) * (1.d0 - flattening)
     &     * dtan(geodetic) )
      geocentric = geocentric * 1.8d2 / pi
c      if(geocentric .lt. 0.d0 ) geocentric = 1.8d2 + geocentric
      if(flag .eq. 1) then
         geocentric = 1.8d2 - geocentric
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calgrid( nzone,vrmin,vrmax,vs,rmin,rmax,
     &     imax,lmin,tlen,vmin,gridpar,dzpar )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c
      integer nzone,imax,lmin
      real*8 vrmin(*),vrmax(*),vs(4,*)
      real*8 rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
      integer izone,i,j
      real*8 coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp
c     
      do 130 izone=1,nzone
c     computing the S-velocity at each zone
         do 110 i=1,4
            v(i) = vs(i,izone)
 110     continue
         vs1 = 0.d0
         vs2 = 0.d0
         do 120 j=1,4
	    if ( j.eq.1 ) then
               coef1 = 1.d0
            else
               coef1 = coef1 * ( vrmin(izone) / rmax )
	    endif
	    if ( j.eq.1 ) then
               coef2 = 1.d0
            else
               coef2 = coef2 * ( vrmax(izone) / rmax )
	    endif
	    vs1 = vs1 + v(j) * coef1
	    vs2 = vs2 + v(j) * coef2
 120     continue
c     computing rh
         rh = vrmax(izone) - vrmin(izone)
c     computing omega,amax
         omega = 2.d0 * pi * dble(imax) / tlen
         if ( vs1.ge.vs2 ) then
	    vmin(izone) = vs2
         else
	    vmin(izone) = vs1
         endif
         amax = vrmax(izone)
         gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) ) 
     &        - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )
     &        / ( amax * amax )
         if ( gtmp.gt.0.d0 ) then
	    dzpar(izone)   = dsqrt( 1.d0/gtmp )
	    gridpar(izone) = rh / dzpar(izone)
         else
	    dzpar(izone)   = 0.d0
	    gridpar(izone) = 0.d0
         endif
 130  continue
c     rearangement of gridpar
      gtmp = 0.d0
      do 140 izone=1,nzone
         gtmp = gtmp + gridpar(izone)
 140  continue
      do 150 izone=1,nzone
         if ( gridpar(izone).gt.0.d0 ) then
	    gridpar(izone) = gridpar(izone) / gtmp
         else
	    rh = vrmax(izone) - vrmin(izone)
	    gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
         endif
 150  continue
c     re-rearangement of gridpar
      gtmp = 0.d0
      do 160 izone=1,nzone
         gtmp = gtmp + gridpar(izone)
 160  continue
      do 170 izone=1,nzone
         gridpar(izone) = gridpar(izone) / gtmp
 170  continue
c     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calra( maxnstack,maxnlay,maxnzone,
     &     nlayer,
     &     gridpar,dzpar,nzone,vrmin,vrmax,
     &     rmin,rmax,nnl,ra,re,nsta,rsta,rrsta,iista,
     &     rs,ciista )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the number and the location of grid points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 pi
      parameter ( pi=3.1415926535897932d0 )
c     
      integer maxnstack,maxnlay,maxnzone
      integer nlayer,inlayer
      integer nzone,nnl(maxnzone)
      real*8 gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,r0
      real*8 ra(maxnlay+maxnzone+1),rs
      integer izone,itmp,i,ntmp
      real*8 ctmp               ! distance betwee source and the nearst
      real*8 rh,re
c  
      integer nsta
      real*8 rsta(maxnstack),rrsta(3,maxnstack)
      integer iista(3,maxnstack),ista,ciista
c
      ctmp = 7000.d0
c Initializing the data
      inlayer = 0
      do 100 i=1,maxnlay+maxnzone+1
         ra(i) = 0.d0
 100  continue
      do 110 izone=1,nzone
         nnl(izone) = 0
 110  continue
c
c computing the number and the location of the grid points
      ra(1) = rmin
      itmp = 1
      do 140 izone=1,nzone
         rh = vrmax(izone) - vrmin(izone)
         if(dzpar(izone).eq.0.d0) then
            ntmp = 1
         else
            ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) 
     &           / 2.d0 / pi  / 7.d-1 + 1 )
         endif
c                             ! ntmp (see Geller & Takeuchi 1995 6.2)
         nnl(izone) = ntmp
         if ( nnl(izone).lt.5 ) nnl(izone)=5
         do 130 i=1,nnl(izone)
	    itmp = itmp + 1
	    ra(itmp) = vrmin(izone)
     &           + rh * dble(i) / dble( nnl(izone) )
 130     continue
 140  continue
c
      itmp = 1
      do 240 izone=1,nzone
         do 230 i=1,nnl(izone)
            do 200 ista=1,nsta
               if( (ra(itmp).lt.rsta(ista))
     &              .and.(rsta(ista).le.ra(itmp+1)))
     &              then
                  if(i.ne.nnl(izone)) then
                     rrsta(1,ista) = ra(itmp)
                     rrsta(2,ista) = ra(itmp+1)
                     rrsta(3,ista) = ra(itmp+2)
c     
                     iista(1,ista) = itmp
                     iista(2,ista) = itmp + 1
                     iista(3,ista) = itmp + 2
                  else
                     rrsta(1,ista) = ra(itmp-1)
                     rrsta(2,ista) = ra(itmp)
                     rrsta(3,ista) = ra(itmp+1)
c     
                     iista(1,ista) = itmp - 1
                     iista(2,ista) = itmp
                     iista(3,ista) = itmp + 1
                  endif
                  if(dabs(rs-rsta(ista)).lt.ctmp) then
                     ciista = itmp
	             ctmp = dabs(rs-rsta(ista))
                  endif
               endif
 200        continue
            itmp = itmp + 1
 230     continue
 240  continue
c recouting the total number of grid points
      inlayer = 0
      do 150 izone=1,nzone
         inlayer = inlayer + nnl(izone)
 150  continue
      nlayer = inlayer
c     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calgra( isp,ra,r0,spn,spo,gra )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer isp(*),spn,itmp
      real*8 ra(*),r0,spo,gra(*)
c     
      itmp = isp(spn) + dint( spo )
      gra(1) = ra(itmp)
      gra(2) = r0
      gra(3) = ra(itmp+1)
c     
      return
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calsp( ndc,nlayer,isp,jsp )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the stack points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer ndc,nlayer(*)
      integer isp(*),jsp(*)
      integer i
c
c computation of isp,jsp,ksp,lsp
      isp(1) = 1
      jsp(1) = 1
      do 100 i=1,ndc
         isp(i+1) = isp(i) + nlayer(i)
         jsp(i+1) = jsp(i) + 4 * nlayer(i)
 100  continue
c     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calspo( ndc,rdc,nlayer,r0,rmin,rmax,ra,
     &     isp,spo,spn )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the source location.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer ndc,nlayer,isp(*),spn
      real*8 rdc(*),r0,rmin,rmax,ra(*),spo
      integer itmp
c
c checking the parameter
      if ( (r0.lt.rmin).or.(r0.gt.rmax) )
     &     pause 'The source location is improper.(calspo)'
c computing 'spo'
      if ( r0.eq.rmax ) then
         spo = dble(nlayer) - 0.01d0
         r0 = ra(nlayer)
     &        + (spo-dble(nlayer-1)) * ( ra(nlayer+1)-ra(nlayer) )
c	  write(6,*) 'r0 is changed to ',r0,spo
      else
         itmp = 2
 110     continue
         if ( r0.lt.ra(itmp) ) then
            continue
         else
            itmp = itmp + 1
            goto 110
         endif
         spo = dble(itmp-2)
     &        + ( r0-ra(itmp-1) ) / ( ra(itmp)-ra(itmp-1) )
c     temporal handling
         if ( (spo-dble(itmp-2)).lt.0.01d0 ) then
	    spo = dble(itmp-2) + 0.01d0
	    r0 = ra(itmp-1)
     &	         + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
c	    write(6,*) 'r0 is changed to ',r0,spo
         endif
         if ( (spo-dble(itmp-2)).gt.0.99d0 ) then
	    spo = dble(itmp-2) + 0.99d0
	    r0 = ra(itmp-1)
     &	         + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
c	    write(6,*) 'r0 is changed to ',r0,spo
         endif
c     
      endif
c computing 'spn'
      spn = 0
      itmp = 1
 120  continue
      spn = spn + 1
      if ( r0.le.rdc(itmp) ) then
         continue
      else
         itmp = itmp + 1
         goto 120
      endif
c changing 'spo'
      spo = spo - dble( isp(spn) - 1 )
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calstg( nzone,rrho,vsv,vsh,
     &     nlayer,nnl,ra,rmax,vnp,vra,rho,
     &     ecL,ecN )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the structure grid points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer nzone,nlayer,nnl(*),vnp
      real*8 rrho(4,*),vsv(4,*),vsh(4,*),ra(*),rmax
      real*8 vra(*),rho(*),ecL(*),ecN(*)
      real*8 trho,tvsv,tvsh,coef
      integer izone,i,j,itmp,jtmp
c     
c     initializing the data
      do 100 i=1,nlayer+nzone+1
         vra(i) = 0.d0
         rho(i) = 0.d0
         ecL(i) = 0.d0
         ecN(i) = 0.d0
 100  continue
c computing the structure grid points
      itmp = 0
      jtmp = 0
      do 130 izone=1,nzone
         do 120 i=1,nnl(izone)+1
	    itmp = itmp + 1
	    jtmp = jtmp + 1
	    vra(itmp) = ra(jtmp)
c --- evaluating the density and elastic constants at this point
	    trho = 0.d0
	    tvsv = 0.d0
	    tvsh = 0.d0
	    do 110 j=1,4
               if ( j.eq.1 ) then
                  coef = 1.d0
               else
                  coef = coef * ( vra(itmp) / rmax )
               endif
               trho = trho + rrho(j,izone) * coef
               tvsv  = tvsv  + vsv(j,izone)   * coef
               tvsh  = tvsh  + vsh(j,izone)   * coef
 110        continue
	    rho(itmp) = trho
	    ecL(itmp)  = rho(itmp) * tvsv * tvsv
	    ecN(itmp)  = rho(itmp) * tvsh * tvsh
 120     continue
         jtmp = jtmp - 1
 130  continue
      vnp = itmp
c     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calgstg( spn,rrho,vsv,vsh,
     &     ra,vra,rmax,rho,ecL,ecN,r0,mu0 )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Computing the structure grid points.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer spn
      real*8 rrho(4,*),vsv(4,*),vsh(4,*)
      real*8 ra(*),rmax
      real*8 vra(*),rho(*),ecL(*),ecN(*),r0,mu0
      real*8 trho,tvsv,tvsh,coef
      integer i,j
c
c initializing the data
      do 100 i=1,3
         vra(i) = 0.d0
         rho(i) = 0.d0
         ecL(i) = 0.d0
         ecN(i) = 0.d0
 100  continue
c     computing the structure grid points
      do 120 i=1,3
         vra(i) = ra(i)
c --- evaluating the density and elastic constants at this point
         trho = 0.d0
         tvsv = 0.d0
         tvsh = 0.d0
         do 110 j=1,4
            if ( j.eq.1 ) then
               coef = 1.d0
            else
               coef = coef * ( vra(i) / rmax )
            endif
            trho = trho + rrho(j,spn) * coef
            tvsv  = tvsv  + vsv(j,spn)   * coef
            tvsh  = tvsh  + vsh(j,spn)   * coef
 110     continue
         rho(i) = trho
         ecL(i)  = rho(i) * tvsv * tvsv
         ecN(i)  = rho(i) * tvsh * tvsh
 120  continue
c     
      mu0 = ecL(2)
c     
      return
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcoef( nzone,omega,q,coef )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 pi
      parameter ( pi = 3.1415926535897932d0 )
c     
      integer izone,nzone
      real*8 omega,q(*)
      complex*16 coef(*)
      real*8 aa,bb
      real*8 omega0
c
      omega0=1.d0 * 2.d0 * pi
c     
      do 100 izone=1,nzone
         if ( omega.eq.0.d0 ) then
	    aa = 1.d0
         else
	    aa = 1.d0
     &	         + dlog( omega / ( omega0 ) ) / ( pi * Q(izone) )
         endif
         bb = 1.d0 / ( 2.d0 * Q(izone) )
         coef(izone) = dcmplx( aa, bb ) * dcmplx( aa, bb )
 100  continue
c     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calu( c0,lsq,bvec,u )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 lsq
      complex*16 c0,bvec(3),u(3)
c     
      u(1) = dcmplx( 0.d0 )
      u(2) = u(2) + c0 * bvec(2) / dcmplx(lsq)
      u(3) = u(3) + c0 * bvec(3) / dcmplx(lsq)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calulcd( c0,c0der,lsq,rsta,theta,
     &     bvec,bvecdt,bvecdp,ulcd )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	
      implicit none
      real*8 lsq,rsta,theta
      complex*16 c0,c0der,bvec(3),bvecdt(3),bvecdp(3),ulcd(9)
c
      complex*16 u1,u2,u3
      complex*16 uder11,uder12,uder13
      complex*16 uder21,uder22,uder23
      complex*16 uder31,uder32,uder33
c     
      u1 = dcmplx(0.d0)
      u2 = c0 * bvec(2) / dcmplx(lsq)
      u3 = c0 * bvec(3) / dcmplx(lsq)
c     
      uder11 = dcmplx(0.d0)
      uder12 = dcmplx(0.d0)
      uder13 = dcmplx(0.d0)
      uder21 = c0der * bvec(2) / dcmplx(lsq)
      uder22 = c0 * bvecdt(2) / dcmplx(lsq)
      uder23 = c0 * bvecdp(2) / dcmplx(lsq)
      uder31 = c0der * bvec(3) / dcmplx(lsq)
      uder32 = c0 * bvecdt(3) / dcmplx(lsq)
      uder33 = c0 * bvecdp(3) / dcmplx(lsq)
c     
      
c      ulcd(1) = ulcd(1) + uder11
      ulcd(1) = dcmplx(0.d0)
      ulcd(2) = ulcd(2) + ( uder12 - u2 ) / rsta 
      ulcd(3) = ulcd(3) + ( uder13 / dsin(theta) - u3 ) / rsta
      ulcd(4) = ulcd(4) + uder21
      ulcd(5) = ulcd(5) + ( uder22 + u1 ) / rsta
      ulcd(6) = ulcd(6) + ( uder23 - u3 * dcos(theta) ) 
     &     / rsta / dsin(theta)
      ulcd(7) = ulcd(7) + uder31
      ulcd(8) = ulcd(8) + uder32 / rsta
      ulcd(9) = ulcd(9) 
     &     + ( ( uder33 + u2 * dcos(theta) ) / dsin(theta) + u1 )
     &     / rsta
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calushlcd( c0,c0der,lsq,rsta,theta,
     &     bvec,bvecdt,bvecdp,ushlcd )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      real*8 lsq,rsta,theta
      complex*16 c0,c0der
      complex*16 bvec(3),bvecdt(3),bvecdp(3)
      complex*16 ushlcd(9)
c
      complex*16 u1,u2,u3
      complex*16 uder11,uder12,uder13
      complex*16 uder21,uder22,uder23
      complex*16 uder31,uder32,uder33
c
      real*8 sintm,cost,rstam
      complex*16 lsqm
c
      cost=dcos(theta)
      sintm=1.d0/dsin(theta)
      rstam=1.d0/rsta
      lsqm=dcmplx(1.d0/lsq)
c
      u1 = dcmplx(0.d0)
      u2 = c0 * bvec(2) * lsqm
      u3 = c0 * bvec(3) * lsqm
c
c uderij is (\phi_i),j with i,j=(r,\theta,\phi)
      uder11 = dcmplx(0.d0)
      uder12 = dcmplx(0.d0)
      uder13 = dcmplx(0.d0)
      uder21 = c0der * bvec(2) * lsqm
      uder22 = c0 * bvecdt(2) * lsqm
      uder23 = c0 * bvecdp(2) * lsqm
      uder31 = c0der * bvec(3) * lsqm
      uder32 = c0 * bvecdt(3) * lsqm
      uder33 = c0 * bvecdp(3) * lsqm
c

c     1:r,r 2:r,t 3:r,p 4:t,r 5:t,t 6:t,p 7:p,r 8:p,t 9:p,p
      ushlcd(1) = dcmplx(0.d0)
      ushlcd(2) = ushlcd(2) + ( uder12 - u2 ) * rstam
      ushlcd(3) = ushlcd(3) + ( uder13 * sintm - u3 ) * rstam
      ushlcd(4) = ushlcd(4) + uder21
      ushlcd(5) = ushlcd(5) + ( uder22 + u1 ) * rstam
      ushlcd(6) = ushlcd(6) + ( uder23 - u3 * cost )
     &     * rstam * sintm
      ushlcd(7) = ushlcd(7) + uder31
      ushlcd(8) = ushlcd(8) + uder32 * rstam
      ushlcd(9) = ushlcd(9)
     &     + ( ( uder33 + u2 * cost ) * sintm + u1 )
     &     * rstam
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matinit( n1,n2,a )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer n1,n2,i,j
      real*8 a(n1,*)
c     
      do 110 j=1,n2
         do 100 i=1,n1
	    a(i,j) = 0.d0
 100     continue
 110  continue
c     
      return
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmatinit( n1,n2,a )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer n1,n2,i,j
      complex*16 a(n1,*)
c     
      do 110 j=1,n2
         do 100 i=1,n1
	    a(i,j) = dcmplx( 0.d0 )
 100     continue
 110  continue
c     
      return
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cvecinit( nn,b )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Filling zero to the vector 'g'.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nn,i
      complex*16 b(*)
c     
      do 100 i=1,nn
         b(i) = dcmplx( 0.d0 )
 100  continue
c     
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpolate( ncomp,nderiv,rsta,rrsta,g,u )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ncomp,nderiv
      real*8 rsta,rrsta(3)
      complex*16 g(3*ncomp),u(ncomp)
c     
      real*8 dh(3)
c     
      integer ip(3),ier,i,itmp,icomp
      complex*16 a(3,3),b(3),wk(3)
      real*8 eps
c     
      data eps / -1.d0 /
c     
      do 210 icomp=1,ncomp
         u(icomp) = dcmplx(0.d0)
 210  continue
c     
      do 100 i=1,3
         dh(i) = rrsta(i) - rsta
 100  continue
c     
      if( (dh(2).eq.0.d0).and.(nderiv.eq.0)) then
         itmp = ncomp + 1
         do 110 icomp=1,ncomp
	    u(icomp) = g(itmp)
	    itmp = itmp + 1
 110     continue
         return
      endif
c     
      do 120 i=1,3
         a(1,i) = dcmplx( 1.d0 )
         a(2,i) = dcmplx( dh(i) )
         a(3,i) = dcmplx( dh(i) * dh(i) / 2.d0 )
 120  continue
c     
      call fillinpb(nderiv,b)
      call glu(a,3,3,b,eps,wk,ip,ier)
c     
c     
      do 130 icomp=1,ncomp
         do 140 i=1,3
            u(icomp) = u(icomp) + b(i) * g( ncomp * (i-1) + icomp )
 140     continue
 130  continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fillinpb( nderiv,b )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nderiv
      complex*16 b(3)
c     
      if( (nderiv.ne.0).and.(nderiv.ne.1).and.(nderiv.ne.2) )
     &     pause 'invalid argument (fillinpb)'
      if(nderiv.eq.0) then
         b(1) = dcmplx( 1.d0 )
         b(2) = dcmplx( 0.d0 )
         b(3) = dcmplx( 0.d0 )
      elseif(nderiv.eq.1) then
         b(1) = dcmplx( 0.d0 )
         b(2) = dcmplx( 1.d0 )
         b(3) = dcmplx( 0.d0 )
      elseif(nderiv.eq.2) then
         b(1) = dcmplx( 0.d0 )
         b(2) = dcmplx( 0.d0 )
         b(3) = dcmplx( 1.d0 )
      endif
c     
      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makeoutputfilempi(dir,obs,nr,nri,indexmap,outu,outlcd)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      character*80 dir,outu(*),outlcd(*)
      integer nr,nri
      integer indexmap(nri)
      character*80 tstr(2)
      character*80 tmpchar1,tmpchar4
      character*80 obs
      integer lstr(2),ir
      integer lchar1,lchar4
c
      tmpchar1 = adjustL(dir)
      lchar1 = len_trim(dir)
      tmpchar1 = tmpchar1(1:lchar1)//'XY'
      lchar1 = len_trim(tmpchar1)
c
      do ir=1,nri
         call makenumchar(indexmap(ir),nr,lstr(1),tstr(1))
         tmpchar4 = tmpchar1(1:lchar1)//tstr(1)(1:lstr(1))//'.'//obs
         lchar4 = len_trim(tmpchar4)
         outu(ir) = tmpchar4(1:lchar4)//'.UB...SH.spc'
         outlcd(ir) = tmpchar4(1:lchar4)//'.PB...SH.spc'
      enddo
c      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makeoutputfile(dir,obs,nr,outu,outlcd)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      character*80 dir,outu(*),outlcd(*)
      integer nr
      character*80 tstr(2)
      character*80 tmpchar1,tmpchar4
      character*80 obs
      integer lstr(2),ir
      integer lchar1,lchar4
c
      tmpchar1 = adjustL(dir)
      lchar1 = len_trim(dir)
      tmpchar1 = tmpchar1(1:lchar1)//'XY'
      lchar1 = len_trim(tmpchar1)
c
      do ir=1,nr
      call makenumchar(ir,nr,lstr(1),tstr(1))
      tmpchar4 = tmpchar1(1:lchar1)//tstr(1)(1:lstr(1))//'.'//obs
      lchar4 = len_trim(tmpchar4)
      outu(ir) = tmpchar4(1:lchar4)//'.UB...SH.spc'
      outlcd(ir) = tmpchar4(1:lchar4)//'.PB...SH.spc'
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makenumchar(i,n,length,str)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i,n,length
      character*80 str
c
      if(n.le.9) then
         write(str,'(i1.1)') i
         length = 1
      elseif(n.le.99) then
         write(str,'(i2.2)') i
         length = 2
      elseif(n.le.999) then
         write(str,'(i3.3)') i
         length = 3
      elseif(n.le.9999) then
         write(str,'(i4.4)') i
         length = 4
      elseif(n.le.99999) then
         write(str,'(i5.5)') i
         length = 5
      elseif(n.le.999999) then
         write(str,'(i6.6)') i
         length = 6
      endif
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer nzone,lsuf
      real*8 omega,vrmax(*),vsv(4,*)
c
      real*8 tvs,coef
      integer i
c
      tvs = 0.d0
      do 100 i=1,4
         if(i.eq.1) then
            coef = 1.d0
         else
            coef = coef 
         endif
         tvs = tvs + ( vsv(i,nzone) ) * coef
 100  continue
c     
      lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calamp(g,l,lsuf,maxamp,ismall,ratl)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer l,lsuf,ismall
      real*8 maxamp,ratl
      complex*16 g
c     
      real*8 amp,ampratio
c     
      ampratio = 0.d0
      amp = cdabs(g)
      if( amp.gt.maxamp ) maxamp = amp
      if ( (amp.ne.0.d0).and.(maxamp.ne.0.d0) ) then
         ampratio = amp / maxamp
      endif
      if( (ampratio.lt.ratl).and.(l.gt.lsuf) ) then
         ismall = ismall + 1
      else
         ismall = 0
      endif
c     
      return
      end
c
