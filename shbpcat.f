      program shbp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ************** shbp.f ****************
c Computation of back-propagated wavefield 
c including displacements and their locally Cartesian derivatives
c for SH synthetic seismograms in transversely isotropic media 
c                                                      2004.5 K.Kawai
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ----------------------------<<constants>>----------------------------
      implicit none
      real*8 pi
      integer maxnlay,maxlmax,ilog,iup
      integer maxnzone,maxnr,maxnstack
      parameter ( pi=3.1415926535897932d0 )
      parameter ( maxnlay = 8300 )
      parameter ( maxnzone = 15 )
      parameter ( maxnr = 400 )     ! number of xy-perturbations
      parameter ( maxnstack = 20 )  ! number of z-perturbations
      parameter ( maxlmax = 80000 )
      parameter ( ilog = 0 )
      parameter ( iup = 1)      ! 0-both; 1-onlyLCD
c ----------------------------<<variables>>----------------------------
c variable for the trial function
      integer nnlayer,nlayer(maxnzone)
      integer l,m
      real*8 ra(maxnlay+maxnzone+1),gra(3),plm(3,0:3,maxnr)
      complex*16 bvec(3,-2:2,maxnr)
      complex*16 bvecdt(3,-2:2,maxnr),bvecdp(3,-2:2,maxnr)
c variable for the structure
      integer nzone
      integer ndc,vnp
      real*8 rmin,rmax
      real*8 vrmin(maxnzone),vrmax(maxnzone)
      real*8 rrho(4,maxnzone),vsv(4,maxnzone),vsh(4,maxnzone)
      real*8 qmu(maxnzone)
      real*8 vra(maxnlay+2*maxnzone+1)
      real*8 rho(maxnlay+2*maxnzone+1)
      real*8 ecL(maxnlay+2*maxnzone+1)
      real*8 ecN(maxnlay+2*maxnzone+1)
c      real*8 gvra(3),grho(3),gecL(3),gecN(3)
      complex*16 coef(maxnzone)
c variable for the periodic range
      integer np,imin,imax
      real*8 tlen,omega,omegai
      complex*16 u(9,maxnstack,maxnr)
c locally Cartesian derivatives
      complex*16 ulcd(27,maxnstack,maxnr)
c variable for the source
c      integer spn,ns
c      real*8 spo,mu0
      real*8 eqlat,eqlon
c variable for the station
      integer nr,ir
      real*8 theta(maxnr),phi(maxnr)
      real*8 lat(maxnr),lon(maxnr)
c variable for the matrix elements
      complex*16 a0( 2,maxnlay+1 ), a2( 2,maxnlay+1 )
      complex*16  a( 2,maxnlay+1 )
      real*8 t( 4*maxnlay )
      real*8 h1( 4*maxnlay ),h2( 4*maxnlay )
      real*8 h3( 4*maxnlay ),h4( 4*maxnlay )
      complex*16 g( maxnlay+1 )
      complex*16 gtmp(1),gdertmp(1)
c variable for the file
      character*80 outu(maxnr),outlcd(maxnr)
      character*80 dir
      character*80 outfile
c variable for grid spacing
      complex*16 tmpc(maxnlay+1)
      real*8 gridpar(maxnzone),dzpar(maxnzone),vmin(maxnzone)
      real*8 re,ratc,ratl,maxamp
      integer kc,lsuf,ismall,llog
c variable for the stack point
      integer isp(maxnzone),jsp(maxnzone)
c variable for the output stack point
      integer nsta
      real*8 rsta(maxnstack),rrsta(3,maxnstack)
      integer iista(3,maxnstack),ciista
c variable for ths observer
      real*8 obslat,obslon
      character*80 obs
c other variables
      integer i,j,jj,nn,lda,ier,ista
      real*8 eps,work( 4*maxnlay ),lsq
      complex*16 dr(maxnlay+1),z(maxnlay+1)
      complex*16 cwork( 4*maxnlay )
c
      data lda/ 2 /
      data eps/ -1.d0 /
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *************** Inputting and computing the parameters ***************
c --- inputting parameter ---
      call pinput( maxnlay,maxnzone,maxnr,maxnstack,
     &     re,ratc,ratl,
     &     tlen,np,omegai,imin,imax,
     &     nzone,vrmin,vrmax,rrho,vsv,vsh,qmu,
     &     nr,theta,phi,lat,lon,dir,
     &     obslat,obslon,obs, nsta,rsta ) 
c --- computing the required parameters ---
c computing and checking the parameters
      rmin = vrmin(1)
      rmax = vrmax(nzone)
      ndc = nzone - 1
      do ir=1,nr
         theta(ir)= theta(ir) / 1.8d2 * pi
         phi(ir)= phi(ir)   / 1.8d2 * pi
      enddo
c ************************** Files Handling **************************
      call makeoutputfile(dir,obs,nr,outu,outlcd)
      do ir=1,nr
         if(iup.eq.0) then
            open( unit=11,file=outu(ir) )
            write(11,*) tlen
            write(11,*) np,nsta
            write(11,*) omegai,lat(ir),lon(ir)
c     write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
            write(11,*) obslat,obslon
            do ista=1,nsta
               write(11,*) rsta(ista)
            enddo
            close(11)
         endif
         open( unit=11,file=outlcd(ir) )
         write(11,*) tlen
         write(11,*) np,nsta,27
         write(11,*) omegai,lat(ir),lon(ir)
c         write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
         write(11,*) obslat,obslon
         do ista=1,nsta
            write(11,*) rsta(ista)
         enddo
         close(11)
      enddo
      if(ilog.eq.1) then
         open(unit=11,file='llog.log',status='unknown')
         close(11)
      endif
c ************************** Files Handling end **********************
c
c computing of the number and the location of grid points
      call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,
     &     imax,1,tlen,
     &     vmin,gridpar,dzpar )
      call calra ( maxnstack,maxnlay,maxnzone,
     &     nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,
     &     rmin,rmax,nlayer,ra,re,nsta,rsta,rrsta,iista,
     &     rmax,ciista )
c --- checking the parameter
      if ( nnlayer.gt.maxnlay )
     &     pause 'The number of grid points is too large.'
c computing the stack points
      call calsp( ndc,nlayer,isp,jsp )
c ******************* Computing the matrix elements *******************
c computing the structure grid points
      call calstg( nzone,rrho,vsv,vsh,
     &     nnlayer,nlayer,ra,rmax,
     &     vnp,vra,rho,ecL,ecN)
      do i=1,ndc+1
         call calmatc( nlayer(i),vnp,vra,rho,2,0,0,
     &        ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
         call calmatc( nlayer(i),vnp,vra,ecL ,2,1,1,
     &        ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
         call calmatc( nlayer(i),vnp,vra,ecL ,1,1,0,
     &        ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
         call calmatc( nlayer(i),vnp,vra,ecL ,0,0,0,
     &        ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
         call calmatc( nlayer(i),vnp,vra,ecN ,0,0,0,
     &        ra( isp(i) ),h4( jsp(i) ),work( jsp(i) ) )
         call caltl( nlayer(i),vnp,vra,rho,
     &        ra( isp(i) ),work( jsp(i) ) )
         call calt( nlayer(i),  t( jsp(i) ),  work( jsp(i) ),
     &        t( jsp(i) ) )
         call calhl( nlayer(i),vnp,vra,ecL,
     &        ra( isp(i) ),work( jsp(i) ) )
         call calt( nlayer(i), h3( jsp(i) ), work( jsp(i) ),
     &        h3( jsp(i) ) )
         call calhl( nlayer(i),vnp,vra,ecN,
     &        ra( isp(i) ),work( jsp(i) ) )
         call calt( nlayer(i), h4( jsp(i) ), work( jsp(i) ),
     &        h4( jsp(i) ) )
      enddo
c ******************** Computing the displacement *********************
      nn = nnlayer + 1
c
      llog = 0
      do i=imin,imax            ! f-loop start
         do ir=1,nr
            call cmatinit( 9,maxnstack,u(1,1,ir))
            call cmatinit( 27,maxnstack,ulcd(1,1,ir))
         enddo
         if ( i.ne.0 ) then
            omega = 2.d0 * pi * dble(i) / tlen
            call callsuf(omega,nzone,vrmax,vsv,lsuf)
            do ir=1,nr
               call matinit( 3,4,plm(1,0,ir) )
            enddo
c
            call calcoef( nzone,omega,qmu,coef )
c     
            call cmatinit( lda,nn,a0 )
            call cmatinit( lda,nn,a2 )
            do j=1,ndc+1
               call cala0( nlayer(j),omega,omegai,
     &              t(jsp(j)), h1(jsp(j)),
     &              h2(jsp(j)), h3(jsp(j)),
     &              h4(jsp(j)),
     &              coef(j), cwork(jsp(j)) )
               call overlap( nlayer(j),cwork(jsp(j)),
     &              a0( 1,isp(j) ) )
               call cala2( nlayer(j),h4(jsp(j)),
     &              coef(j), cwork(jsp(j)) )
               call overlap( nlayer(j),cwork(jsp(j)),
     &              a2( 1,isp(j) ) )
            enddo
c     
            kc = 1
            ismall = 0
            maxamp = -1.d0
            llog = maxlmax
            do l=0,maxlmax      ! l-loop start
               lsq = dsqrt( dble(l)*dble(l+1) )
               if( ismall.gt.20 ) then
                  if(llog.gt.l) llog = l
                  exit
               endif
c     
               do jj=1,maxnlay+1 ! initialize
                  tmpc(jj) = dcmplx(0.d0)
               enddo
c
c ***** Computing the trial function *****
               do ir=1,nr
                  call calbvec( l,theta(ir),phi(ir),
     &                 plm(1,0,ir),bvec(1,-2,ir),
     &                 bvecdt(1,-2,ir),bvecdp(1,-2,ir) )
               enddo
c computing the coefficient matrix elements
c --- Initializing the matrix elements
               call cmatinit( lda,nn,a )
               call cala( nn,l,lda,a0,a2,a )
c
               do m=-2,2        ! m-loop start
                  if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
c     ---f_r---
                     call cvecinit( nn,g )
                     call calgpfsh(l,m,1,g(nn))
                     if ( (m.eq.-2).or.(m.eq.-l) ) then
                        call dclisb0( a,nn,1,lda,g,eps,dr,z,ier)
                     else
                        call dcsbsub0( a,nn,1,lda,g,eps,dr,z,ier)
                     endif
                     do ir=1,nr
                        do ista=1,nsta
                           call interpolate(1,0,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),gtmp )
                           call interpolate(1,1,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),
     &                          gdertmp )
                           call calu( gtmp(1),lsq,bvec(1,m,ir),
     &                          u(1,ista,ir) )
                           call calulcd( gtmp(1),gdertmp(1),lsq,
     &                          rsta(ista),theta(ir),
     &                          bvec(1,m,ir),bvecdt(1,m,ir),
     &                          bvecdp(1,m,ir),ulcd(1,ista,ir) )
                        enddo
                     enddo
c     ---f_\theta---
                     call cvecinit( nn,g )
                     call calgpfsh(l,m,2,g(nn))
                     call dcsbsub0( a,nn,1,lda,g,eps,dr,z,ier)
                     do ir=1,nr
                        do ista=1,nsta
                           call interpolate(1,0,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),gtmp )
                           call interpolate(1,1,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),
     &                          gdertmp )
                           call calu( gtmp(1),lsq,bvec(1,m,ir),
     &                          u(4,ista,ir) )
                           call calulcd( gtmp(1),gdertmp(1),lsq,
     &                          rsta(ista),theta(ir),
     &                          bvec(1,m,ir),bvecdt(1,m,ir),
     &                          bvecdp(1,m,ir),ulcd(10,ista,ir) )
                        enddo
                     enddo
c     ---f_\phi---
                     call cvecinit( nn,g )
                     call calgpfsh(l,m,3,g(nn))
                     call dcsbsub0( a,nn,1,lda,g,eps,dr,z,ier)
                     call calamp(g(ciista),l,lsuf,maxamp,ismall,ratl)
                     do ir=1,nr
                        do ista=1,nsta
                           call interpolate(1,0,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),gtmp )
                           call interpolate(1,1,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),
     &                          gdertmp )
                           call calu( gtmp(1),lsq,bvec(1,m,ir),
     &                          u(7,ista,ir) )
                           call calulcd( gtmp(1),gdertmp(1),lsq,
     &                          rsta(ista),theta(ir),
     &                          bvec(1,m,ir),bvecdt(1,m,ir),
     &                          bvecdp(1,m,ir),ulcd(19,ista,ir) )
                        enddo
                     enddo
                  endif
               enddo            ! m-loop end
            enddo               ! l-loop end
         endif
c ************************** Files Handling **************************
         do ir=1,nr
            if(iup.eq.0) then
               open( unit=11,file=outu(ir),
     &              position='append',status='old')
            endif
            open( unit=12,file=outlcd(ir),
     &           position='append',status='old')
            do ista=1,nsta
               if(iup.eq.0) then
                  write(11,*) i,
     &                 dble(u(1,ista,ir)),dimag(u(1,ista,ir))
                  do jj=2,9
                     write(11,*) dble(u(jj,ista,ir)),
     &                    dimag(u(jj,ista,ir))
                  enddo
               endif
               write(12,*) i,
     &              dble(ulcd(1,ista,ir)),dimag(ulcd(1,ista,ir))
               do jj=2,27
                  write(12,*) dble(ulcd(jj,ista,ir)),
     &                 dimag(ulcd(jj,ista,ir))
               enddo
            enddo
            if(iup.eq.0) close(11)
            close(12)
         enddo
         if(ilog.eq.1) then
            open(unit=11,file='llog.log',position='append',status='old')
            write(11,*) i,llog,nnlayer
            close(11)
         endif
c************************** Files Handling end **********************
c
      enddo                     ! f-loop end
c     

      end
