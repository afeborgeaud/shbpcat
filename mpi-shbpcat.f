      program mpishbpcat
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ************** shbp.f ****************
c Computation of back-propagated wavefield 
c including displacements and their locally Cartesian derivatives
c for SH synthetic seismograms in transversely isotropic media 
c                                                      2004.5 K.Kawai
c parallel MPI 
c        2012.03 K.Konishi
c
c Catalogue (epicentral distance grid) for efficient
c computation of partial derivatives by exploiting
c the spherical symmetry of the input model.
c Modification of the parallelization
c                               2018.9 A.F.E.Borgeaud
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ----------------------------<<constants>>----------------------------
      implicit none
      real*8 pi
      integer maxnlay,maxlmax,ilog,iup
      integer maxnzone,maxnr,maxnstack
      parameter ( pi=3.1415926535897932d0 )
      parameter ( maxnlay = 8300 )
      parameter ( maxnzone = 35 )
      parameter ( maxnr = 15000 )     ! number of xy-perturbations
      parameter ( maxnstack = 150 )  ! number of z-perturbations
      parameter ( maxlmax = 80000 )
      parameter ( ilog = 0)
      parameter ( iup = 1)      ! 0-both; 1-onlyLCD
c ----------------------------<<variables>>----------------------------
c variable for the trial function
      integer nnlayer,nlayer(maxnzone)
      integer l,m,mm
      real*8 ra(maxnlay+maxnzone+1),gra(3)
      real*8, allocatable, dimension(:,:,:) :: plm
      complex*16, allocatable, dimension(:,:,:) ::  bvec
      complex*16, allocatable, dimension(:,:,:) :: bvecdt, bvecdp
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
      integer np,imin,imax,ns
      real*8 tlen,omega,omegai
c non zero u components for SH case ignoring dependence on \phi
      complex*16, allocatable, dimension(:,:,:,:) :: u
c non zero locally cartesian derivatives for the SH case ignoring dependence on \phi
      complex*16, allocatable, dimension(:,:,:,:) :: ushlcd

c variable for the source
c      integer spn,ns
c      real*8 spo,mu0
      real*8 eqlat,eqlon
c variable for the station
      integer nr,ir,nri
      integer, allocatable, dimension(:) :: indexmap
      real*8 theta0(maxnr),phi0(maxnr)
      real*8 lat0(maxnr),lon0(maxnr)
      real*8, allocatable, dimension(:) :: theta,phi,lat,lon
c variable for the matrix elements
      complex*16 a0( 2,maxnlay+1 ), a2( 2,maxnlay+1 )
      complex*16  a( 2,maxnlay+1 )
      real*8 t( 4*maxnlay )
      real*8 h1( 4*maxnlay ),h2( 4*maxnlay )
      real*8 h3( 4*maxnlay ),h4( 4*maxnlay )
      complex*16 g( maxnlay+1 )
      complex*16 gtmp(1),gdertmp(1)
c variable for the file
      character*80, allocatable, dimension(:) :: outu,outlcd
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
c
      real*8 t1i,t1f,t2i,t2f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **************MPI***********************************
      include 'mpif.h'
      integer :: mpii
      integer :: petot, my_rank, ierr
      integer :: filenum, mpios
      integer :: outputmemory   ! MB
      integer :: outputinterval
      real(8) :: memoryperomega ! MB
c     memoryperomega = (9+27)*16*maxnstack * maxnr*0.000001
      
      integer :: outputindex
      integer, allocatable, dimension (:) :: outputi
      complex*16, allocatable, dimension (:,:,:,:,:) :: outputu
      complex*16, allocatable, dimension(:,:,:,:,:) :: outputulcd
c       when the values to be output use memory over outputmemory MB,
c       they are written in output files. The interval is outputinterval
c       memoryperomega is the quantity of memory used for one omega step
      character *2 :: char_rank
      data outputmemory /2000/
      call mpi_init (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
c     write(*,*) "myrank",my_rank
c     complex*16 u(9,maxnstack,maxnr)

c      complex*16 ulcd(27,maxnstack,maxnr) ! locally Cartesian derivatives
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      filenum = 10+my_rank
      write (char_rank, '(I2.2)') my_rank
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *************** Inputting and computing the parameters ***************
      if (my_rank.eq.0) then
c --- inputting parameter ---
         call pinput( maxnlay,maxnzone,maxnr,maxnstack,
     &        re,ratc,ratl,
     &        tlen,np,omegai,imin,imax,
     &        nzone,vrmin,vrmax,rrho,vsv,vsh,qmu,
     &        nr,theta0,phi0,lat0,lon0,dir,
     &        obslat,obslon,obs, nsta,rsta ) 
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call MPI_BCAST (re, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,
     $     ierr)
      call MPI_BCAST (ratc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,
     $     ierr)
      call MPI_BCAST (ratl, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,
     $     ierr)
      call MPI_BCAST (tlen, 1, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (np, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (omegai, 1, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (imin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (imax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (nzone, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (vrmin, maxnzone, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (vrmax, maxnzone, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (rrho, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (vsv, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (vsh, 4*maxnzone, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (qmu, maxnzone, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
c     
      call MPI_BCAST (obslat, 1, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (obslon, 1, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (nr, 1, MPI_INTEGER, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (theta0, maxnr, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (phi0, maxnr, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (lat0(1), maxnr, MPI_DOUBLE_PRECISION,
     $     0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (lon0(1), maxnr, MPI_DOUBLE_PRECISION,
     $     0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (dir, 80, MPI_CHARACTER,
     $     0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (obs, 80, MPI_CHARACTER,
     $     0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (rsta, maxnstack, MPI_DOUBLE_PRECISION, 0,
     $     MPI_COMM_WORLD, ierr)
      call MPI_BCAST (nsta, 1, MPI_INTEGER, 0,
     $     MPI_COMM_WORLD, ierr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c------------ define and allocate local variables on each processor ------------
      if (my_rank.eq.PETOT-1) then
            nri=nr - (PETOT-1) * (nr/PETOT)
      else
            nri=nr/PETOT
      endif
c      write(*,*) my_rank,nri

      allocate(indexmap(nri), theta(nri), phi(nri),
     & lat(nri), lon(nri))

      do ir =1,nri
            indexmap(ir)=my_rank*(nr/PETOT) + ir
c            write(*,*) my_rank,ir,indexmap(ir)
      enddo

      do ir =1,nri
            theta(ir)=theta0(indexmap(ir)) / 1.8d2 * pi
            phi(ir)=phi0(indexmap(ir)) / 1.8d2 * pi
            lat(ir)=lat0(indexmap(ir))
            lon(ir)=lon0(indexmap(ir))

c            write(*,*) my_rank,ir,indexmap(ir),theta(ir)
      enddo

      allocate(plm(3,0:3,nri))
      allocate(bvec(3,-2:2,nri), bvecdt(3,-2:2,nri),
     & bvecdp(3,-2:2,nri))
      allocate(u(6,-1:1,nsta,nri))
      allocate(ushlcd(18,-1:1,nsta,nri))
      allocate(outu(nri), outlcd(nri))

c
      memoryperomega = 9*3*16*nsta*nri*0.000001
      outputinterval = outputmemory/memoryperomega != integer*nr
c
      if(outputinterval.gt.(np+1)) outputinterval=np+1
c
      allocate (outputi(outputinterval))
      allocate (outputulcd(18,-1:1,nsta,nri,outputinterval))
      if (iup.eq.0)then
         memoryperomega = (3*18+9)*16*nsta*nri*0.000001
         outputinterval = outputmemory/memoryperomega
         allocate(outputu(9,-1:1,nsta,nri,outputinterval))
      endif
c
      write(*,*) "outputinterval ",outputinterval
      write(*,*) "allocate outputulcd ",18*3*nsta*nri*outputinterval*
     &  16*1e-6,"Mb"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- computing the required parameters ---
c     computing and checking the parameters
      rmin = vrmin(1)
      rmax = vrmax(nzone)
      ndc = nzone - 1
c
c ************************** Files Handling **************************
      call makeoutputfilempi(dir,obs,nr,nri,indexmap,outu,outlcd)
c     *** mimiwo sumashite mewokorasebahora tobiraga hirakeba subetegamieruwo
c
      do ir=1,nri
            if(iup.eq.0) then
c               call unlink(outu(ir))
               open( unit=11,file=trim(outu(ir)) )
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
c            call unlink(outlcd(ir))

            open( unit=filenum,file=trim(outlcd(ir)),form='binary',
     &         convert='big_endian' )
            write(filenum) tlen
c 8 used to identify PBCAT in Kibrary
            write(filenum) np,nsta,8
c            write(11) omegai,lat(ir),lon(ir)
            write(filenum) omegai,theta(ir)*1.8d2/pi,
     &        phi(ir)*1.8d2/pi
c phi(ir) is actually not used. Should be deleted in a future (cleaner) release...
c     write(11,*) theta(ir)*1.8d2/pi,phi(ir)*1.8d2/pi
            write(filenum) obslat,obslon
            do ista=1,nsta
               write(filenum) rsta(ista)
            enddo
            close(filenum)
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
cccccccccccccccccccccccc MPI ccccccccccccccccccccccccccccccccc
      outputindex =1
c
      t2i=MPI_Wtime()
      do i =0,np ! f-loop start
c         write(*,*) my_rank,i
	   do ir =1,nri
	     do ista =1,nsta
	       do mm =-1,1
		   do jj =1,6
		     u(jj,mm,ista,ir)=dcmplx(0.d0)
		   enddo
		 enddo
	     enddo
         enddo
c
	   do ir =1,nri
	     do ista =1,nsta
	       do mm =-1,1
		   do jj =1,18
		     ushlcd(jj,mm,ista,ir)=dcmplx(0.d0)
		   enddo
	       enddo
	     enddo
	   enddo
c
         if ( i.ne.0 ) then
            omega = 2.d0 * pi * dble(i) / tlen
            call callsuf(omega,nzone,vrmax,vsv,lsuf)
            do ir=1,nri
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
c               write(*,*) my_rank,l
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
               t1i = MPI_Wtime()
               do ir=1,nri
                  call calbvecreal2( l,theta(ir),
     &                 plm(1,0,ir),bvec(1,-2,ir),
     &                 bvecdt(1,-2,ir),bvecdp(1,-2,ir) )
               enddo
c               t1f = MPI_Wtime()
c               write(*,'("calbvec in ",e20.11," s")') t1f-t1i
c computing the coefficient matrix elements
c --- Initializing the matrix elements
               call cmatinit( lda,nn,a )
               call cala( nn,l,lda,a0,a2,a )
c
c---- for the SH case, it seems only m=1 is necessary (m=-1 is calculated in a later step) ---
               do m=-1,1        ! m-loop start
                  if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then

c                        call cvecinit( nn,g )

c---- we do not compute f_r since it is always 0 for the SH case -------
c
c---- f_\theta -------------------------------------------------------------
                     call cvecinit( nn,g )
                     call calgpfsh(l,m,2,g(nn))
                     if((m.eq.-1).or.(m.eq.-l)) then
                        call dclisb0( a,nn,1,lda,g,eps,dr,z,ier)
                     else
                        call dcsbsub0( a,nn,1,lda,g,eps,dr,z,ier)
                     endif
                     do ista=1,nsta
                           call interpolate(1,0,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),gtmp )
                           call interpolate(1,1,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),
     &                          gdertmp )
                     do ir=1,nri
                           call calu( gtmp(1),lsq,bvec(1,m,ir),
     &                          u(1,m,ista,ir) )
                           call calushlcd( gtmp(1),gdertmp(1),lsq,
     &                          rsta(ista),theta(ir),
     &                          bvec(1,m,ir),bvecdt(1,m,ir),
     &                          bvecdp(1,m,ir),ushlcd(1,m,ista,ir) )
                        enddo
                     enddo
c
c      mm=m
c      if (m.eq.mm) then
c      write(*,*) "f_\theta",i,l,mm
c      write(*,'(e20.11,e20.11,e20.11,e20.11,e20.11,e20.11' //
c     & ',e20.11,e20.11,e20.11)')
c     &      dble(ushlcd(1,mm,1,1)),dble(ushlcd(2,mm,1,1)),
c     &      dble(ushlcd(3,mm,1,1)),dble(ushlcd(4,mm,1,1)),
c     &      dble(ushlcd(5,mm,1,1)),dble(ushlcd(6,mm,1,1)),
c     &      dble(ushlcd(7,mm,1,1)),dble(ushlcd(8,mm,1,1)),
c     &      dble(ushlcd(9,mm,1,1))
c      write(*,'(e20.11,e20.11,e20.11,e20.11,e20.11,e20.11' //
c     & ',e20.11,e20.11,e20.11)')
c     &      dimag(ushlcd(1,mm,1,1)),dimag(ushlcd(2,mm,1,1)),
c     &      dimag(ushlcd(3,mm,1,1)),dimag(ushlcd(4,mm,1,1)),
c     &      dimag(ushlcd(5,mm,1,1)),dimag(ushlcd(6,mm,1,1)),
c     &      dimag(ushlcd(7,mm,1,1)),dimag(ushlcd(8,mm,1,1)),
c     &      dimag(ushlcd(9,mm,1,1))
c      endif

c---- f_\phi ----------------------------------------------------------------
c                     if (m.ge.0) then
                     t1i=MPI_Wtime()
                     call cvecinit( nn,g )
                     call calgpfsh(l,m,3,g(nn))
                     call dcsbsub0(a,nn,1,lda,g,eps,dr,z,ier)
                     if (m.eq.1) then
                        call calamp(g(ciista),l,lsuf,maxamp,ismall,ratl)
c                        write(*,*) "(amp, maxamp) =", i,l,m,
c     &                    g(ciista),maxamp
                     endif
                     do ista=1,nsta
                           call interpolate(1,0,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),gtmp )
                           call interpolate(1,1,rsta(ista),
     &                          rrsta(1,ista),g(iista(1,ista)),
     &                          gdertmp )
c                           write(*,*) gtmp,gdertmp
                         do ir=1,nri
c                           if(m.ge.0) then
                               call calu( gtmp(1),lsq,bvec(1,m,ir),
     &                           u(4,m,ista,ir) )
                               call calushlcd( gtmp(1),gdertmp(1),lsq,
     &                           rsta(ista),theta(ir),
     &                           bvec(1,m,ir),bvecdt(1,m,ir),
     &                           bvecdp(1,m,ir),ushlcd(10,m,ista,ir) )
c                           endif
                        enddo
                     enddo
c                     endif
c
c      mm=m
c      if (m.eq.mm) then
c      write(*,*) "f_\phi",i,l,mm
c      write(*,'(e20.11,e20.11,e20.11,e20.11,e20.11,e20.11' //
c     & ',e20.11,e20.11,e20.11)')
c     &      dble(ushlcd(10,mm,1,1)),dble(ushlcd(11,mm,1,1)),
c     &      dble(ushlcd(12,mm,1,1)),dble(ushlcd(13,mm,1,1)),
c     &      dble(ushlcd(14,mm,1,1)),dble(ushlcd(14,mm,1,1)),
c     &      dble(ushlcd(16,mm,1,1)),dble(ushlcd(17,mm,1,1)),
c     &      dble(ushlcd(18,mm,1,1))
c      write(*,'(e20.11,e20.11,e20.11,e20.11,e20.11,e20.11' //
c     & ',e20.11,e20.11,e20.11)')
c     &      dimag(ushlcd(10,mm,1,1)),dimag(ushlcd(11,mm,1,1)),
c     &      dimag(ushlcd(12,mm,1,1)),dimag(ushlcd(13,mm,1,1)),
c     &      dimag(ushlcd(14,mm,1,1)),dimag(ushlcd(14,mm,1,1)),
c     &      dimag(ushlcd(16,mm,1,1)),dimag(ushlcd(17,mm,1,1)),
c     &      dimag(ushlcd(18,mm,1,1))
c      endif
c
c                     t1f = MPI_Wtime()
c                     write(*,'("interp + calu in ",e20.11," s")') t1f-t1i
                  endif
               enddo            ! m-loop end
            enddo               ! l-loop end
c            t2f=MPI_Wtime()
c            write(*,'(i5," i= ",i5," took ",e20.11," s")')
c     &  my_rank,i,t2f-t2i
         endif
c ************************** Files Handling **************************
         outputi(outputindex)=i
         do ir =1,nri
            do ista =1,nsta
		   do mm =-1,1,2
		     do jj= 1,18
                   outputulcd(jj,mm,ista,ir,outputindex) =
     &               ushlcd(jj,mm,ista,ir)
                 enddo
               enddo
            enddo
         enddo
         if (iup .eq.0) then
            do ir =1,nri
               do ista =1,nsta
		     do mm =-1,1,2
                   do jj= 1,6
                     outputu(jj,mm,ista,ir,outputindex)=
     &                 u(jj,mm,ista,ir)
                   enddo
		     enddo
               enddo
            enddo
         endif
         
         if (outputindex.ge.outputinterval .or.
     $        i.eq.np) then
c            write(*,*) my_rank, "kakikomimasu"
            
            do ir=1,nri
ccccccccccccccccccccccccccccc
c               if(iup.eq.0) then
c                  open( unit=filenum+30,file=trim(outu(ir)),
c     $                 position='append',status='unknown', 
c     $                 share ='denyrw')
c               endif
cccccccccccccccccccccccccccccc

c120            continue
               open( unit=filenum,file=trim(outlcd(ir)),
     &              position='append',status='unknown',
     &                  share='denyrw',form='binary',
     &                 convert='big_endian',iostat=mpios)
c                goto 555
c 100           continue
      if(mpios.ne.0) then
        write(*,*) "Error while opening the file ",
     &     trim(outlcd(ir)), ", io status = ", mpios
        call mpi_finalize(ierr)
        stop
      endif
c write(*,*) my_rank,"waiting for 100s", outlcd(ir)
c TODO check if sleep is necessary
c               call system(" sleep 100")
c               goto 120
c555             continue
c                write(*,*) my_rank, "starts writing",outlcd(ir)
               do mpii=1, outputindex
                  do ista=1,nsta
cccccccccccccccccccccccccccccccc
		      if(iup.eq.0) then
			  write(filenum+30,*) outputi(mpii),
     &                       dble(outputu(1,-1,ista,ir, mpii)),
     &                       dimag(outputu(1,-1,ista,ir, mpii))
			  write(filenum+30,*)
     &                 dble(outputu(1,1,ista,ir, mpii)),
     &                       dimag(outputu(1,1,ista,ir, mpii))
			do mm=-1,1,2
			  do jj=2,6
				write(filenum+30,*)
     &			  dble(outputu(jj,mm,ista,ir
     &                          ,mpii)),dimag(outputu(jj,mm,ista,ir
     &                          ,mpii))
			  enddo
			  enddo
			endif
cccccccccccccccccccccccccccccccccccc
                  write(filenum) outputi(mpii),
     &                    dble(outputulcd(1,-1,ista,ir, mpii)),
     &                    dimag(outputulcd(1,-1,ista,ir, mpii))
                  write(filenum) dble(outputulcd(1,1,ista,ir, mpii)),
     &                    dimag(outputulcd(1,1,ista,ir, mpii))
c jj loop before mm loop is not really optimum, but DO NOT CHANGE IT OR KIBRARY WON'T WORK
              do jj=2,18
			  do mm=-1,1,2
                         write(filenum) dble(outputulcd(jj,mm,ista,ir,
     &                       mpii)), dimag(outputulcd(jj,mm,ista,ir,
     &                       mpii))
                       enddo
                     enddo
                  enddo
                  
               enddo
ccccccccccccccccccccccccccccc
c               if(iup.eq.0) close(filenum+30)
cccccccccccccccccccccccccccccc
            close(filenum)
            enddo 
            outputindex=0
         endif

         outputindex = outputindex+1       
         if(ilog.eq.1) then
            open(unit=11,file='llog.log',position='append',status='old')
            write(11,*) my_rank,i,llog,nnlayer
            close(11)
         endif
c**************************Files Handling end **********************
c
      enddo                     ! f-loop end
c     
c       stop      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      t2f=MPI_Wtime()
      write(*,*) my_rank,"finished in",t2f-t2i,"s"
c
      write(*,*) my_rank, "Ivalice looks to the horizon!"
c      call mpi_barrier(MPI_COMM_WORLD,ierr)
      
      call mpi_finalize(ierr)
      stop
      
      end program mpishbpcat
      

