c     Fortran library subroutine entry point

      SUBROUTINE cdmlib(
c     input file cdm.in
     $     lambda,beam,object,trope,
     $     materiaumulti,nnnr,tolinit,methodeit,polarizability,
     $     nquad,nlecture,filereread,nmat,
c     output file cdm.out
     $     nlocal,nmacro,nsection,nsectionsca,nquickdiffracte,nrig,
     $     nforce,nforced,ntorque,ntorqued,nproche,nlentille,nquicklens,
     $     nenergie,nobjet,
c     cube, sphere (includes multiple)
     $     density,side, sidex, sidey, sidez, hauteur,
     $     numberobjet, rayonmulti, xgmulti, ygmulti, zgmulti,
     $     epsmulti, epsanimulti,lc,hc,ng,
c     ellipsoid+arbitrary
     $     demiaxea,demiaxeb,demiaxec,thetaobj,phiobj,psiobj,
     $     namefileobj,
c     planewavecircular.in / planewavelinear.in files
     $     theta, phi, pp, ss, P0, w0, xgaus, ygaus, zgaus,namefileinc,
c     ondeplane multiple
     $     thetam, phim, ppm, ssm,E0m,nbinc,
c     return info stringf
     $     infostr, nstop,
c     return scalar results
     $     nbsphere, ndipole, aretecube,
     $     lambda10n, k0, tol1, ncompte, nloop,
     $     efficacite,efficaciteref,efficacitetrans,
     $     Cext, Cabs, Csca, Cscai, gasym, irra, E0,
     $     forcet, forcem,
     $     couplet, couplem,
     $     nxm, nym, nzm,
     $     incidentfield, localfield, macroscopicfield,
     $     xs, ys, zs, xswf, yswf, zswf,
     $     ntheta, nphi, thetafield, phifield, poyntingfield,
     $     forcex,forcey,forcez,forcexmulti,forceymulti,forcezmulti,
     $     torquex,torquey,torquez,torquexmulti,torqueymulti,
     $     torquezmulti,
     $     incidentfieldx, incidentfieldy, incidentfieldz,
     $     localfieldx, localfieldy, localfieldz,
     $     macroscopicfieldx, macroscopicfieldy, macroscopicfieldz,
     $     polarisa,epsilon,
     $     nfft2d,Eimagex,Eimagey,Eimagez,Eimageincx,Eimageincy,
     $     Eimageincz,Efourierx,Efouriery,Efourierz,Efourierincx,
     $     Efourierincy,Efourierincz,kxy,xy,numaper,gross,
c****************************************************
c     tableaux utilises que dans cdmlib
c****************************************************
c     taille double complex (3*nxm*nym*nzm)
     $     FF,FF0,FFloc,xr,xi,
c     taille double complex (3*nxm*nym*nzm,12)
     $     wrk,
c     taille double complex (8*nxm*nym*nzm)
     $     FFTTENSORxx, FFTTENSORxy,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,
     $     FFTTENSORzz,vectx,vecty,vectz,
c     taille double complex (nfft2d,nfft2d,3)
     $     Ediffkzpos,Ediffkzneg,
c     taille entier (nxm*nym*nzm)
     $     Tabdip,Tabmulti)

      implicit none
      
c     integer
      integer ii,jj,kk,ll,i,j,k,l,cnt,nstop,cntwf,kkm,jjm,iim
      integer  nlocal,nmacro,nsection,nsectionsca,nforce ,nforced
     $     ,nquickdiffracte,ntorque,ntorqued,nsens,nproche,nprochefft
     $     ,nlecture,nlecture1,long,long1,ierror,nrig,nquicklens,nquad
     $     ,nenergie,nmat

c     variables for the object
      integer nbsphere3,nbsphere,ndipole,test ,numberobjet,is,ng
      integer nx,ny,nz,nx2,ny2,nxy2,nz2,nxm,nym,nzm,ntotal,nxm2,nym2
     $     ,nzm2,nxym2
      integer subunit,nsubunit,comparaison
c     definition of the size for the code
      INTEGER nmax, ntotalm

c     variables for the positions
      double precision x,y,z,xmin,xmax,ymin,ymax ,zmin,zmax ,rayon,side
     $     ,sidex,sidey,sidez,hauteur,density ,xgmulti(numberobjet)
     $     ,ygmulti(numberobjet) ,zgmulti(numberobjet)
     $     ,rayonmulti(numberobjet),demiaxea ,demiaxeb,demiaxec,thetaobj
     $     ,phiobj,psiobj
      double precision aretecube
      integer iphi,itheta,nphi,ntheta
      DOUBLE PRECISION,DIMENSION(nxm*nym*nzm)::xs,ys,zs
      DOUBLE PRECISION,DIMENSION(nxm*nym*nzm)::xswf,yswf,zswf
      
      double precision pi,lambda,lambda10n,k0,k02,k03,epi,epr,c

c     variables for the material
      double precision eps0,quatpieps0
      double complex, dimension(nxm*nym*nzm,3,3) :: polarisa,epsilon
      double complex eps,epsani(3,3),epsmulti(numberobjet)
     $     ,epsanimulti(3,3,numberobjet)
      character(2) polarizability
      character (64), DIMENSION(numberobjet) :: materiaumulti
      character(64) materiau,object,beam,namefileobj,namefileinc
     $     ,filereread,filereread1
      character(3) trope,file1
c     variables for the incident field and local field
      DOUBLE PRECISION, DIMENSION(nxm*nym*nzm) :: incidentfield,
     $     localfield,macroscopicfield,forcex,forcey,forcez, torquex
     $     ,torquey,torquez
      double precision forcexmulti(numberobjet)
     $     ,forceymulti(numberobjet),forcezmulti(numberobjet)
     $     ,torquexmulti(numberobjet),torqueymulti(numberobjet)
     $     ,torquezmulti(numberobjet)
      integer nbinc
      double precision thetam(10), phim(10), ppm(10), ssm(10)
      double complex E0m(10)
      double precision ss,pp,theta,phi,I0,Emod,tmp
      double complex Eloc(3),Em(3),E0,uncomp,icomp,zzero,Emx,Emy,Emz
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldx
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldy
      double complex, dimension(nxm*nym*nzm) :: macroscopicfieldz
      double complex, dimension(nxm*nym*nzm) :: localfieldx
      double complex, dimension(nxm*nym*nzm) :: localfieldy
      double complex, dimension(nxm*nym*nzm) :: localfieldz
      double complex, dimension(nxm*nym*nzm) :: incidentfieldx
      double complex, dimension(nxm*nym*nzm) :: incidentfieldy
      double complex, dimension(nxm*nym*nzm) :: incidentfieldz
      double complex propaesplibre(3,3)
      double complex, dimension(3*nxm*nym*nzm) :: FF,FF0,FFloc

c     Green function
      integer, dimension(nxm*nym*nzm) :: Tabdip,Tabmulti
      integer indice,indicex,indicey
      double complex, dimension(8*nxm*nym*nzm) :: FFTTENSORxx,
     $     FFTTENSORxy,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,
     $     FFTTENSORzz,vectx,vecty,vectz

      double precision forcet(3),forcem,forcemie
      double precision couplet(3),couplem,xg,yg,zg,lc,hc,couplemie
      double complex Eder(3,3)
      
c     computation of the cross section
      integer imaxk0
      double precision Cext,normal(3),deltatheta,deltaphi,Csca,Cscai
     $     ,Cabs,gasym,thetas,phis,u1,u2,MIECEXT,MIECABS,MIECSCA,GSCA
      double complex ctmp,ctmp1
      
c     variables for the iterative method
      INTEGER ldabi, nlar
      integer nnnr,ncompte,nt,IP(3)
      integer NLIM,nou,nstat,nloop,STEPERR
      DOUBLE PRECISION  NORM,TOL,tolinit,tol1,tole
   
      double complex, dimension(3*nxm*nym*nzm) :: xr,xi
      double complex, dimension(3*nxm*nym*nzm,12) :: wrk
      
c     Poynting vector and energy conservation
      integer i2,j2,ii2,jj2
      double precision Poyntinginc,efficacite,efficaciteref
     $     ,efficacitetrans

c     Info string
      character(64) infostr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     nouvelle variable a passer en argument d'entree
c     power et diametre
      double precision P0,w0,xgaus,ygaus,zgaus,xdip,ydip,zdip
      character(12) methodeit

c     nouvelle variable de sortie Irra     
      double precision irra
      
c     nouvelle variable
      integer nloin
      
c     variable pour avoir l'image a travers la lentille
      integer nlentille,nobjet,nfft2d,nfft2d2
      double precision kx,ky,kz,deltakx,deltaky,numaper,deltax,gross

      DOUBLE PRECISION,DIMENSION(max((ntheta+1)*nphi,nfft2d
     $     *nfft2d))::thetafield,phifield,poyntingfield
      double complex Ediffkzpos(nfft2d,nfft2d,3),Ediffkzneg(nfft2d
     $     ,nfft2d,3)
      double precision kxy(nfft2d),xy(nfft2d)
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d),Eimageincx(nfft2d*nfft2d)
     $     ,Eimageincy(nfft2d *nfft2d) ,Eimageincz(nfft2d*nfft2d)
     $     ,Efourierx(nfft2d*nfft2d) ,Efouriery(nfft2d*nfft2d)
     $     ,Efourierz(nfft2d*nfft2d),Efourierincx(nfft2d*nfft2d)
     $     ,Efourierincy(nfft2d*nfft2d) ,Efourierincz(nfft2d*nfft2d)
      double precision t1,t2

c     variables pour le temps
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer,dimension(8) :: values,values2,valuesi,valuesf
      integer iret,omp_get_max_threads
      integer*8 planf,planb,plan2f,plan2b,planfn,planbn
      integer FFTW_FORWARD,FFTW_ESTIMATE,FFTW_BACKWARD

      call date_and_time(date,time,zone,valuesi)
      
      call dfftw_init_threads(iret)
      if (iret.eq.0) then
         write(*,*) 'iret',iret
         infostr='Unlikely error during thread initialization'
         nstop=1
         return
      endif
      CALL dfftw_plan_with_nthreads(omp_get_max_threads())
      
c     FF0 : champ incident
c     FF  : dipole
c     FFloc : champ local
      write(*,*) '******************************************'
      write(*,*) '*************** INPUT DATA ***************'
      write(*,*) '******************************************'
      write(*,*) 'Wavelength       :',lambda
      write(*,*) 'Beam             : ',beam
      write(*,*) 'Object           : ',object,'iso',trope
      write(*,*) 'Discretization   : ',nnnr
      write(*,*) 'Size of the box  : ',nxm,nym,nzm
      write(*,*) 'Iterative method : ',methodeit,'tolerance asked'
     $     ,tolinit
      write(*,*) 'Rigorous  approx : ',nrig
      write(*,*) 'Local field      : ',nlocal
      write(*,*) 'Macroscopic field: ',nmacro
      write(*,*) 'Cross section    : ',nsection,'Csca',nsectionsca
      write(*,*) 'Quick cross section : ',nquickdiffracte
      write(*,*) 'Emissivity       : ',nenergie
      write(*,*) 'Lens             : ',nlentille
      write(*,*) 'Quick Lens       : ',nquicklens
      write(*,*) 'Optical force    : ',nforce ,nforced
      write(*,*) 'Optical torque   : ',ntorque ,ntorqued
      write(*,*) 'Near field       : ',nproche
      write(*,*) 'Polarizability   : ',polarizability
      write(*,*) 'Quadrature Green : ',nquad
      write(*,*) 'Write mat file   : ',nmat


c     arret de suite si nnnr trop petit par rapport a n*m
      if (nnnr.gt.min(nxm,nym,nzm)) then
         infostr='nxm nym or nzm smaller than discretization'
         nstop = 1;
         return         
      endif
      
      
      materiau = materiaumulti(1)
c     Liste des arrets si le choix des options est mal faits.
c     si rien n'est demandé
      if (nobjet.eq.0.and.nlocal.eq.0.and.nmacro.eq.0
     $     .and.nsection.eq.0.and.nsectionsca.eq.0
     $     .and.nforce.eq.0.and.nforced.eq.0
     $     .and.ntorque.eq.0.and.ntorqued.eq.0
     $     .and.nlentille.eq.0.and.nenergie.eq.0) then 
         infostr='No calculation requested!'
         nstop = 1;
         return
      endif
c     liste des arrets pour emissivité
      if (nenergie.eq.1) then
         if (beam(1:7).eq.'antenna'.or.beam(1:9).eq.'arbitrary') then
            infostr='Can not compute emissivity with this beam'
            nstop=1
            return
         endif
         if (dabs(theta).ge.90.d0) then
            nstop=1
            infostr='illumination on the wrong side'
            return
         endif
         if (dabs(theta).gt.90.d0.and.dabs(theta).ge.45.d0) then
            nstop=1
            infostr='Angle too strong: Imprecise calculation'
            return
         endif
      endif
c     liste des arrets pour lentille
      if (nlentille.eq.1) then
         if (beam(1:7).eq.'antenna'.or.beam(1:9).eq.'arbitrary') then
            infostr='Can not compute total field with this beam'
            nstop=1
            return
         endif
         if (numaper.le.0.d0.or.numaper.gt.1.d0) then
            nstop=1
            infostr='problem with numerical aperture!'
            return
         endif

         if (dabs(gross).lt.1.d0) then
            nstop=1
            infostr='problem with manification!'
            return
         endif

         if (dabs(theta).ge.90.d0) then
            nstop=1
            infostr='illumination on the wrong side'
            return
         endif
         
      endif

c     liste des arrets pour section efficace

      if (nsectionsca.eq.1) then
         if (dabs(theta).ge.45.d0) then
            nstop=1
            infostr='Angle too strong: Imprecise calculation with FFT'
            return
         endif
      endif

      if (nsection.eq.1) then

         if (beam(1:11).ne.'pwavelinear'.and.beam(1:13).ne
     $        .'pwavecircular') then
            nstop=1
            infostr='Cross section only with plane wave'
            return
         endif
      endif




      
c     calculation size parameter initialization
      nmax = nxm*nym*nzm
      ntotalm = 8*nxm*nym*nzm
      ldabi = 3*nxm*nym*nzm
      nlar = 12
      nstop=0

      xg = xgmulti(1)
      yg = ygmulti(1)
      zg = zgmulti(1)
      rayon = rayonmulti(1)
      eps = epsmulti(1)
      epsani = epsanimulti(:,:,1)
c     pas assez discrétisé
      if (nmax.lt.8) then
         nstop=1
         infostr='Check Nxm Nym and Nzm!'
         return
      endif

      if (nnnr.lt.2) then
         nstop=1
         write(*,*) 'There is no discretization!'
         infostr='There is no discretization!'
         return
      endif
c     arret pour tolerance dans la méthode itérative
      if (tolinit.lt.1.d-12) then 
         infostr='Tolerance too small!'
         nstop = 1;
         return
      endif
      if (tolinit.gt.0.1d0) then 
         infostr='Tolerance too large!'
         nstop = 1;
         return
      endif

c     initialization of the data
      icomp=(0.d0,1.d0)
      uncomp=(1.d0,0.d0)
      pi=dacos(-1.d0)
      zzero=(0.d0,0.d0)
      c=299792458.d0
      eps0=1.d0
      quatpieps0=1.d0/(c*c*1.d-7)
      infostr='STARTED!'
      
c     transform the input data in meter
      lambda=lambda*1.d-9
      w0=w0*1.d-9

      if (w0.le.0.d0) then
         nstop=1
         infostr='waist=0!'
         return
      endif
      if (P0.le.0.d0) then
         nstop=1
         infostr='Power=0!'
         return
      endif
c     open the output file:
      open(99,file='output')
      write(99,*) '************* OUTPUT FILE ***************'

c     Intensity of the incident field
      open(36,file='incidentfieldx.mat')
      open(37,file='incidentfieldy.mat')
      open(38,file='incidentfieldz.mat')
      open(39,file='incidentfield.mat')
c     Intensity of the local field
      open(40,file='localfieldx.mat')
      open(41,file='localfieldy.mat')
      open(42,file='localfieldz.mat')
      open(43,file='localfield.mat')
c     Intensity of the macroscopic field
      open(44,file='macroscopicfieldx.mat')
      open(45,file='macroscopicfieldy.mat')
      open(46,file='macroscopicfieldz.mat')
      open(47,file='macroscopicfield.mat')
c     Intensity of the incident field
      open(136,file='incidentfieldxwf.mat')
      open(137,file='incidentfieldywf.mat')
      open(138,file='incidentfieldzwf.mat')
      open(139,file='incidentfieldwf.mat')
c     Intensity of the local field
      open(140,file='localfieldxwf.mat')
      open(141,file='localfieldywf.mat')
      open(142,file='localfieldzwf.mat')
      open(143,file='localfieldwf.mat')
c     Intensity of the macroscopic field
      open(144,file='macroscopicfieldxwf.mat')
      open(145,file='macroscopicfieldywf.mat')
      open(146,file='macroscopicfieldzwf.mat')
      open(147,file='macroscopicfieldwf.mat')
c     Far field discretization
      open(67,file='xwf.mat')
      open(68,file='ywf.mat')
      open(69,file='zwf.mat')
c     save the Poynting vecteur
      open(50,file='poynting.mat')
      open(51,file='theta.mat')
      open(52,file='phi.mat')
c     pour FFT
      open(53,file='poyntingpos.mat')
      open(54,file='poyntingneg.mat')
      open(55,file='kx.mat')
      open(56,file='ky.mat')

c     save the density of the optical force
      open(60,file='forcex.mat')
      open(61,file='forcey.mat')
      open(62,file='forcez.mat')
c     save the density of optical torque
      open(63,file='torquex.mat')
      open(64,file='torquey.mat')
      open(65,file='torquez.mat')

c     compute the relative permittivity fron a database versus the
c     wavelength of illumination
c     write(*,*) 'Relative permittivity',eps,materiau(1:2),lambda

      do k=1,numberobjet
         materiau = materiaumulti(k)
         if (materiau(1:2).ne.'xx') then
            
            call interpdielec(lambda,materiau,epr,epi,infostr,nstop)
            if (nstop.eq.1) return
            epsmulti(k)=(epr*uncomp+icomp*epi)
            write(99,*) 'Relative permittivity',eps
            write(99,*) 'Relative permittivity',eps,materiau(1:2),lambda
         else
            if (trope(1:3).eq.'iso') then
               write(*,*) 'Relative permittivity',eps
               write(99,*) 'Relative permittivity',eps
            else 
               do i=1,3
                  do j=1,3
                     write(*,*) 'Relative
     $permittivity',epsanimulti(i,j,k),i,j,k
                     write(99,*) 'Relative permittivity',epsani(i,j),i,j
                  enddo
               enddo
            endif
         endif
      enddo

c     test pour voir si la permittivite a un sens
      do k=1,numberobjet
         if (trope(1:3).eq.'iso') then
            if (dreal(epsmulti(k)).eq.1.d0 .and.
     $           dimag(epsmulti(k)).eq.0.d0) then
               infostr='Object is in vacuum'
               nstop=1
            elseif (cdabs(epsmulti(k)).eq.0.d0) then
               infostr='No relative permittivity given'
               nstop=1
            endif
         else
            tmp=0.d0
            do ii=1,3
               do jj=1,3
                  tmp=tmp+cdabs(epsani(ii,jj))
               enddo
            enddo
            if (tmp.eq.0.d0) then
               infostr='No relative permittivity given'
               nstop=1
            endif
c     Pour la pola tensorielle uniquement CR et CM disponible
            if (polarizability.ne.'RR'.and.polarizability.ne.'CM') then
               infostr
     $              ='Polarizability should be CM or RR with anisotropy'
               nstop=1
            endif
         endif
      enddo
      if (nstop.eq.1) return

      eps = epsmulti(1)
      materiau = materiaumulti(1)

c     wavenumber
      k0=2.d0*pi/lambda
      k03=k0*k0*k0
c     look  for compute near field with FFT
c     décrément de 1 de nproche pour faciliter le code C
      nproche=nproche-1
      write(*,*) 'number maximum of subunit :nmax = ',nmax     
      write(*,*) 'number of layer for the object :nnnr = ',nnnr

      write(*,*) '*************END INPUT DATA *******************'
      write(*,*) ' '
      
c     passe objet dans boite si utilise FFT
      if (nlecture.eq.1.and.nproche.eq.-1) nproche=0
      if (nenergie.eq.1.and.nproche.eq.-1) nproche=0
      if (nquickdiffracte.eq.1.and.nproche.eq.-1) nproche=0
      if (nquicklens.eq.1.and.nproche.eq.-1) nproche=0
      if((nrig.eq.6.or.nrig.eq.7).and.nproche.eq.-1) nproche=0
      if (beam(1:8).eq.'gfftwave'.and.nproche.eq.-1) nproche=0

      
      nprochefft=0
      if (nproche.ge.1) then
         nprochefft=nproche
         nproche=0
      endif

      write(*,*) '************** BEGIN OBJECT *******************'
      write(*,*) 'Built the object : ',object(1:20)
c     Built the object
      if (object(1:6).eq.'sphere') then
         numberobjet=1
         call objetsphere(trope,eps,epsani,eps0,xs,ys,zs,k0,aretecube
     $        ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,polarizability
     $        ,nproche,epsilon,polarisa,rayon,xg,yg,zg,nmat,infostr
     $        ,nstop)
         write(99,*) 'sphere',rayon
      elseif (object(1:12).eq.'inhomosphere') then
         numberobjet=1
         write(*,*) 'coucou4',nproche
         if (trope.ne.'iso') then
            nstop=1
            infostr='Permittivity not scalar for inhomogenous sphere'
            write(99,*)
     $           'Permittivity not scalar for inhomogenous sphere'
            return
         endif
         
         
         call objetsphereinhomo(eps,eps0,xs,ys,zs,k0 ,aretecube ,tabdip
     $        ,nnnr,nmax,nbsphere,ndipole,nx,ny,nz ,polarizability
     $        ,nproche,epsilon,polarisa,rayon,lc,hc ,ng
     $        ,macroscopicfieldx ,nmat,infostr ,nstop)

         write(99,*) 'inhomogenous sphere',rayon,lc,hc,ng

         macroscopicfieldx=0.d0
      elseif (object(1:12).eq.'inhomocuboid') then
         numberobjet=1
         write(*,*) 'coucou4'
         if (trope.ne.'iso') then
            nstop=1
            infostr='Permittivity not scalar for inhomogenous cuboid'
            write(99,*)
     $           'Permittivity not scalar for inhomogenous cuboid'
            return
         endif
         
         write(*,*) 'prout',sidex,sidey,sidez,xg,yg ,zg,lc,hc,ng
         call objetparainhomo(eps,eps0,xs,ys,zs,k0 ,aretecube ,tabdip
     $        ,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym,nzm
     $        ,polarizability ,epsilon ,polarisa,sidex,sidey,sidez,xg,yg
     $        ,zg,lc,hc,ng,macroscopicfieldx,nproche,nmat,infostr,nstop)

         write(99,*) 'inhomogenous cuboid',sidex,sidey,sidez,xg,yg ,zg
     $        ,lc,hc,ng

         macroscopicfieldx=0.d0
         
      elseif (object(1:12).eq.'randomsphere') then
         numberobjet=1

         call objetsphererandom(trope,eps,epsani,eps0,xs,ys,zs,k0
     $        ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm
     $        ,nym ,nzm,polarizability,epsilon,polarisa,sidex,sidey
     $        ,sidez,xg,yg,zg ,xswf,yswf,zswf,rayon,density,ng,nproche
     $        ,nmat,infostr ,nstop)
         
         
         write(99,*) 'random sphere',rayon,density,ng

         

         
      elseif (object(1:4).eq.'cube') then
         write(99,*) 'cube:side',side
         write(*,*) 'cube:side',side
         numberobjet=1
         call objetcube(trope,eps,epsani ,eps0,xs,ys,zs,k0,aretecube
     $        ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,polarizability
     $        ,epsilon,polarisa,side,xg,yg,zg,nmat,infostr,nstop)
      elseif(object(1:6).eq.'cuboid') then
         write(99,*) 'cuboid:sidex,sidey,sizez ',sidex,sidey,sidez
         write(*,*) 'cuboid:sidex,sidey,sizez ',sidex,sidey,sidez
         numberobjet=1
         call objetpara(trope,eps,epsani,eps0,xs,ys,zs,k0,aretecube
     $        ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym,nzm
     $        ,polarizability ,epsilon,polarisa,sidex,sidey,sidez,xg,yg
     $        ,zg,phiobj ,thetaobj,psiobj,nproche,nmat,infostr,nstop)
         write(99,*) 'side',sidex,sidey,sidez,nx,ny,nz
         write(*,*) trope,eps,epsani,eps0,k0,aretecube ,nnnr,nmax
     $        ,nbsphere,ndipole,nx,ny,nz,polarizability,sidex,sidey
     $        ,sidez,xg,yg,zg,phiobj ,thetaobj ,psiobj,nproche,infostr
     $        ,nstop
      elseif(object(1:9).eq.'ellipsoid') then
         numberobjet=1
         call objetellipse(trope,eps,epsani,eps0,xs,ys,zs,k0 ,aretecube
     $        ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym,nzm
     $        ,polarizability ,nproche,epsilon,polarisa,demiaxea
     $        ,demiaxeb,demiaxec,xg ,yg,zg,phiobj,thetaobj,psiobj
     $        ,nmat,infostr,nstop)
         write(99,*) 'ellipse',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:8).eq.'nspheres') then
         call objetnspheres(trope,epsmulti,epsanimulti,numberobjet
     $        ,numberobjet,xgmulti,ygmulti,zgmulti,rayonmulti,eps0,xs
     $        ,ys,zs,k0 ,aretecube,tabdip,tabmulti,nnnr,nmax,nbsphere
     $        ,ndipole,nx,ny,nz,polarizability,nproche ,epsilon,polarisa
     $        ,nmat,infostr,nstop)
         write(99,*) 'multisphere',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:8).eq.'cylinder') then
         numberobjet=1
         call objetcylindre(trope,eps,epsani,eps0,xs,ys,zs,k0 ,aretecube
     $        ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym,nzm
     $        ,polarizability,nproche,epsilon,polarisa,rayon ,hauteur,xg
     $        ,yg,zg,phiobj,thetaobj,psiobj,nmat,infostr,nstop)
         write(99,*) 'cylindre',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:16).eq.'concentricsphere') then
         call objetsphereconcentric(trope,epsmulti,epsanimulti
     $        ,numberobjet,numberobjet,xg,yg,zg,rayonmulti,eps0,xs,ys
     $        ,zs,k0,aretecube,tabdip,tabmulti,nnnr,nmax,nbsphere
     $        ,ndipole,nx,ny,nz,polarizability,nproche
     $        ,epsilon,polarisa,nmat,infostr,nstop)
         write(99,*) 'concentricsphere',nbsphere,ndipole,nx,ny,nz
      elseif(object(1:9).eq.'arbitrary') then
         numberobjet=1
         call objetarbitrary(trope,eps,epsani,eps0,xs,ys,zs,k0
     $        ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz
     $        ,polarizability,namefileobj,nproche,epsilon ,polarisa
     $        ,infostr,nstop)
         write(99,*) 'arbitrary',nbsphere,ndipole,nx,ny,nz
      else
         write(99,*) 'Object unknown'
         infostr='Object unknown!'
         nstop=1
         return
      endif
      write(*,*) 'Object under study created : ',object     
      write(*,*) 'Error (O right, 1 problem)',nstop
      
      
      if (nstop.eq.1) return
c     if (nnnr*nnnr*nnnr.gt.nmax) then
c     nstop=1
c     write(*,*) 'nxm nym and nzm too small'
c     infostr
c     $    ='nxm nym and nzm too small compare to the discretization'
c     return
c     endif

      
      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box : Box too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif
      if (nx.le.1.or.ny.le.1.or.nz.le.1) then
         nstop=1
         infostr='Dimension Problem nx ny nz too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif
      if (nstop.eq.1) return
      if (nstop == -1) then
         infostr = 'Calculation cancelled after object created!'
         return
      endif

c     ecriture dans fichiers du epsilon
      if (nmat.eq.0) then
         open(66,file='epsilon.mat')
         if (trope.eq.'iso') then
            do i=1,ndipole
               k=tabdip(i)
               if (k.ne.0) then
                  write(66,*) dreal(epsilon(k,2,2)),dimag(epsilon(k,2
     $                 ,2))
               else
                  write(66,*) 1.d0 , 0.d0
               endif
            enddo
         else
            do i=1,ndipole
               k=tabdip(i)
               if (k.ne.0) then
                  do ii=1,3
                     do jj=1,3
                        write(66,*) dreal(epsilon(k,ii,jj))
     $                       ,dimag(epsilon(k,ii,jj))
                     enddo
                  enddo
               else
                  do ii=1,3
                     do jj=1,3
                        write(66,*) 1.d0 , 0.d0
                     enddo
                  enddo
               endif
            enddo          
         endif
c     fin ecriture du epsilon
         close(66)
      endif


      
      lambda10n = lambda/10.d0/cdabs(cdsqrt(eps))
      if (aretecube.gt.2.5d0*lambda10n) then
         nstop=1
         infostr='meshsize larger than lambda/4!'
         return
      endif

      write(*,*) 'Number of subunit for the object:',nbsphere
      write(*,*) 'Number of subunit for the mesh  :',ndipole
      write(*,*) 'Size of the mesh size           :',aretecube
      write(*,*) 'number of subunit along x y and z',nx,ny,nz
      write(99,*) 'number of subunit for the object',nbsphere
      write(99,*) 'number of subunit for the mesh ',ndipole
      write(99,*) 'number of subunit along x y and z',nx,ny,nz
      write(99,*) 'mesh size',aretecube
      write(99,*) 'lambda/(10n)',lambda/10.d0/cdabs(cdsqrt(eps))

      write(*,*) '*************** END OBJECT *****************'
      write(*,*) ' '

      if (beam(1:11).eq.'pwavelinear') then
         write(99,*) 'Beam pwavelinear'
      elseif (beam(1:13).eq.'pwavecircular') then
         write(99,*) 'Beam pwavecircular'
      elseif (beam(1:7).eq.'antenna') then
         write(99,*) 'Beam antenna'
      elseif (beam(1:11).eq.'gwavelinear') then
         write(99,*) 'Beam gwavelinear'
      elseif (beam(1:13).eq.'gwavecircular') then
         write(99,*) 'Beam gwavecircular'
      elseif (beam(1:14).eq.'gfftwavelinear') then
         write(99,*) 'Beam gwavelinear fft'
      elseif (beam(1:16).eq.'gfftwavecircular') then
         write(99,*) 'Beam gwavecircular fft'
      elseif (beam(1:15).eq.'gparawavelinear') then
         write(99,*) 'Beam gparawavelinear'
      elseif (beam(1:17).eq.'gparawavecircular') then
         write(99,*) 'Beam gparawavecircular'
      elseif (beam(1:15).eq.'wavelinearmulti') then 
         write(99,*) 'Beam wavelinearmulti'
      elseif (beam(1:9).eq.'arbitrary') then
         write(99,*) 'Beam arbitrary'
      else
         write(99,*) 'Beam unknown'
         infostr='Beam unknown!'
         nstop=1
         return
      endif

c     cré le fichier de data pour connaitre les options pour matlab
      if (nmat.eq.0) then
         open(900,file='inputmatlab.mat')
         write(900,*) nproche
         write(900,*) nlocal
         write(900,*) nmacro
         write(900,*) nsection
         write(900,*) nsectionsca    
         write(900,*) nquickdiffracte
         write(900,*) nforce
         write(900,*) nforced
         write(900,*) ntorque
         write(900,*) ntorqued
         write(900,*) nlentille
         write(900,*) nquicklens
         write(900,*) nphi
         write(900,*) ntheta+1
         if (trope.eq.'iso') write(900,*) 0
         if (trope.eq.'ani') write(900,*) 1
         write(900,*) nfft2d
         write(900,*) k0
         write(900,*) numaper
         write(900,*) nprochefft
         write(900,*) nobjet
         close(900)
      endif

c     ne fait que l'objet
      if (nobjet.eq.1) then 
         infostr='Dipole calculation completed!'
         return
      endif

      
c     multiplication by a factor 2: Toeplitz matrix transformed in a
c     circulant matrix with a doble size.
      nbsphere3=3*nbsphere
      nx2=2*nx
      ny2=2*ny
      nz2=2*nz
      nxy2=nx2*ny2
      ntotal=8*nx*ny*nz      

c     calcul les plans pour la FFT3D
      FFTW_FORWARD=-1
      FFTW_BACKWARD=+1
      FFTW_ESTIMATE=64
      call dfftw_plan_dft_3d(planb, nx2,ny2,nz2, vectx,vectx
     $     ,FFTW_BACKWARD,FFTW_ESTIMATE)
      write(*,*) 'planb',planb
      call dfftw_plan_dft_3d(planf, nx2,ny2,nz2, vectx,vectx
     $     ,FFTW_FORWARD,FFTW_ESTIMATE)

c     calcul les plans pour la FFT2D
      call dfftw_plan_dft_2d(plan2b, nfft2d,nfft2d, vectx,vectx
     $     ,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(plan2f, nfft2d,nfft2d, vectx,vectx
     $     ,FFTW_FORWARD,FFTW_ESTIMATE)
      
c     if the computation asked is rigourous then compute the Green
c     function
      if (nrig.eq.0.or.nrig.eq.3.or.nrig.eq.4.or.nrig.eq.5) then
       

         write(*,*)
     $        '**** Begin computation of the Green function ****'
c     Compute the Green fonction in free space
         write(*,*) 'nquqd',nquad
         if (nquad.eq.0) then
            call greencalculfft(nx,ny,nz,nx2,ny2,nz2,nxy2,nxm,nym,nzm
     $           ,aretecube,k0,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $           ,FFTTENSORyy ,FFTTENSORyz,FFTTENSORzz,planb)
         else
            call greencalculquadfft(nx,ny,nz,nx2,ny2,nz2,nxy2,nxm,nym
     $           ,nzm,nquad,tolinit,aretecube,k0,FFTTENSORxx,FFTTENSORxy
     $           ,FFTTENSORxz ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz
     $           ,planb) 
         endif

         write(*,*) 'End FFT computation of the Green function'
         write(*,*) '*****************************************'
         write(*,*) ' '
c     Compute the incident field at each subunit of the object
         if (nstop == -1) then
            infostr = 'Calculation cancelled after FFT Green function!'
            return
         endif
      endif
      
      write(*,*) '********* Compute the incident field ********'
      write(*,*) 'Beam used',beam
      write(99,*) '******* Compute the incident field *******'
      write(99,*) 'Beam used',beam
      write(99,*) 'k0=',k0     
      if (beam(1:11).eq.'pwavelinear') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol1=',pp
         write(99,*) 'pol2=',ss
c     compute E0
         call irradiance(P0,w0,E0,irra)
         I0=cdabs(E0)**2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere           
            call ondeplane(xs(i),ys(i),zs(i),k0,E0,ss,pp,theta,phi,
     $           FF0(3*i-2),FF0(3*i-1),FF0(3*i),nstop,infostr)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL  
      elseif (beam(1:13).eq.'pwavecircular') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol=',ss
c     compute E0
         call irradiance(P0,w0,E0,irra)
         I0=cdabs(E0)**2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            call ondecirce(xs(i),ys(i),zs(i),k0,E0,ss,theta,phi,
     $           FF0(3*i-2),FF0(3*i-1),FF0(3*i))         
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL        
      elseif (beam(1:15).eq.'wavelinearmulti') then
         call irradiance(P0,w0,E0,irra)
         I0=cdabs(E0)**2
         tmp=0.d0
         do i=1,nbinc
            tmp=tmp+cdabs(E0m(i))**2.d0
            if (cdabs(E0m(i)).eq.0) then
               infostr='Magnitude equal to zero in incident beam'
               nstop=1
               return
            endif
         enddo
         tmp=dsqrt(cdabs(E0)**2.d0/tmp)
         do i=1,nbinc
            E0m(i)=E0m(i)*tmp
         enddo
         
         do i=1,nbinc
            write(99,*) 'theta=',thetam(i)
            write(99,*) 'phi=',phim(i)
            write(99,*) 'pol1=',ppm(i)
            write(99,*) 'pol2=',ssm(i)
            write(99,*) 'E0=',E0m(i)
         enddo
c     compute E0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,nbsphere
            call ondeplanemulti(xs(i),ys(i),zs(i),k0,E0m,ssm,ppm,thetam
     $           ,phim,nbinc,FF0(3*i-2),FF0(3*i-1),FF0(3*i),nstop
     $           ,infostr)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
      elseif (beam(1:7).eq.'antenna') then
         E0=dsqrt(3.d0*quatpieps0/(k0**4.d0)*P0)/quatpieps0*uncomp
         xdip=xgaus*1.d-9
         ydip=ygaus*1.d-9
         zdip=zgaus*1.d-9
c     write(*,*) 'coucou',xdip,ydip,zdip,E0
         call dipoleinside(xdip,ydip,zdip,xs,ys,zs,aretecube,nmax
     $        ,nbsphere)
c     write(*,*) 'coucou',xdip,ydip,zdip,theta,phi,E0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)             
         do i=1,nbsphere
            call dipoleinc(xdip,ydip,zdip,theta,phi,xs(i),ys(i),zs(i)
     $           ,aretecube,k0,E0,FF0(3*i-2),FF0(3*i-1),FF0(3*i),nstop
     $           ,infostr)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL             
c     compute the intenisty at the first object location
         call dipoleinc(xdip,ydip,zdip,theta,phi,xgmulti(1),ygmulti(1)
     $        ,zgmulti(1),aretecube,k0,E0,Em(1),Em(2),Em(3),nstop
     $        ,infostr)
         I0=cdabs(Em(1))**2.d0+cdabs(Em(2))**2.d0+cdabs(Em(3))**2.d0
         write(*,*) 'I0',I0
      elseif (beam(1:11).eq.'gwavelinear') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol1=',pp
         write(99,*) 'pol2=',ss

         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         
         nloin=0
         tol=tolinit
         call gaussianpuissance(P0,irra,w0,k0,E0)
         I0=cdabs(E0)**2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)    
         do i=1,nbsphere
            call gaussianchamp(xs(i),ys(i),zs(i),xgaus,ygaus,zgaus,theta
     $           ,phi,w0,k0,ss,pp,E0,FF0(3*i-2),FF0(3*i-1),FF0(3*i),tol
     $           ,nloin,nstop,infostr)
            
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
c         do i=1,nbsphere
c            write(*,*) 'champ',i,FF0(3*i-2),FF0(3*i-1),FF0(3*i)
c         enddo
c         stop
      elseif (beam(1:13).eq.'gwavecircular') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol=',ss
         
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         
         nloin=0
         tol=tolinit
         call gaussianpuissance(P0,irra,w0,k0,E0)
         I0=cdabs(E0)**2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)           
         do i=1,nbsphere
            call gaussianchampcirc(xs(i),ys(i),zs(i),xgaus,ygaus,zgaus
     $           ,theta,phi,w0,k0,ss,E0,FF0(3*i-2),FF0(3*i-1),FF0(3*i)
     $           ,tol,nloin,nstop,infostr)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      elseif (beam(1:8).eq.'gfftwave') then
         write(*,*) 'incident wave with FFT'
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol1=',pp
         write(99,*) 'pol2=',ss
         
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         
         call gaussianpuissance(P0,irra,w0,k0,E0)
         I0=cdabs(E0)**2
         call gaussianchampfft(xs,ys,zs,aretecube,k0,w0,E0,ss,pp ,theta
     $        ,phi,xgaus,ygaus,zgaus,beam,ndipole,nx,ny,nz,nxm,nym,nzm
     $        ,nmax,nfft2d,Eimagex,Eimagey,Eimagez,FF0,tolinit,plan2f
     $        ,plan2b,nstop ,infostr)
         
      elseif (beam(1:15).eq.'gparawavelinear') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol1=',pp
         write(99,*) 'pol2=',ss
         call gaussianpuissancepara(w0,P0,irra,E0)

         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         
         I0=cdabs(E0)**2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,nbsphere
            call gaussianparalinear(xs(i),ys(i),zs(i),xgaus,ygaus,zgaus
     $           ,theta,phi,w0,k0,ss,pp,E0,FF0(3*i-2),FF0(3*i-1),FF0(3
     $           *i),nstop,infostr)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL    
      elseif (beam(1:17).eq.'gparawavecircular') then
         write(99,*) 'theta=',theta
         write(99,*) 'phi=',phi
         write(99,*) 'pol=',ss
         
         xgaus=xgaus*1.d-9
         ygaus=ygaus*1.d-9
         zgaus=zgaus*1.d-9
         
         call gaussianpuissancepara(w0,P0,irra,E0)
         I0=cdabs(E0)**2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)          
         do i=1,nbsphere
            call gaussianparacirc(xs(i),ys(i),zs(i),xgaus,ygaus,zgaus
     $           ,theta,phi,w0,k0,ss,E0,FF0(3*i-2),FF0(3*i-1),FF0(3*i)
     $           ,nstop,infostr)           
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL             
      elseif  (beam(1:9).eq.'arbitrary') then
         call incidentarbitrary(xs,ys,zs,aretecube,FF0,nxm,nym,nzm
     $        ,nbsphere,nstop,namefileinc,infostr)
         if (nstop.eq.1) return
         
      endif

c     ecriture dans .mat et creation de incidentifeld et composantes.
      subunit=0
      if (nmat.eq.0) then
         do i=1,ndipole
            k=tabdip(i)
            if (k.ne.0) then
               subunit=subunit+1
               incidentfieldx(subunit)=FF0(3*k-2)
               incidentfieldy(subunit)=FF0(3*k-1)
               incidentfieldz(subunit)=FF0(3*k)
               incidentfield(subunit)=dsqrt(dreal(FF0(3*k-2)
     $              *dconjg(FF0(3*k-2))+FF0(3*k-1)*dconjg(FF0(3*k-1))
     $              +FF0(3*k)*dconjg(FF0(3*k))))
               write(36,*)dreal(incidentfieldx(subunit))
     $              ,dimag(incidentfieldx(subunit))
               write(37,*)dreal(incidentfieldy(subunit))
     $              ,dimag(incidentfieldy(subunit))
               write(38,*)dreal(incidentfieldz(subunit))
     $              ,dimag(incidentfieldz(subunit))
               write(39,*) incidentfield(subunit)
            else
               write(36,*) 0.d0,0.d0
               write(37,*) 0.d0,0.d0
               write(38,*) 0.d0,0.d0
               write(39,*) 0.d0
            endif            
         enddo
      else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC)         
         do i=1,ndipole
            k=tabdip(i)
            if (k.ne.0) then
               incidentfieldx(k)=FF0(3*k-2)
               incidentfieldy(k)=FF0(3*k-1)
               incidentfieldz(k)=FF0(3*k)
               incidentfield(k)=dsqrt(dreal(FF0(3*k-2) *dconjg(FF0(3*k
     $              -2))+FF0(3*k-1)*dconjg(FF0(3*k-1)) +FF0(3*k)
     $              *dconjg(FF0(3*k))))            
            endif            
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          
      endif
c     cLose intensity of the incident field
      close(36)
      close(37)
      close(38)
      close(39)
      
      write(99,*) 'Field modulus',cdabs(E0)
      if (nstop.eq.1) return
      if (nstop == -1) then
         infostr = 'Calculation cancelled after incident field!'
         return
      endif
      
      write(*,*) 'Power      : ',P0
      write(*,*) 'Irradiance : ',irra
      write(*,*) 'Intensity  : ',I0      
      write(*,*) '********** End incident field ***********'
      write(*,*) ' '

      
      if (nlecture.eq.1) then
         write(*,*) 'Reread or write the dipoles',nlecture
         call relecture(lambda,nfft2d,nrig,beam,object,trope,nnnr
     $        ,tolinit,density,side, sidex, sidey, sidez, hauteur,
     $        numberobjet,rayonmulti ,xgmulti, ygmulti, zgmulti,
     $        epsmulti,epsanimulti,lc,hc,ng ,demiaxea ,demiaxeb,demiaxec
     $        ,thetaobj,phiobj,psiobj ,namefileobj, theta, phi, pp, ss
     $        ,thetam,phim, ppm, ssm ,E0m,nbinc, P0, w0, xgaus, ygaus
     $        ,zgaus,namefileinc ,numberobjet,polarizability,nquad
     $        ,filereread,nlecture1 ,nstop ,infostr)
         write(*,*) 'reread',nlecture1,nstop,infostr
         if (nstop.eq.1) return
c     reread the local field         
         if (nlecture1.eq.1) then
            write(*,*) 'reread the dipoles'
            file1='.lf'
            long = len( trim(filereread  ) )
            long1 = len( trim( file1 ) )
            filereread1=filereread(1:long)//file1(1:long1)
            open(1000,file=filereread1,status='old',iostat=ierror,form
     $           ='unformatted')
            if (ierror.ne.0) then
               nstop=1
               infostr='file for read local field do not exist!'
               return
            else
               do i=1,nbsphere3
                  read(1000) FF(i),FFloc(i)
               enddo
            endif
            close(1000)
            goto 1000
         endif
      endif
      write(*,*) '**** Begin to solve the linear system ****'
      call cpu_time(t1)
      call date_and_time(date,time,zone,values)
      
      if (nrig.eq.0) then
         write(99,*) '***** Solve the linear system *****'
c     Compute the local field at each subunit position by solving Ax=b
c     as an initial guess the incident field
c     initilization for solve Ax=b
         tol=tolinit
         ncompte=0
         nloop=0

         write(*,*) 'nproche',nproche
         if (nproche.eq.-1) then
            write(*,*) 'rrrrrrrr'
            call inverserig(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $           ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $           ,Tabdip,ntotalm,ntotal,ldabi,nlar,nmax,ndipole,nxm,nym
     $           ,nzm,nx,ny,nz,nx2,ny2,nxy2,nz2,nbsphere,nbsphere3,XI
     $           ,XR,wrk,FF,FF0,FFloc,polarisa,methodeit,tol,tol1,nloop
     $           ,ncompte ,planf,planb,nstop,infostr)

         else
            write(*,*) 'ddddd'
            call inverserigopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $           ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $           ,ntotalm,ntotal,ldabi,nlar,nmax,nxm,nym,nzm,nx,ny,nz
     $           ,nx2,ny2,nxy2,nz2,nbsphere,nbsphere3,XI,XR,wrk,FF,FF0
     $           ,FFloc,polarisa,methodeit,tol,tol1,nloop,ncompte
     $           ,planf,planb,nstop,infostr)

         endif
         if (nstop.eq.1) return
         
         if (nlecture.eq.1.and.nlecture1.eq.0) then
            file1='.lf'
            long = len( trim(filereread  ) )
            long1 = len( trim( file1 ) )
            filereread1=filereread(1:long)//file1(1:long1)
            open(1000,file=filereread1,status='new',form='unformatted')
            do i=1,nbsphere3
               write(1000) FF(i),FFloc(i)
            enddo
            close(1000)
         endif



         write(99,*) 'methode iterative used',methodeit
         write(99,*) 'Tolerance asked for the iterative method',tolinit
         write(99,*) 'Tolerance obtained for the iterative method',tol1
         write(99,*) 'Number of product Ax for the iterative method'
     $        ,ncompte,nloop
         
         if (tol1.ge.tolinit) then
            nstop=1
            infostr='Converge do not reach!'
            write(*,*) 'Converge do not reach!',tol1,'Tol Asked',tolinit
            write(99,*) tol1,tolinit
            return
         endif
         

         if (nstop == -1) then
            infostr = 'Calculation cancelled after iterative method!'
            return
         endif
         
      elseif (nrig.eq.1) then 
c     **************************************        
c     Renormalized Born approximation field
c     **************************************
         write(*,*) 'Renormalized Born approximation'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere3
            FFloc(i)=FF0(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii,jj)    
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            do ii=1,3
               FF(k+ii)=0.d0
               do jj=1,3
                  FF(k+ii)=FF(k+ii)+polarisa(i,ii,jj)*FFloc(k+jj)
               enddo
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      elseif (nrig.eq.2) then 
c     **************************************        
c     Born approximation field
c     **************************************
         write(*,*) 'Born approximation'
         nsens=-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,kk,Em,Eloc,ii,jj,epsani)   
!$OMP DO SCHEDULE(STATIC) 
         do k=1,nbsphere
            kk=3*(k-1)
            Em(1)=FF0(kk+1)
            Em(2)=FF0(kk+2)
            Em(3)=FF0(kk+3)
            do ii=1,3
               do jj=1,3
                  epsani(ii,jj)=epsilon(k,ii,jj)
               enddo
            enddo 
            call local_macro(Eloc,Em,epsani,aretecube,k0,nsens)
            FFloc(kk+1)=Eloc(1)
            FFloc(kk+2)=Eloc(2)
            FFloc(kk+3)=Eloc(3)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL          

c     dipole
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii,jj)    
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            do ii=1,3
               FF(k+ii)=0.d0
               do jj=1,3
                  FF(k+ii)=FF(k+ii)+polarisa(i,ii,jj)*FFloc(k+jj)
               enddo
            enddo
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      elseif (nrig.eq.3) then
c     renormalized Born ordre 1
         write(*,*) 'Born series order 1'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            xr(i)=FF0(i)           
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)  
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=-polarisa(i,1,1)*xr(k+1)-polarisa(i,1,2)*xr(k+2)
     $           -polarisa(i,1,3)*xr(k+3)
            xi(k+2)=-polarisa(i,2,1)*xr(k+1)-polarisa(i,2,2)*xr(k+2)
     $           -polarisa(i,2,3)*xr(k+3)
            xi(k+3)=-polarisa(i,3,1)*xr(k+1)-polarisa(i,3,2)*xr(k+2)
     $           -polarisa(i,3,3)*xr(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
         
         call produitfftmatvect3(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,Tabdip,ntotalm,ntotal,nmax ,ndipole,nxm ,nym,nzm,nx,ny
     $        ,nz,nx2,ny2,nxy2,nz2,XI,XR,planf,planb)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k) 
!$OMP DO SCHEDULE(STATIC)       
         do i=1,nbsphere
            k=3*(i-1)
            FFloc(k+1)=xr(k+1)
            FFloc(k+2)=xr(k+2)
            FFloc(k+3)=xr(k+3)
            FF(k+1)=polarisa(i,1,1)*xr(k+1)+polarisa(i,1,2)*xr(k+2)
     $           +polarisa(i,1,3)*xr(k+3)
            FF(k+2)=polarisa(i,2,1)*xr(k+1)+polarisa(i,2,2)*xr(k+2)
     $           +polarisa(i,2,3)*xr(k+3)
            FF(k+3)=polarisa(i,3,1)*xr(k+1)+polarisa(i,3,2)*xr(k+2)
     $           +polarisa(i,3,3)*xr(k+3)
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      

         
      elseif (nrig.eq.4) then
c     Rytov renormalize
         write(*,*) 'Rytov approximation on local field'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere3
            xr(i)=0.d0       
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            xi(k+1)=-polarisa(i,1,1)*FF0(k+1)-polarisa(i,1,2)*FF0(k+2)
     $           -polarisa(i,1,3)*FF0(k+3)
            xi(k+2)=-polarisa(i,2,1)*FF0(k+1)-polarisa(i,2,2)*FF0(k+2)
     $           -polarisa(i,2,3)*FF0(k+3)
            xi(k+3)=-polarisa(i,3,1)*FF0(k+1)-polarisa(i,3,2)*FF0(k+2)
     $           -polarisa(i,3,3)*FF0(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
         
         call produitfftmatvect3(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,Tabdip,ntotalm,ntotal,nmax ,ndipole,nxm ,nym,nzm,nx,ny
     $        ,nz,nx2,ny2,nxy2,nz2,XI,XR,planf,planb)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,tmp)
!$OMP DO SCHEDULE(STATIC)      
         do i=1,nbsphere
            k=3*(i-1)
            tmp=dsqrt(dreal(FF0(k+1)*dconjg(FF0(k+1))+FF0(k+2)
     $           *dconjg(FF0(k+2))+FF0(k+3)*dconjg(FF0(k+3))))*1.d-6
            if (cdabs(FF0(k+1)).le.tmp) then
               FFloc(k+1)=0.d0
            else
               FFloc(k+1)=FF0(k+1)*cdexp(xr(k+1)/FF0(k+1))
            endif
            if (cdabs(FF0(k+2)).le.tmp) then
               FFloc(k+2)=0.d0
            else
               FFloc(k+2)=FF0(k+2)*cdexp(xr(k+2)/FF0(k+2))
            endif
            if (cdabs(FF0(k+3)).le.tmp) then
               FFloc(k+3)=0.d0
            else
               FFloc(k+3)=FF0(k+3)*cdexp(xr(k+3)/FF0(k+3))
            endif
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)       
         do i=1,nbsphere
            k=3*(i-1)  
            FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)*FFloc(k
     $           +2)+polarisa(i,1,3)*FFloc(k+3)
            FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)*FFloc(k
     $           +2)+polarisa(i,2,3)*FFloc(k+3)
            FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)*FFloc(k
     $           +2)+polarisa(i,3,3)*FFloc(k+3)
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      

      elseif (nrig.eq.5) then
c     Rytov
         write(*,*) 'Rytov approximation on macroscopic field'
c     a3=aretecube*aretecube*aretecube
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)   
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere
            k=3*(i-1)
            xr(k+1)=-FF0(k+1)*(epsilon(i,1,1)-eps0)/3.d0
            xr(k+2)=-FF0(k+2)*(epsilon(i,2,2)-eps0)/3.d0
            xr(k+3)=-FF0(k+3)*(epsilon(i,3,3)-eps0)/3.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         nsens=-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,Em,Eloc,ii,jj,epsani)
!$OMP DO SCHEDULE(STATIC)
         do i=1,nbsphere
            k=3*(i-1)
            Em(1)=FF0(k+1)
            Em(2)=FF0(k+2)
            Em(3)=FF0(k+3)
            do ii=1,3
               do jj=1,3
                  epsani(ii,jj)=epsilon(i,ii,jj)
               enddo
            enddo 
            call local_macro(Eloc,Em,epsani,aretecube,k0,nsens)
            FF(k+1)=Eloc(1)
            FF(k+2)=Eloc(2)
            FF(k+3)=Eloc(3)
            xi(k+1)=-polarisa(i,1,1)*FF(k+1)-polarisa(i,1,2)*FF(k+2)
     $           -polarisa(i,1,3)*FF(k+3)
            xi(k+2)=-polarisa(i,2,1)*FF(k+1)-polarisa(i,2,2)*FF(k+2)
     $           -polarisa(i,2,3)*FF(k+3)
            xi(k+3)=-polarisa(i,3,1)*FF(k+1)-polarisa(i,3,2)*FF(k+2)
     $           -polarisa(i,3,3)*FF(k+3)
         enddo        
!$OMP ENDDO 
!$OMP END PARALLEL
         
         call produitfftmatvect3(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,Tabdip,ntotalm,ntotal,nmax ,ndipole,nxm ,nym,nzm,nx,ny
     $        ,nz,nx2,ny2,nxy2,nz2,XI,XR,planf,planb)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,tmp)
!$OMP DO SCHEDULE(STATIC)      
         do i=1,nbsphere
            k=3*(i-1)
            tmp=dsqrt(dreal(FF(k+1)*dconjg(FF(k+1))+FF(k+2)
     $           *dconjg(FF(k+2))+FF(k+3)*dconjg(FF(k+3))))*1.d-6
            if (cdabs(FF(k+1)).le.tmp) then
               FFloc(k+1)=0.d0
            else
               FFloc(k+1)=FF(k+1)*cdexp(xr(k+1)/FF(k+1))
            endif
            if (cdabs(FF(k+2)).le.tmp) then
               FFloc(k+2)=0.d0
            else
               FFloc(k+2)=FF(k+2)*cdexp(xr(k+2)/FF(k+2))
            endif
            if (cdabs(FF(k+3)).le.tmp) then
               FFloc(k+3)=0.d0
            else
               FFloc(k+3)=FF(k+3)*cdexp(xr(k+3)/FF(k+3))
            endif
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)      
         do i=1,nbsphere
            k=3*(i-1)  
            FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)*FFloc(k
     $           +2)+polarisa(i,1,3)*FFloc(k+3)
            FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)*FFloc(k
     $           +2)+polarisa(i,2,3)*FFloc(k+3)
            FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)*FFloc(k
     $           +2)+polarisa(i,3,3)*FFloc(k+3)
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL      
         
      elseif (nrig.eq.6) then
c     Beam propagation method
         write(*,*) 'Beam propagation method'

         if (beam(1:7).eq.'antenna'.or.beam(1:15).eq.'wavelinearmulti'
     $        .or. beam(1:15).eq.'gparawavelinear'.or. beam(1:17).eq
     $        .'gparawavecircular' )then
            nstop=1.d0
            infostr='BPM not possible with this beam'
            return
         endif
         
         

         call beampropagationmacro(xs,ys,zs,aretecube,k0,w0,E0,ss,pp
     $        ,theta,phi,xgaus,ygaus,zgaus,beam,epsilon,ndipole ,nx ,ny
     $        ,nz,nxm,nym ,nzm,nmax,nfft2d,Eimagex,Eimagey,Eimagez,FF0
     $        ,FFloc,FF,plan2f,plan2b ,nstop,infostr)

         if (nstop.eq.1) return

         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)       
         do i=1,nbsphere
            k=3*(i-1)  
            FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)*FFloc(k
     $           +2)+polarisa(i,1,3)*FFloc(k+3)
            FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)*FFloc(k
     $           +2)+polarisa(i,2,3)*FFloc(k+3)
            FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)*FFloc(k
     $           +2)+polarisa(i,3,3)*FFloc(k+3)
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL    
      elseif (nrig.eq.7) then
c     Beam propagation method
         write(*,*) 'renormalized Beam propagation method'

         if (beam(1:7).eq.'antenna'.or.beam(1:15).eq.'wavelinearmulti')
     $        then
            nstop=1.d0
            infostr='BPM not possible with this beam'
            return
         endif
         
         

         call beampropagation(xs,ys,zs,aretecube,k0,w0,E0,ss,pp,theta
     $        ,phi,xgaus,ygaus,zgaus,beam,epsilon,ndipole ,nx ,ny,nz,nxm
     $        ,nym ,nzm,nmax,nfft2d,Eimagex,Eimagey,Eimagez,FF0 ,FFloc
     $        ,FF,plan2f,plan2b,nstop,infostr)

         if (nstop.eq.1) return

         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)       
         do i=1,nbsphere
            k=3*(i-1)  
            FF(k+1)=polarisa(i,1,1)*FFloc(k+1)+polarisa(i,1,2)*FFloc(k
     $           +2)+polarisa(i,1,3)*FFloc(k+3)
            FF(k+2)=polarisa(i,2,1)*FFloc(k+1)+polarisa(i,2,2)*FFloc(k
     $           +2)+polarisa(i,2,3)*FFloc(k+3)
            FF(k+3)=polarisa(i,3,1)*FFloc(k+1)+polarisa(i,3,2)*FFloc(k
     $           +2)+polarisa(i,3,3)*FFloc(k+3)
         enddo            
!$OMP ENDDO 
!$OMP END PARALLEL    
      endif





      
      call cpu_time(t2)
      call date_and_time(date,time,zone,values2)
      
      write(*,*) 'CPU time to solve Ax=b in second  : ',t2-t1
      write(*,*) 'Real time to solve Ax=b in second : ',(values2(5)
     $     -values(5))*3600.d0+(values2(6)-values(6))*60.d0+(values2(7)
     $     -values(7))+(values2(8)-values(8))/1000.d0

      write(*,*) 'Method iterative used             : ',methodeit
      write(*,*) 'Tolerance obtained                : ',tol1
      write(*,*) 'Tolerance asked                   : ',tolinit
      write(*,*) 'Number of product Ax needs        : ',ncompte
      write(*,*) '******** End to solve the linear system ********'
      write(*,*) ' '
c     compute the near field with FFT
 1000 if (nprochefft.ge.1) then
c     compute the FFT in all the box
         write(*,*) 'Compute Near field with FFT'

         call factor235(nxm,IP,test)        
         if (test.eq.1) then
            infostr='nxm problem!'
            nstop=1
            return
         endif
         call factor235(nym,IP,test)        
         if (test.eq.1) then
            infostr='nym problem!'
            nstop=1
            return
         endif
         call factor235(nzm,IP,test)        
         if (test.eq.1) then
            infostr='nzm problem!'
            nstop=1
            return
         endif

   
         nxm2=2*nxm
         nym2=2*nym
         nzm2=2*nzm
         nxym2=nxm2*nym2
         call dfftw_plan_dft_3d(planbn, nx2,ny2,nz2, vectx,vectx
     $        ,FFTW_BACKWARD,FFTW_ESTIMATE)
         call dfftw_plan_dft_3d(planfn, nx2,ny2,nz2, vectx,vectx
     $        ,FFTW_FORWARD,FFTW_ESTIMATE)

         if (nquad.eq.0) then
            call greencalculfft(nxm,nym,nzm,nxm2,nym2,nzm2,nxym2,nxm,nym
     $           ,nzm,aretecube,k0,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $           ,FFTTENSORyy ,FFTTENSORyz,FFTTENSORzz,planbn)
         else
            call greencalculquadfft(nxm,nym,nzm,nxm2,nym2,nzm2,nxym2,nxm
     $           ,nym,nzm,nquad,tolinit,aretecube,k0,FFTTENSORxx
     $           ,FFTTENSORxy,FFTTENSORxz ,FFTTENSORyy,FFTTENSORyz
     $           ,FFTTENSORzz,planbn) 
         endif

            
         
c     compte the size and edge of the box
         xmax=-1.d300
         xmin=1.d300
         ymax=-1.d300
         ymin=1.d300
         zmax=-1.d300
         zmin=1.d300      
!$OMP PARALLEL 
!$OMP DO SCHEDULE(STATIC) REDUCTION(max:xmax,ymax,zmax)
!$OMP& REDUCTION(min:xmin,ymin,zmin)  
         do i=1,nbsphere
            xmax=max(xmax,xs(i))
            xmin=min(xmin,xs(i))
            ymax=max(ymax,ys(i))
            ymin=min(ymin,ys(i))
            zmax=max(zmax,zs(i))
            zmin=min(zmin,zs(i))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL   
         
         write(*,*) 'Size of the box for the near field'
         write(*,*) 'xmin and xmax =',xmin,xmax
         write(*,*) 'ymin and ymax =',ymin,ymax
         write(*,*) 'zmin and zmax =',zmin,zmax
         write(*,*) 'meshsize = ',aretecube
c     compute new position of computation, dipole and incident field at
c     the new grid
         write(*,*) 'Discretization of the box for the near field'
         write(*,*) 'nx =',nx
         write(*,*) 'ny =',ny
         write(*,*) 'nz =',nz
         kkm=(nxm-nx)/2
         jjm=(nym-ny)/2
         iim=(nzm-nz)/2
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,3*nmax
            xi(i)=0.d0
            xr(i)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
         l=1
         subunit=0
c     initialize and reuse polarisa : pola= background
        
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
         do i=1,nmax            
            polarisa(i,1,1)=1.d0
            polarisa(i,2,2)=1.d0
            polarisa(i,3,3)=1.d0
            polarisa(i,1,2)=0.d0
            polarisa(i,1,3)=0.d0
            polarisa(i,2,1)=0.d0
            polarisa(i,2,3)=0.d0
            polarisa(i,3,1)=0.d0
            polarisa(i,3,2)=0.d0
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         cntwf = 0
         
         if (beam(1:8).eq.'gfftwave') then
            do i=1,nzm
               do j=1,nym
                  do k=1,nxm
                     cntwf = cntwf + 1
                     zswf(cntwf) = zmin+dble(i-iim-1)*aretecube
                     yswf(cntwf) = ymin+dble(j-jjm-1)*aretecube
                     xswf(cntwf) = xmin+dble(k-kkm-1)*aretecube
                  enddo
               enddo
            enddo
            call  gaussianchampfft(xswf,yswf,zswf,aretecube,k0,w0,E0,ss
     $           ,pp,theta,phi,xgaus,ygaus,zgaus,beam,ndipole,nxm,nym
     $           ,nzm,nxm,nym ,nzm,nmax,nfft2d,Eimagex,Eimagey,Eimagez
     $           ,xr,tolinit,plan2f,plan2b,nstop,infostr)

            subunit=0
            cntwf=0
            do i=1,nzm
               do j=1,nym
                  do k=1,nxm
                     cntwf = cntwf + 1
                     x=xmin+dble(k-kkm-1)*aretecube
                     y=ymin+dble(j-jjm-1)*aretecube
                     z=zmin+dble(i-iim-1)*aretecube
                     zswf(cntwf) = z
                     yswf(cntwf) = y
                     xswf(cntwf) = x
                     
                     if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(69,*) z
                     if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(68,*) y
                     if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(67,*) x
                     
                     subunit=subunit+1              
                     nsubunit=3*(subunit-1)
                 
                     if (comparaison(x,y,z,xs(l),ys(l),zs(l)
     $                    ,lambda).eq.1.and.l.le.nbsphere)  then
                        ll=3*(l-1)        
c     signe - pour compenser le signe moins dans la routibe AX: E=E0-(-Ap)
                        xi(nsubunit+1)=-FF(ll+1)
                        xi(nsubunit+2)=-FF(ll+2)
                        xi(nsubunit+3)=-FF(ll+3)    
                        do ii=1,3
                           do jj=1,3
                              polarisa(subunit,ii,jj)=epsilon(l,ii,jj)
                           enddo
                        enddo
                        l=l+1
                     endif
                     
                  enddo
               enddo
            enddo

         else
            
            do i=1,nzm
               do j=1,nym
                  do k=1,nxm
                     cntwf = cntwf + 1
                     x=xmin+dble(k-kkm-1)*aretecube
                     y=ymin+dble(j-jjm-1)*aretecube
                     z=zmin+dble(i-iim-1)*aretecube
                     zswf(cntwf) = z
                     yswf(cntwf) = y
                     xswf(cntwf) = x

                     if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(69,*) z
                     if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(68,*) y
                     if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(67,*) x
                     
                     subunit=subunit+1              
                     nsubunit=3*(subunit-1)
c     write(*,*) 'coucou', x,y,z,xs(l),ys(l),zs(l),l
c     $                 ,comparaison(x,y,z,xs(l),ys(l),zs(l),lambda)
                     if (comparaison(x,y,z,xs(l),ys(l),zs(l)
     $                    ,lambda).eq.1.and.l.le.nbsphere)  then
                        ll=3*(l-1)        
c     signe - pour compenser le signe moins dans la routibe AX: E=E0-(-Ap)
                        xi(nsubunit+1)=-FF(ll+1)
                        xi(nsubunit+2)=-FF(ll+2)
                        xi(nsubunit+3)=-FF(ll+3)
                        xr(nsubunit+1)=FF0(ll+1)
                        xr(nsubunit+2)=FF0(ll+2)
                        xr(nsubunit+3)=FF0(ll+3)
                        do ii=1,3
                           do jj=1,3
                              polarisa(subunit,ii,jj)=epsilon(l,ii,jj)
                           enddo
                        enddo
                        l=l+1
                     else
c     comptute the incident field not computed yet
                        if (beam(1:11).eq.'pwavelinear') then
                           call ondeplane(x,y,z,k0,E0,ss,pp,theta,phi
     $                          ,xr(nsubunit+1),xr(nsubunit+2)
     $                          ,xr(nsubunit +3),nstop,infostr)
                        elseif (beam(1:13).eq.'pwavecircular') then
                           call ondecirce(x,y,z,k0,E0,ss,theta,phi
     $                          ,xr(nsubunit+1),xr(nsubunit+2)
     $                          ,xr(nsubunit +3))
                        elseif (beam(1:15).eq.'wavelinearmulti') then
                           call ondeplanemulti(x,y,z,k0,E0m,ssm,ppm
     $                          ,thetam,phim,nbinc,xr(nsubunit+1)
     $                          ,xr(nsubunit+2),xr(nsubunit+3),nstop
     $                          ,infostr)
                        elseif (beam(1:7).eq.'antenna') then
                           call dipoleinc(xdip,ydip,zdip,theta,phi,x ,y
     $                          ,z,aretecube,k0,E0,xr(nsubunit+1)
     $                          ,xr(nsubunit+2),xr(nsubunit +3),nstop
     $                          ,infostr)
                        elseif (beam(1:11).eq.'gwavelinear') then
                           call gaussianchamp(x,y,z,xgaus,ygaus,zgaus
     $                          ,theta,phi,w0,k0,ss,pp,E0,xr(nsubunit+1)
     $                          ,xr(nsubunit+2),xr(nsubunit+3),tol,nloin
     $                          ,nstop,infostr)
                        elseif (beam(1:13).eq.'gwavecircular') then
                           call gaussianchampcirc(x,y,z,xgaus,ygaus
     $                          ,zgaus,theta,phi,w0,k0,ss,E0,xr(nsubunit
     $                          +1),xr(nsubunit+2),xr(nsubunit+3),tol
     $                          ,nloin)
                        elseif (beam(1:15).eq.'gparawavelinear') then
                           call gaussianparalinear(x,y,z,xgaus,ygaus
     $                          ,zgaus,theta,phi,w0,k0,ss,pp,E0
     $                          ,xr(nsubunit+1),xr(nsubunit+2)
     $                          ,xr(nsubunit+3),nstop,infostr)
                        elseif (beam(1:17).eq.'gparawavecircular') then
                           call gaussianparacirc(x,y,z,xgaus,ygaus,zgaus
     $                          ,theta,phi,w0,k0,ss,E0,xr(nsubunit+1)
     $                          ,xr(nsubunit+2),xr(nsubunit+3),nstop
     $                          ,infostr)
                        elseif (beam(1:13).eq.'arbitrary') then
                           call incidentarbitrarypos(x,y,z,aretecube
     $                          ,xr(nsubunit+1),xr(nsubunit+2)
     $                          ,xr(nsubunit +3),nstop,namefileinc
     $                          ,infostr)
                        endif
                     endif
                  enddo
               enddo
            enddo
         endif
c     Close Far field discretization
         close(67)
         close(68)
         close(69)         
         
c     save the incident field for wide field
         if (nmat.eq.0) then
            do i=1,subunit
               incidentfieldx(i) = xr(3*i-2)
               incidentfieldy(i) = xr(3*i-1)
               incidentfieldz(i) = xr(3*i)
               incidentfield(i) = dsqrt(dreal(xr(3*i-2)*dconjg(xr(3*i
     $              -2))+xr(3*i-1)*dconjg(xr(3*i-1))+xr(3*i)*dconjg(xr(3
     $              *i))))
               write(136,*) dreal(xr(3*i-2)),dimag(xr(3*i-2))
               write(137,*) dreal(xr(3*i-1)),dimag(xr(3*i-1))
               write(138,*) dreal(xr(3*i)),dimag(xr(3*i))     
               write(139,*) incidentfield(i)
            enddo
         else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)            
            do i=1,subunit
               incidentfieldx(i) = xr(3*i-2)
               incidentfieldy(i) = xr(3*i-1)
               incidentfieldz(i) = xr(3*i)
               incidentfield(i) = dsqrt(dreal(xr(3*i-2)*dconjg(xr(3*i
     $              -2))+xr(3*i-1)*dconjg(xr(3*i-1))+xr(3*i)*dconjg(xr(3
     $              *i))))          
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL              
         endif
c     Close Intensity of the incident field wide field
         close(136)
         close(137)
         close(138)
         close(139)         
c     calcul du champ la position local par le produit matrice vecteur
c     FFT:

         call produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,ntotalm,ntotalm,nmax,nxm,nym,nzm,nxm,nym,nzm,nxm2,nym2
     $        ,nxym2,nzm2,xi,xr,planfn,planbn)
         write(*,*) 'nlocal',nlocal,subunit,xr(1),xr(2),xr(3)
         if (nlocal.eq.1) then
            if (nmat.eq.0) then
               do i=1,subunit
                  ii=3*(i-1)
                  localfieldx(i) = xr(ii+1)
                  localfieldy(i) = xr(ii+2)
                  localfieldz(i) = xr(ii+3)
                  localfield(i) = dsqrt(dreal(xr(ii+1)*dconjg(xr(ii+1))
     $                 +xr(ii+2)*dconjg(xr(ii+2))+xr(ii+3)*dconjg(xr(ii
     $                 +3))))
                  write(140,*) dreal(xr(ii+1)),dimag(xr(ii+1))
                  write(141,*) dreal(xr(ii+2)),dimag(xr(ii+2))
                  write(142,*) dreal(xr(ii+3)),dimag(xr(ii+3))
                  write(143,*) localfield(i)
               enddo
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)   
!$OMP DO SCHEDULE(STATIC)
               do i=1,subunit
                  ii=3*(i-1)
                  localfieldx(i) = xr(ii+1)
                  localfieldy(i) = xr(ii+2)
                  localfieldz(i) = xr(ii+3)
                  localfield(i) = dsqrt(dreal(xr(ii+1)*dconjg(xr(ii+1))
     $                 +xr(ii+2)*dconjg(xr(ii+2))+xr(ii+3)*dconjg(xr(ii
     $                 +3))))
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL                
            endif
c     Close Intensity of the local field wide field
            close(140)
            close(141)
            close(142)
            close(143)
         endif
         write(*,*) 'nlocal',nmacro,subunit
         if (nmacro.eq.1) then
            nsens=1
            if (nmat.eq.0) then
               do i=1,subunit              
                  ii=3*(i-1)
                  Eloc(1)= xr(ii+1)
                  Eloc(2)= xr(ii+2)
                  Eloc(3)= xr(ii+3)
c     pour l'instant faut au niveau du epsilonpour le champ macro
                  do ii=1,3
                     do jj=1,3
                        epsani(ii,jj)=polarisa(i,ii,jj)
                     enddo
                  enddo 
                  call local_macro(Eloc,Em,epsani,aretecube,k0,nsens)
                  macroscopicfieldx(i) = Em(1)
                  macroscopicfieldy(i) = Em(2)
                  macroscopicfieldz(i) = Em(3)
                  macroscopicfield(i) = dsqrt(dreal(Em(1)*dconjg(Em(1))
     $                 +Em(2)*dconjg(Em(2))+Em(3)*dconjg(Em(3))))
                  write(144,*) dreal(Em(1)),dimag(Em(1))
                  write(145,*) dreal(Em(2)),dimag(Em(2))
                  write(146,*) dreal(Em(3)),dimag(Em(3))             
                  write(147,*) macroscopicfield(i)
               enddo
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii,Eloc,epsani,Em)   
!$OMP DO SCHEDULE(STATIC)               
               do i=1,subunit              
                  ii=3*(i-1)
                  Eloc(1)= xr(ii+1)
                  Eloc(2)= xr(ii+2)
                  Eloc(3)= xr(ii+3)
c     pour l'instant faut au niveau du epsilonpour le champ macro
                  do ii=1,3
                     do jj=1,3
                        epsani(ii,jj)=polarisa(i,ii,jj)
                     enddo
                  enddo 
                  call local_macro(Eloc,Em,epsani,aretecube,k0,nsens)
                  macroscopicfieldx(i) = Em(1)
                  macroscopicfieldy(i) = Em(2)
                  macroscopicfieldz(i) = Em(3)
                  macroscopicfield(i) = dsqrt(dreal(Em(1)*dconjg(Em(1))
     $                 +Em(2)*dconjg(Em(2))+Em(3)*dconjg(Em(3))))
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL 
            endif
c     Intensity of the macroscopic field wide field
            close(144)
            close(145)
            close(146)
            close(147)            
         endif
         write(*,*) 'End  computation of Near field'
         write(*,*) ' '
      else
c     ***************************************************
c     fin de condition sur le champ proche large
c     debut calcul du champ local et macro dans juste l'objet      
         if (nlocal.eq.1) then
            if (nmat.eq.0) then
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     localfieldx(k)=FFloc(ii+1)
                     localfieldy(k)=FFloc(ii+2)
                     localfieldz(k)=FFloc(ii+3)
                     localfield(k) = dsqrt(dreal(FFloc(ii+1)
     $                    *dconjg(FFloc(ii+1))+FFloc(ii+2)
     $                    *dconjg(FFloc(ii +2))+FFloc(ii+3)
     $                    *dconjg(FFloc(ii+3))))
                     write(40,*) dreal(FFloc(ii+1)),dimag(FFloc(ii+1))
                     write(41,*) dreal(FFloc(ii+2)),dimag(FFloc(ii+2))
                     write(42,*) dreal(FFloc(ii+3)),dimag(FFloc(ii+3))
                     write(43,*) localfield(k)
                  else
                     write(40,*) 0.d0,0.d0
                     write(41,*) 0.d0,0.d0
                     write(42,*) 0.d0,0.d0
                     write(43,*) 0.d0
                  endif
               enddo
            else

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii)   
!$OMP DO SCHEDULE(STATIC)
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     localfieldx(k)=FFloc(ii+1)
                     localfieldy(k)=FFloc(ii+2)
                     localfieldz(k)=FFloc(ii+3)
                     localfield(k) = dsqrt(dreal(FFloc(ii+1)
     $                    *dconjg(FFloc(ii+1))+FFloc(ii+2)
     $                    *dconjg(FFloc(ii +2))+FFloc(ii+3)
     $                    *dconjg(FFloc(ii+3))))
                  endif
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL                 
            endif
c     Close Intensity of the local field
            close(40)
            close(41)
            close(42)
            close(43)            
         endif
c     compute and save the macroscopic field
         
         if (nmacro.eq.1) then
            nsens=1
            if (nmat.eq.0) then
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     Eloc(1)= FFloc(ii+1)
                     Eloc(2)= FFloc(ii+2)
                     Eloc(3)= FFloc(ii+3)
                     do ii=1,3
                        do jj=1,3
                           epsani(ii,jj)=epsilon(k,ii,jj)
                        enddo
                     enddo 
                     
                     call local_macro(Eloc,Em,epsani,aretecube,k0,nsens)
                     macroscopicfieldx(k)=Em(1)
                     macroscopicfieldy(k)=Em(2)
                     macroscopicfieldz(k)=Em(3)
                     write(44,*) dreal(Em(1)),dimag(Em(1))
                     write(45,*) dreal(Em(2)),dimag(Em(2))
                     write(46,*) dreal(Em(3)),dimag(Em(3))
                     macroscopicfield(k)= dsqrt(dreal(Em(1)
     $                    *dconjg(Em(1))+Em(2)*dconjg(Em(2))+Em(3)
     $                    *dconjg(Em(3))))
                     write(47,*) macroscopicfield(k)
                  else
                     write(44,*) 0.d0,0.d0
                     write(45,*) 0.d0,0.d0
                     write(46,*) 0.d0,0.d0
                     write(47,*) 0.d0
                  endif
               enddo
            else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,ii,Eloc,epsani,Em)   
!$OMP DO SCHEDULE(STATIC)
               do i=1,ndipole
                  k=tabdip(i)
                  if (k.ne.0) then
                     ii=3*(k-1)
                     Eloc(1)= FFloc(ii+1)
                     Eloc(2)= FFloc(ii+2)
                     Eloc(3)= FFloc(ii+3)
                     do ii=1,3
                        do jj=1,3
                           epsani(ii,jj)=epsilon(k,ii,jj)
                        enddo
                     enddo                      
                     call local_macro(Eloc,Em,epsani,aretecube,k0,nsens)
                     macroscopicfieldx(k)=Em(1)
                     macroscopicfieldy(k)=Em(2)
                     macroscopicfieldz(k)=Em(3)
                     macroscopicfield(k)= dsqrt(dreal(Em(1)
     $                    *dconjg(Em(1))+Em(2)*dconjg(Em(2))+Em(3)
     $                    *dconjg(Em(3))))
                  endif
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL                 
            endif
c     Close Intensity of the macroscopic field
            close(44)
            close(45)
            close(46)
            close(47)            
         endif
      endif

      if (nsection.eq.1.or.nsectionsca.eq.1) then
         if (object(1:6).eq.'sphere') then
            CALL CALLBHMIE(1.d0,eps,rayon,lambda,MIECEXT,MIECABS,MIECSCA
     $           ,GSCA)          
         endif
      endif

      
      if (nsection.eq.1) then
c     Compute the extinction cross section and absorbing cross section
         write(*,*) '************* Compute cross section ********'
         Cext=0.d0   
         Cabs=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,kk)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Cext,Cabs)               
         do i=1,nbsphere
            kk=3*(i-1)
            Cext=Cext+dimag(dconjg(FF0(kk+1))*FF(kk+1)+dconjg(FF0(kk
     $           +2))*FF(kk+2)+dconjg(FF0(kk+3))*FF(kk+3))
            Cabs=Cabs+dimag(FF(kk+1)*dconjg(FFloc(kk+1))+FF(kk+2)
     $           *dconjg(FFloc(kk+2))+FF(kk+3)*dconjg(FFloc(kk+3)))
     $           -2.d0/3.d0*k03*dreal(FF(kk+1)*dconjg(FF(kk+1))+FF(kk+2)
     $           *dconjg(FF(kk+2))+FF(kk+3)*dconjg(FF(kk+3)))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL         
         Cext=4.d0*pi*k0*Cext/I0
         Cabs=4.d0*pi*k0*Cabs/I0

c     compute the scattering cross section
         Csca=Cext-Cabs
         write(99,*) 'extinction cross section',Cext
         write(99,*) 'scattering cross section',Csca
         if (Cabs.le.Cext*1.d-12) then
            Cabs=0.d0
            write(99,*) 'absorbing cross section ',Cabs
         else
            write(99,*) 'absorbing cross section ',Cabs
         endif
         if (nstop == -1) then
            infostr = 'Calculation cancelled after cross section!'
            return
         endif
         write(*,*) 'Cext = ',Cext         
         write(*,*) 'Csca = ',Csca
         if (Cabs.le.Cext*1.d-12) then
            Cabs=0.d0
            write(*,*) 'Cabs = ',Cabs
         else
            write(*,*) 'Cabs = ',Cabs
         endif
         
         if (object(1:6).eq.'sphere') then
            write(*,*) 'Comparison with Mie s series'
            write(*,*) 'Cext error in %',(Cext-MIECEXT)/MIECEXT*100.d0
            write(*,*) 'Csca error in %',(Csca-MIECSCA)/MIECSCA*100.d0
            if (Cabs.ne.0.d0) then
              write(*,*) 'Cabs error in %',(Cabs-MIECABS)/MIECABS*100.d0
            endif

            write(*,*) '********* End cross section ************'
            write(*,*) ' '
         endif
      endif

c     compute the scattering cross section by integration in far field
c     and compute the Poynting vector along the normal in CGS n.P=c/(8
c     pi) Re(ExH^*)=c/(8 pi)|E|^2
      
      if (nsectionsca.eq.1) then
         write(*,*) '***** Compute Csca g and Poynting *******'
         if (nquickdiffracte.eq.0) then
c     compute the diffracted field in summing the dipole
            write(*,*) 'Slow method with discretization:',ntheta,nphi
            if (ntheta.le.10.or.nphi.le.20) then
               infostr='ntheta and nphi too small!'
               nstop=1
               return
            endif

            deltatheta=pi/dble(ntheta)
            deltaphi=2.d0*pi/dble(nphi)
            Cscai=0.d0
            Poyntinginc=0.d0
            gasym=0.d0
            cnt = 0



            do itheta=0,ntheta
               thetas=deltatheta*dble(itheta)
               do iphi=0,nphi-1
                  cnt = cnt + 1
                  phis=deltaphi*dble(iphi)
                  normal(1)=dsin(thetas)*dcos(phis)
                  normal(2)=dsin(thetas)*dsin(phis)
                  normal(3)=dcos(thetas)
                  Emx=0.d0
                  Emy=0.d0
                  Emz=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,kk,ctmp,ctmp1)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Emx,Emy,Emz)      
                  do i=1,nbsphere
                     kk=3*(i-1)
                     ctmp=cdexp(-icomp*k0*(normal(1)*xs(i)+normal(2)
     $                    *ys(i)+normal(3)*zs(i)))
                     ctmp1=normal(1)*FF(kk+1)+normal(2)*FF(kk+2)
     $                    +normal(3)*FF(kk+3)
                     Emx=Emx+ctmp*(FF(kk+1)-ctmp1*normal(1))
                     Emy=Emy+ctmp*(FF(kk+2)-ctmp1*normal(2))
                     Emz=Emz+ctmp*(FF(kk+3)-ctmp1*normal(3))
                  enddo
!$OMP ENDDO 
!$OMP END PARALLEL                  
                  Emod=dreal(Emx*dconjg(Emx)+Emy*dconjg(Emy)
     $                 +Emz*dconjg(Emz))
                  if (nmat.eq.0) then
                     write(50,*) Emod*k0*k0*k0*k0*c/8/pi*quatpieps0
                     write(51,*) thetas*180.d0/pi
                     write(52,*) phis*180.d0/pi
                  endif
                  poyntingfield(cnt) = Emod*k0*k0*k0*k0*c/8/pi
     $                 *quatpieps0
                  thetafield(cnt) = thetas*180.d0/pi
                  phifield(cnt) = phis*180.d0/pi
                  Cscai=Cscai+deltaphi*deltatheta*dsin(thetas)*Emod
                  gasym=gasym+deltaphi*deltatheta*dsin(thetas)*Emod
     $                 *(normal(1)*dsin(theta*pi/180.d0)*dcos(phi*pi
     $                 /180.d0)+normal(2)*dsin(theta*pi/180.d0)*
     $                 dsin(phi *pi/180.d0)+normal(3)*dcos(theta*pi
     $                 /180.d0))
               enddo
            enddo

c     save the Poynting vecteur
            close(50)
            close(51)
            close(52)
            
            gasym=gasym/cscai
            Cscai=Cscai*k03*k0/I0
            write(99,*) 'scattering cross section with integration'
     $           ,Cscai
            write(99,*) 'scattering asymetric parameter',gasym

            if (beam(1:11).ne.'pwavelinear'.and.beam(1:13).ne
     $           .'pwavecircular') then
               gasym=0.d0
               Cscai=0.d0
            else
               write(*,*) 'Csca =',Cscai
               write(*,*) 'g    =',gasym
               if (object(1:6).eq.'sphere') then             
                  write(*,*) 'Comparison with Mie s series'
                  write(*,*) 'Csca error in %',(Cscai-MIECSCA)/MIECSCA
     $                 *100.d0
                  write(*,*) 'g error in %',(gasym-GSCA)/GSCA*100.d0
               endif
            endif
         else
c     compute the diffracted field with FFT
            write(*,*) 'Diffracted field: Quick method with FFT'
            call diffractefft2d(nx,ny,nz,nxm,nym,nzm,nfft2d,k0,xs,ys,zs
     $           ,aretecube,Efourierx,Efouriery,Efourierz,FF,imaxk0
     $           ,deltakx,deltaky,Ediffkzpos,Ediffkzneg,plan2f,plan2b
     $           ,nstop,infostr)
            write(*,*) 'End Quick method with FFT'
            if (nstop.eq.1) return
c     put the field in file with the right angles
            write(*,*) 'Number of point in the numerical aperture'
     $           ,2*imaxk0+1
            if (nmat.eq.0) then
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky
                     if (j.eq.0) write(55,*) kx
                     if (i.eq.0) write(56,*) ky
                     if (k0*k0-kx*kx-ky*ky.gt.0.d0) then 
                        ii=imaxk0+i+1
                        jj=imaxk0+j+1
                        Emod=cdabs(Ediffkzpos(ii,jj,1))**2
     $                       +cdabs(Ediffkzpos(ii,jj,2))**2
     $                       +cdabs(Ediffkzpos(ii,jj,3))**2
                        write(53,*) Emod*k0*k0*k0*k0*c/8/pi*quatpieps0
                        Emod=cdabs(Ediffkzneg(ii,jj,1))**2
     $                       +cdabs(Ediffkzneg(ii,jj,2))**2
     $                       +cdabs(Ediffkzneg(ii,jj,3))**2
                        write(54,*) Emod*k0*k0*k0*k0*c/8/pi*quatpieps0
                     else
                        Emod=0.d0
                        write(53,*) Emod
                        write(54,*) Emod
                     endif
                  enddo
               enddo
            endif
c     save  the Poynting vecteur in kx ky
            close(53)
            close(54)
            close(55)
            close(56)            
            
            call computegcfft2d(imaxk0,deltakx,deltaky,k0,I0,theta,phi
     $           ,nfft2d,Ediffkzpos,Ediffkzneg,gasym,Cscai)

            if (beam(1:11).ne.'pwavelinear'.and.beam(1:13).ne
     $           .'pwavecircular') then
               gasym=0.d0
               Cscai=0.d0
            else
               if (object(1:6).eq.'sphere') then             
                  write(*,*) 'Comparison with Mie s series'
                  write(*,*) 'Csca error in %',(Cscai-MIECSCA)/MIECSCA
     $                 *100.d0
                  write(*,*) 'g error in %',(gasym-GSCA)/GSCA*100.d0
               endif
            endif
            
            deltatheta=pi/dble(ntheta)
            deltaphi=2.d0*pi/dble(nphi)
            cnt = 0

            do itheta=0,ntheta
               thetas=deltatheta*dble(itheta)
               do iphi=0,nphi-1
                  cnt = cnt + 1
                  phis=deltaphi*dble(iphi)
                  if (nmat.eq.0) then
                     write(51,*) thetas*180.d0/pi
                     write(52,*) phis*180.d0/pi
                  endif
                  
                  thetafield(cnt) = thetas*180.d0/pi
                  phifield(cnt) = phis*180.d0/pi

                  kx=k0*dsin(thetas)*dcos(phis)
                  ky=k0*dsin(thetas)*dsin(phis)
                  i=int(kx/deltakx)
                  j=int(ky/deltakx)
                  
                  ii=imaxk0+i+1
                  jj=imaxk0+j+1

                  if (kx.ge.0.d0) i2=i+1
                  if (kx.lt.0.d0) i2=i-1
                  if (ky.ge.0.d0) j2=j+1
                  if (ky.lt.0.d0) j2=j-1

                  if ((dble(i2)**2+dble(j2)**2)*deltakx*deltakx.le.k0
     $                 *k0) then
                     ii2=imaxk0+i2+1
                     jj2=imaxk0+j2+1
                     if (thetas.le.pi/2.d0) then
                        Emod=(cdabs(Ediffkzpos(ii,jj,1))**2
     $                       +cdabs(Ediffkzpos(ii,jj,2))**2
     $                       +cdabs(Ediffkzpos(ii,jj,3))**2)
     $                       *dabs((dble(i2) *deltakx-kx)*(dble(j2)
     $                       *deltakx-ky))+ (cdabs(Ediffkzpos(ii2,jj,1))
     $                       **2 +cdabs(Ediffkzpos(ii2,jj,2))**2
     $                       +cdabs(Ediffkzpos(ii2,jj,3))**2)*dabs((
     $                       -dble(i) *deltakx+kx)*(dble(j2)*deltakx
     $                       -ky))+ (cdabs(Ediffkzpos(ii,jj2,1))**2
     $                       +cdabs(Ediffkzpos(ii,jj2,2))**2
     $                       +cdabs(Ediffkzpos(ii,jj2,3))**2)
     $                       *dabs((dble(i2) *deltakx-kx)*(-dble(j)
     $                       *deltakx+ky))+ (cdabs(Ediffkzpos(ii2,jj2
     $                       ,1))**2 +cdabs(Ediffkzpos(ii2,jj2,2))**2
     $                       +cdabs(Ediffkzpos(ii2,jj2,3))**2)*dabs((
     $                       -dble(i) *deltakx+kx)*(-dble(j)*deltakx
     $                       +ky))
                        
                     else

                        Emod=(cdabs(Ediffkzneg(ii,jj,1))**2
     $                       +cdabs(Ediffkzneg(ii,jj,2))**2
     $                       +cdabs(Ediffkzneg(ii,jj,3))**2)
     $                       *dabs((dble(i2) *deltakx-kx)*(dble(j2)
     $                       *deltakx-ky))+ (cdabs(Ediffkzneg(ii2,jj,1))
     $                       **2 +cdabs(Ediffkzneg(ii2,jj,2))**2
     $                       +cdabs(Ediffkzneg(ii2,jj,3))**2)*dabs((
     $                       -dble(i) *deltakx+kx)*(dble(j2)*deltakx
     $                       -ky))+ (cdabs(Ediffkzneg(ii,jj2,1))**2
     $                       +cdabs(Ediffkzneg(ii,jj2,2))**2
     $                       +cdabs(Ediffkzneg(ii,jj2,3))**2)
     $                       *dabs((dble(i2) *deltakx-kx)*(-dble(j)
     $                       *deltakx+ky))+ (cdabs(Ediffkzneg(ii2,jj2
     $                       ,1))**2 +cdabs(Ediffkzneg(ii2,jj2,2))**2
     $                       +cdabs(Ediffkzneg(ii2,jj2,3))**2)*dabs((
     $                       -dble(i) *deltakx+kx)*(-dble(j)*deltakx
     $                       +ky))
                        
                     endif

                     

                     Emod=Emod/deltakx/deltakx                     
                     poyntingfield(cnt) = Emod*k0*k0*c/8/pi *quatpieps0
                     if (nmat.eq.0) write(50,*) Emod*k0*k0*c/8/pi
     $                    *quatpieps0
                     
                  else
                     if (thetas.le.pi/2.d0) then
                        Emod=cdabs(Ediffkzpos(ii,jj,1))**2
     $                       +cdabs(Ediffkzpos(ii,jj,2))**2
     $                       +cdabs(Ediffkzpos(ii,jj,3))**2
                     else
                        Emod=cdabs(Ediffkzneg(ii,jj,1))**2
     $                       +cdabs(Ediffkzneg(ii,jj,2))**2
     $                       +cdabs(Ediffkzneg(ii,jj,3))**2
                        
                     endif
                     if (nmat.eq.0) write(50,*) Emod*k0*k0*c/8/pi
     $                    *quatpieps0
                     poyntingfield(cnt) = Emod*k0*k0*c/8/pi
     $                    *quatpieps0
                  endif
               enddo
            enddo
            
c     save the Poynting vecteur
            close(50)
            close(51)
            close(52)
            
         endif
         write(*,*) '****** End computation of Csca g and Poynting ****'
         write(*,*) ' '
      endif

      if (nenergie.eq.1) then
         write(*,*) '********* Emissivity ************'
         
         
         call diffractefft2denergie(nx,ny,nz,nxm,nym,nzm,nfft2d,k0,xs,ys
     $        ,zs,E0,ss,pp,theta,phi,thetam ,phim, ppm, ssm,E0m,nbinc
     $        ,xdip,ydip,zdip,xgaus,ygaus,zgaus ,w0,aretecube,tol
     $        ,Efourierx,Efouriery,Efourierz,FF,imaxk0 ,deltakx,deltaky
     $        ,Ediffkzpos,Ediffkzneg,beam,efficacite ,efficaciteref
     $        ,efficacitetrans,nsectionsca,nquickdiffracte,plan2f,plan2b
     $        ,nstop ,infostr)
ccc         write(*,*) 'ttt',Efouriery
         if (nstop.eq.1) return
         write(*,*) 'Absorptivity   : ',efficacite
         write(*,*) 'Reflextivity   : ',efficaciteref
         write(*,*) 'Transmittivity : ',efficacitetrans
         write(*,*) '********* End Emissivity **********'
         write(*,*) ' '
      endif

      
      if (nlentille.eq.1) then
         
         write(*,*) '********* Microscopy *************'
         write(*,*) 'Microscopy with NA=',numaper    
         
         nfft2d2=nfft2d/2
         numaper=numaper*k0
         
         if (nquicklens.eq.0) then
            write(*,*) 'Slow method'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
            do i=1,nfft2d*nfft2d
               Efourierx(i)=0.d0
               Efouriery(i)=0.d0
               Efourierz(i)=0.d0              
               Eimageincx(i)=0.d0
               Eimageincy(i)=0.d0
               Eimageincz(i)=0.d0           
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL
             
            write(*,*) 'Step size delta k : ',2.d0*pi/(dble(nfft2d)
     $           *aretecube)
            k=0
 222        deltakx=2.d0*pi/(dble(nfft2d)*aretecube)/dble(2**k)
            imaxk0=nint(k0/deltakx)+1
            
            if (imaxk0.le.20) then
               k=k+1
               write(*,*) 'change delta k :',k,dble(nfft2d*(2**k))
     $              ,nfft2d
               goto 222
            endif
            write(*,*) 'Final delta k',deltakx
            
            deltaky=deltakx
            deltax=2.d0*pi/dble(nfft2d)/deltakx
c     appelle routine pour calculer le champ incident.
            call  diffractefft2dinc(nfft2d,k0,E0,ss,pp,theta,phi,thetam
     $           ,phim, ppm, ssm,E0m,nbinc,xdip,ydip,zdip,xgaus,ygaus
     $           ,zgaus,w0,aretecube,tol,Efourierx,Efouriery,Efourierz
     $           ,imaxk0,deltakx,deltaky,Ediffkzneg,numaper,beam,plan2f
     $           ,plan2b,nstop ,infostr)
            if (nstop.eq.1) return
            
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
            do i=1,nfft2d*nfft2d
               Efourierx(i)=0.d0
               Efouriery(i)=0.d0
               Efourierz(i)=0.d0
               Efourierincx(i)=0.d0
               Efourierincy(i)=0.d0
               Efourierincz(i)=0.d0
               Eimageincx(i)=0.d0
               Eimageincy(i)=0.d0
               Eimageincz(i)=0.d0
               Eimagex(i)=0.d0
               Eimagey(i)=0.d0
               Eimagez(i)=0.d0
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL            

            do i=-imaxk0,imaxk0
               
               if (i.ge.0) then
                  indicex=i+1
               else
                  indicex=nfft2d+i+1
               endif

               kx=deltakx*dble(i)
               kxy(i+nfft2d2+1)=kx/k0
               
               do j=-imaxk0,imaxk0

                  if (j.ge.0) then
                     indicey=j+1
                  else
                     indicey=nfft2d+j+1
                  endif

                  ky=deltaky*dble(j)
                  kz=1.d0
                  indice=indicex+nfft2d*(indicey-1)
                  if (dsqrt(kx*kx+ky*ky).le.numaper) then
                     kz=dsqrt(k0*k0-kx*kx-ky*ky)
c     write(*,*) 'coucou',i,j,kz
c     compute theta and phi assuming z axis of the lens
                     normal(1)=kx/k0
                     normal(2)=ky/k0
                     normal(3)=dsqrt(1.d0-normal(1)*normal(1)-normal(2)
     $                    *normal(2))
                     Emx=0.d0
                     Emy=0.d0
                     Emz=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,kk,ctmp,ctmp1)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:Emx,Emy,Emz)
                     do ii=1,nbsphere
                        kk=3*(ii-1)
                        ctmp=cdexp(-icomp*k0*(normal(1)*xs(ii)+normal(2)
     $                       *ys(ii)+normal(3)*zs(ii)))
                        ctmp1=normal(1)*FF(kk+1)+normal(2)*FF(kk+2)
     $                       +normal(3)*FF(kk+3)
                        Emx=Emx+ctmp*(FF(kk+1)-ctmp1*normal(1))
                        Emy=Emy+ctmp*(FF(kk+2)-ctmp1*normal(2))
                        Emz=Emz+ctmp*(FF(kk+3)-ctmp1*normal(3))
                     enddo
!$OMP ENDDO 
!$OMP END PARALLEL                     
                     ctmp=-2.d0*pi*icomp*kz                   
                     kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1
c     write(*,*) 'toto',kk,indice,i,j
                     Efourierx(kk)=Emx*k0*k0/ctmp
                     Efouriery(kk)=Emy*k0*k0/ctmp
                     Efourierz(kk)=Emz*k0*k0/ctmp
                     Efourierincx(kk)=Efourierx(kk)+Ediffkzneg(ii,jj,1)
                     Efourierincy(kk)=Efouriery(kk)+Ediffkzneg(ii,jj,2)
                     Efourierincz(kk)=Efourierz(kk)+Ediffkzneg(ii,jj,3)
c     write(*,*) dsqrt(cdabs(Efourierx(kk))**2
c     $                    +cdabs(Efouriery(kk))**2+cdabs(Efourierz(kk))
c     $                    **2)
c     calcul de l'axe de rotation et effectue la rotation 
                     u1=-normal(2)
                     u2=normal(1)
                     tmp=dsqrt(u1*u1+u2*u2)
                     if (tmp.eq.0.d0) then
                        Eimagex(indice)=Efourierx(kk)
                        Eimagey(indice)=Efouriery(kk)
                        Eimagez(indice)=Efourierz(kk)
                        Eimageincx(indice)=Efourierincx(kk)
                        Eimageincy(indice)=Efourierincy(kk)
                        Eimageincz(indice)=Efourierincz(kk)
                     else
                        u1=u1/tmp
                        u2=u2/tmp

c     rotation de theta plus theta'
                        tmp=dasin(dsin(-dacos(normal(3)))/gross)
     $                       -dacos(normal(3))

                        Eimagex(indice)=(u1*u1+(1.d0-u1*u1) *dcos(tmp))
     $                       *Efourierx(kk)+u1*u2*(1.d0 -dcos(tmp))
     $                       *Efouriery(kk)+u2*dsin(tmp) *Efourierz(kk)
                        Eimagey(indice)=u1*u2*(1.d0-dcos(tmp))
     $                       *Efourierx(kk)+(u2*u2+(1.d0-u2*u2)
     $                       *dcos(tmp))*Efouriery(kk)-u1*dsin(tmp)
     $                       *Efourierz(kk)
                        Eimagez(indice)=-u2*dsin(tmp)*Efourierx(kk)+u1
     $                       *dsin(tmp)*Efouriery(kk)+dcos(tmp)
     $                       *Efourierz(kk)
                        Eimageincx(indice)=(u1*u1+(1.d0-u1*u1)
     $                       *dcos(tmp))*Efourierincx(kk)+u1*u2*(1.d0
     $                       -dcos(tmp))*Efourierincy(kk)+u2*dsin(tmp)
     $                       *Efourierincz(kk)
                        Eimageincy(indice)=u1*u2*(1.d0-dcos(tmp))
     $                       *Efourierincx(kk)+(u2*u2+(1.d0-u2*u2)
     $                       *dcos(tmp)) *Efourierincy(kk)-u1 *dsin(tmp)
     $                       *Efourierincz(kk)
                        Eimageincz(indice)=-u2*dsin(tmp)
     $                       *Efourierincx(kk)+u1*dsin(tmp)
     $                       *Efourierincy(kk)+dcos(tmp)
     $                       *Efourierincz(kk)
                     endif
                     

                  endif

                  if (i.eq.-nfft2d2.or.j.eq.-nfft2d2) then
                     Eimageincx(indice)=0.d0
                     Eimageincy(indice)=0.d0
                     Eimageincz(indice)=0.d0
                  endif
               enddo
            enddo
            write(*,*) 'End computation',imaxk0

            do i=-nfft2d2,nfft2d2-1
               xy(i+nfft2d2+1)=-deltax*dble(i)*gross
            enddo

            if (nmat.eq.0) then
               open(400,file='kxfourier.mat')
               open(401,file='ximage.mat')
               do i=-nfft2d2,nfft2d2-1
                  write(401,*) xy(i+nfft2d2+1)
               enddo
               do i=-imaxk0,imaxk0
                  write(400,*) kxy(i+nfft2d2+1)
               enddo
               close(400)
               close(401)

               open(300,file='fourier.mat')         
               open(301,file='fourierx.mat')         
               open(302,file='fouriery.mat')         
               open(303,file='fourierz.mat')
               open(304,file='fourierinc.mat')         
               open(305,file='fourierincx.mat')         
               open(306,file='fourierincy.mat')         
               open(307,file='fourierincz.mat')        
               do j=-imaxk0,imaxk0
                  do i=-imaxk0,imaxk0
                     k=(i+nfft2d2+1)+nfft2d*(j+nfft2d2)
                     write(300,*) dsqrt(dreal(Efourierx(k)
     $                    *dconjg(Efourierx(k))+Efouriery(k)
     $                    *dconjg(Efouriery(k))+Efourierz(k)
     $                    *dconjg(Efourierz(k))))
                     write(301,*)dreal(Efourierx(k)),dimag(Efourierx(k))
                     write(302,*)dreal(Efouriery(k)),dimag(Efouriery(k))
                     write(303,*)dreal(Efourierz(k)),dimag(Efourierz(k))
                     write(304,*) dsqrt(dreal(Efourierincx(k)
     $                    *dconjg(Efourierincx(k))+Efourierincy(k)
     $                    *dconjg(Efourierincy(k))+Efourierincz(k)
     $                    *dconjg(Efourierincz(k))))
                     write(305,*) dreal(Efourierincx(k))
     $                    ,dimag(Efourierincx(k))
                     write(306,*) dreal(Efourierincy(k))
     $                    ,dimag(Efourierincy(k))
                     write(307,*) dreal(Efourierincz(k))
     $                    ,dimag(Efourierincz(k)) 
                  enddo
               enddo

               close(300)
               close(301)
               close(302)
               close(303)
               close(304)
               close(305)
               close(306)
               close(307)
            endif
            
            write(*,*) '***** End microsocopy *********'

         else
            
            write(*,*) 'Quick FFT method for microscopy'

            if (nsectionsca*nquickdiffracte.eq.0.and.nenergie.eq.0) then
               call diffractefft2dlens(nx,ny,nz,nxm,nym,nzm,nfft2d,k0,xs
     $              ,ys,zs,aretecube,Efourierx,Efouriery,Efourierz,FF
     $              ,imaxk0,deltakx,deltaky,Ediffkzpos,numaper,plan2f
     $              ,plan2b,nstop ,infostr)
               if (nstop.eq.1) return

            elseif (nsectionsca*nquickdiffracte.eq.1.and.nenergie.eq.0)
     $              then
               k02=k0*k0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,ii,jj)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)               
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky                     
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1                   
                     if (kx*kx+ky*ky.ge.numaper*numaper) then
                        Ediffkzpos(ii,jj,1)=0.d0
                        Ediffkzpos(ii,jj,2)=0.d0
                        Ediffkzpos(ii,jj,3)=0.d0
                     else
                        kz=dsqrt(k0*k0-kx*kx-ky*ky)
                        ctmp=-2.d0*pi*icomp*kz
                        Ediffkzpos(ii,jj,1)=Ediffkzpos(ii,jj,1)*k02/ctmp
                        Ediffkzpos(ii,jj,2)=Ediffkzpos(ii,jj,2)*k02/ctmp
                        Ediffkzpos(ii,jj,3)=Ediffkzpos(ii,jj,3)*k02/ctmp
                     endif
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
               
            else
c     met a zero ce qui en dehors de l'AN
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,kx,ky,ii,jj)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)               
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky                     
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1                   
                     if (kx*kx+ky*ky.ge.numaper*numaper) then
                        Ediffkzpos(ii,jj,1)=0.d0
                        Ediffkzpos(ii,jj,2)=0.d0
                        Ediffkzpos(ii,jj,3)=0.d0
                     endif                 
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
               
            endif


            
            deltax=2.d0*pi/dble(nfft2d)/deltakx
            write(*,*) 'Size FFT              : ',nfft2d
            write(*,*) 'Step size delta k     : ',deltakx
            write(*,*) 'Step size delta x     : ',deltax
            write(*,*) 'Number of point in NA : ',imaxk0
c     calcul le champ incident

            if (nenergie.eq.0) then
               call  diffractefft2dinc(nfft2d,k0,E0,ss,pp,theta,phi
     $              ,thetam,phim, ppm, ssm,E0m,nbinc,xdip,ydip,zdip
     $              ,xgaus,ygaus,zgaus,w0,aretecube,tol,Efourierx
     $              ,Efouriery,Efourierz,imaxk0,deltakx,deltaky
     $              ,Ediffkzneg ,numaper,beam,plan2f,plan2b,nstop
     $              ,infostr)
c     stop
               if (nstop.eq.1) return

            else
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,kx,ky,ii,jj,indicex,indicey,indice)   
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
               do i=-imaxk0,imaxk0
                  do j=-imaxk0,imaxk0
                     kx=dble(i)*deltakx
                     ky=dble(j)*deltaky                     
                     ii=imaxk0+i+1
                     jj=imaxk0+j+1
                     Ediffkzneg(ii,jj,1)=0.d0
                     Ediffkzneg(ii,jj,2)=0.d0
                     Ediffkzneg(ii,jj,3)=0.d0
                     if (kx*kx+ky*ky.lt.numaper*numaper) then
                        if (i.ge.0) then
                           indicex=i+1
                        else
                           indicex=nfft2d+i+1
                        endif
                        if (j.ge.0) then
                           indicey=j+1
                        else
                           indicey=nfft2d+j+1
                        endif                        
                        indice=indicex+nfft2d*(indicey-1)
                        Ediffkzneg(ii,jj,1)=Efourierx(indice)
                        Ediffkzneg(ii,jj,2)=Efouriery(indice)
                        Ediffkzneg(ii,jj,3)=Efourierz(indice)
                     endif
                     
                  enddo
               enddo
!$OMP ENDDO 
!$OMP END PARALLEL
               
            endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,kx)   
!$OMP DO SCHEDULE(STATIC)
            do i=-nfft2d2,nfft2d2-1
               kx=deltakx*dble(i)
               xy(i+nfft2d2+1)=-deltax*dble(i)*gross
               kxy(i+nfft2d2+1)=kx/k0
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            
            if (nmat.eq.0) then
               open(400,file='kxfourier.mat')
               open(401,file='ximage.mat')
               do i=-nfft2d2,nfft2d2-1
                  write(401,*) xy(i+nfft2d2+1)
               enddo
               do i=-imaxk0,imaxk0
                  write(400,*) kxy(i+nfft2d2+1)
               enddo
               close(400)
               close(401)
            endif
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)   
!$OMP DO SCHEDULE(STATIC)
            do i=1,nfft2d*nfft2d
               Efourierx(i)=0.d0
               Efouriery(i)=0.d0
               Efourierz(i)=0.d0
               Efourierincx(i)=0.d0
               Efourierincy(i)=0.d0
               Efourierincz(i)=0.d0
               Eimageincx(i)=0.d0
               Eimageincy(i)=0.d0
               Eimageincz(i)=0.d0
               Eimagex(i)=0.d0
               Eimagey(i)=0.d0
               Eimagez(i)=0.d0
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,indicex,indicey,indice)
!$OMP& PRIVATE(kx,ky,kz,u1,u2,normal,tmp,ii,jj,kk)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            do i=-imaxk0,imaxk0 
               do j=-imaxk0,imaxk0
                  if (i.ge.0) then
                     indicex=i+1
                  else
                     indicex=nfft2d+i+1
                  endif
                  if (j.ge.0) then
                     indicey=j+1
                  else
                     indicey=nfft2d+j+1
                  endif
                  indice=indicex+nfft2d*(indicey-1)

                  kx=deltakx*dble(i)
                  ky=deltakx*dble(j)
                  kz=k0*k0-kx*kx-ky*ky
                  if (kz.ge.0.d0) then
                     kz=dsqrt(kz)
                  else
                     kz=0.d0
                  endif
                  u1=-ky/k0
                  u2=kx/k0
                  tmp=dsqrt(u1*u1+u2*u2)

                  ii=imaxk0+i+1
                  jj=imaxk0+j+1
                  kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
                  Efourierx(kk)=Ediffkzpos(ii,jj,1)
                  Efouriery(kk)=Ediffkzpos(ii,jj,2)
                  Efourierz(kk)=Ediffkzpos(ii,jj,3)
                  Efourierincx(kk)=Ediffkzpos(ii,jj,1)+Ediffkzneg(ii,jj
     $                 ,1)
                  Efourierincy(kk)=Ediffkzpos(ii,jj,2)+Ediffkzneg(ii,jj
     $                 ,2)
                  Efourierincz(kk)=Ediffkzpos(ii,jj,3)+Ediffkzneg(ii,jj
     $                 ,3)
                  if (tmp.eq.0.d0) then                    
                     Eimagex(indice)=Efourierx(kk)
                     Eimagey(indice)=Efouriery(kk)
                     Eimagez(indice)=Efourierz(kk)
                     Eimageincx(indice)=Efourierincx(kk)
                     Eimageincy(indice)=Efourierincy(kk)
                     Eimageincz(indice)=Efourierincz(kk)
                  else
                     u1=u1/tmp
                     u2=u2/tmp
                     normal(3)=kz/k0
                     tmp=dasin(dsin(-dacos(normal(3)))/gross)
     $                    -dacos(normal(3))

                     Eimageincx(indice)=(u1*u1+(1.d0-u1*u1) *dcos(tmp))
     $                    *Efourierincx(kk)+u1*u2*(1.d0-dcos(tmp))
     $                    *Efourierincy(kk)+u2*dsin(tmp)
     $                    *Efourierincz(kk)
                     Eimageincy(indice)=u1*u2*(1.d0-dcos(tmp))
     $                    *Efourierincx(kk)+(u2*u2+(1.d0-u2*u2)
     $                    *dcos(tmp)) *Efourierincy(kk)-u1 *dsin(tmp)
     $                    *Efourierincz(kk)
                     Eimageincz(indice)=-u2*dsin(tmp)*Efourierincx(kk)
     $                    +u1*dsin(tmp)*Efourierincy(kk)+dcos(tmp)
     $                    *Efourierincz(kk)
                     Eimagex(indice)=(u1*u1+(1.d0-u1*u1)*dcos(tmp))
     $                    *Efourierx(kk)+u1*u2*(1.d0-dcos(tmp))
     $                    *Efouriery(kk)+u2*dsin(tmp)*Efourierz(kk)
                     Eimagey(indice)=u1*u2*(1.d0-dcos(tmp))
     $                    *Efourierx(kk)+(u2*u2+(1.d0-u2*u2) *dcos(tmp))
     $                    *Efouriery(kk)-u1 *dsin(tmp)*Efourierz(kk)
                     Eimagez(indice)=-u2*dsin(tmp)*Efourierx(kk)
     $                    +u1*dsin(tmp)*Efouriery(kk)+dcos(tmp)
     $                    *Efourierz(kk)

                  endif
               enddo
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL

            if (nmat.eq.0) then
               open(300,file='fourier.mat')         
               open(301,file='fourierx.mat')         
               open(302,file='fouriery.mat')         
               open(303,file='fourierz.mat')
               open(304,file='fourierinc.mat')         
               open(305,file='fourierincx.mat')         
               open(306,file='fourierincy.mat')         
               open(307,file='fourierincz.mat')  
               do j=-imaxk0,imaxk0
                  do i=-imaxk0,imaxk0
                     k=(i+nfft2d2+1)+nfft2d*(j+nfft2d2)
                     write(300,*) dsqrt(dreal(Efourierx(k)
     $                    *dconjg(Efourierx(k))+Efouriery(k)
     $                    *dconjg(Efouriery(k))+Efourierz(k)
     $                    *dconjg(Efourierz(k))))
                     write(301,*)dreal(Efourierx(k)),dimag(Efourierx(k))
                     write(302,*)dreal(Efouriery(k)),dimag(Efouriery(k))
                     write(303,*)dreal(Efourierz(k)),dimag(Efourierz(k))
                     write(304,*) dsqrt(dreal(Efourierincx(k)
     $                    *dconjg(Efourierincx(k))+Efourierincy(k)
     $                    *dconjg(Efourierincy(k))+Efourierincz(k)
     $                    *dconjg(Efourierincz(k))))
                     write(305,*) dreal(Efourierincx(k))
     $                    ,dimag(Efourierincx(k))
                     write(306,*) dreal(Efourierincy(k))
     $                    ,dimag(Efourierincy(k))
                     write(307,*) dreal(Efourierincz(k))
     $                    ,dimag(Efourierincz(k)) 
                  enddo
               enddo      

               close(300)
               close(301)
               close(302)
               close(303)
               close(304)
               close(305)
               close(306)
               close(307)
            endif
         endif

         write(*,*) 'Compute the images through the microscope'

         call fouriertoimage(deltakx,deltaky,Eimagex,Eimagey,Eimagez
     $        ,Eimageincx,Eimageincy,Eimageincz,nfft2D,nfft2d2,plan2b
     $        ,plan2f)
         
c     calcul intensite de l'image
         if (nmat.eq.0) then
            open(301,file='image.mat')
            open(302,file='imagex.mat')
            open(303,file='imagey.mat')
            open(304,file='imagez.mat')
            open(305,file='imageinc.mat')
            open(306,file='imageincx.mat')
            open(307,file='imageincy.mat')
            open(308,file='imageincz.mat')

            do i=1,nfft2D*nfft2D
               write(301,*) dsqrt(dreal(Eimagex(i)*dconjg(Eimagex(i))
     $              +Eimagey(i)*dconjg(Eimagey(i))+Eimagez(i)
     $              *dconjg(Eimagez(i))))
               write(302,*) dreal(Eimagex(i)),dimag(Eimagex(i))
               write(303,*) dreal(Eimagey(i)),dimag(Eimagey(i))
               write(304,*) dreal(Eimagez(i)),dimag(Eimagez(i))
               write(305,*) dsqrt(dreal(Eimageincx(i)
     $              *dconjg(Eimageincx(i))+Eimageincy(i)
     $              *dconjg(Eimageincy(i))+Eimageincz(i)
     $              *dconjg(Eimageincz(i))))
               write(306,*) dreal(Eimageincx(i)),dimag(Eimageincx(i))
               write(307,*) dreal(Eimageincy(i)),dimag(Eimageincy(i))
               write(308,*) dreal(Eimageincz(i)),dimag(Eimageincz(i))
            enddo
            close(301)
            close(302)
            close(303)
            close(304)
            close(305)
            close(306)
            close(307)
            close(308)
         endif
         write(*,*) 'End image through the microsocpe'
         write(*,*) '************* End microscopy **************'
         write(*,*) ' '
      endif

      

      if (nforce.eq.1) then
         write(99,*) '****** Compute the optical force *********'
         write(*,*) '******** Compute the optical force ***********'
c     Begin the computation of the force
c     Compute the FFT of the dipole

c     traite faisceau particulier en premier
         if (beam(1:9).eq.'arbitrary') then
            call incidentarbitrarydercreate(namefileinc)
         elseif (beam(1:8).eq.'gfftwave') then
            test=1
            call derivechamp(nx,ny,nz,nmax,aretecube,FF0,test,xi)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)
!$OMP DO  SCHEDULE(STATIC)            
            do i=1,ndipole
               ii=3*(i-1)
               wrk(i,1)=xi(ii+1)
               wrk(i,2)=xi(ii+2)
               wrk(i,3)=xi(ii+3)
c               write(*,*) 'x',i,xi(ii+1),xi(ii+2),xi(ii+3)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL
            
            test=2
            call derivechamp(nx,ny,nz,nmax,aretecube,FF0,test,xi)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)
!$OMP DO  SCHEDULE(STATIC)
            do i=1,ndipole
               ii=3*(i-1)
               wrk(i,4)=xi(ii+1)
               wrk(i,5)=xi(ii+2)
               wrk(i,6)=xi(ii+3)
c               write(*,*) 'y',i,xi(ii+1),xi(ii+2),xi(ii+3)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL

            test=3
            call derivechamp(nx,ny,nz,nmax,aretecube,FF0,test,xi)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ii)
!$OMP DO  SCHEDULE(STATIC)
            do i=1,ndipole
               ii=3*(i-1)
               wrk(i,7)=xi(ii+1)
               wrk(i,8)=xi(ii+2)
               wrk(i,9)=xi(ii+3)
c               write(*,*) 'z',i,xi(ii+1),xi(ii+2),xi(ii+3)
            enddo
!$OMP ENDDO 
!$OMP END PARALLEL            
         endif
         write(*,*) 'rrrrf',nproche
         if (nproche.eq.-1) then
            call calculforcefft(nx,ny,nz,nx2,ny2,nz2,nxy2,nxm,nym,nzm
     $           ,tabdip,vectx,vecty,vectz,FF,planb)
         else
            call calculforcefftopt(nx,ny,nz,nx2,ny2,nz2,nxy2,nxm,nym,nzm
     $           ,vectx,vecty,vectz,FF,planb)
         endif
         write(*,*) 'rrrrf'
         
c     do i=1,nx2*ny2*nz2
c     write(99,*) 'vect',vectx(i),vecty(i),vectz(i)
c     enddo
c     Compute the x derivative of the local field
         test=1            
         
         call derivativefield2(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
     $        ,ntotal,ntotalm,nmax,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,test,FF0,planb,planf)
         write(*,*) 'rrrrffff'
         do i=1,ndipole
            ii=3*i
            if (Tabdip(i).ne.0) then
               j=Tabdip(i)  
               if (beam(1:11).eq.'pwavelinear') then
                  call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
     $                 ,phi,test,Eder)  
               elseif(beam(1:13).eq.'pwavecircular') then
                  call  ondecirced(xs(j),ys(j),zs(j),k0,E0,ss,theta,phi
     $                 ,test,Eder)
               elseif (beam(1:15).eq.'wavelinearmulti') then
                  call ondeplanedmulti(xs(j),ys(j),zs(j),k0,E0m,ssm,ppm
     $                 ,thetam,phim,nbinc,test,Eder,nstop,infostr)
               elseif (beam(1:7).eq.'antenna') then
                  call dipoleincder(xdip,ydip,zdip,theta,phi,xs(j),ys(j)
     $                 ,zs(j),aretecube,k0,E0,Eder,test,nstop,infostr)
               elseif(beam(1:11).eq.'gwavelinear') then
                  call  gaussianchampd(xs(j),ys(j),zs(j),xgaus,ygaus
     $                 ,zgaus,theta,phi,w0,k0,ss,pp,E0,Eder,tol,nloin)
                  wrk(j,4)=Eder(1,2)
                  wrk(j,5)=Eder(2,2)
                  wrk(j,6)=Eder(3,2)
                  wrk(j,7)=Eder(1,3)
                  wrk(j,8)=Eder(2,3)
                  wrk(j,9)=Eder(3,3)               
               elseif(beam(1:13).eq.'gwavecircular') then
                  call  gaussianchampdcirc(xs(j),ys(j),zs(j),xgaus,ygaus
     $                 ,zgaus,theta,phi,w0,k0,ss,E0,Eder,tol,nloin)
                  wrk(j,4)=Eder(1,2)
                  wrk(j,5)=Eder(2,2)
                  wrk(j,6)=Eder(3,2)
                  wrk(j,7)=Eder(1,3)
                  wrk(j,8)=Eder(2,3)
                  wrk(j,9)=Eder(3,3)          
               elseif(beam(1:15).eq.'gparawavelinear') then
                  call  gaussianparalineard(xs(j),ys(j),zs(j),xgaus
     $                 ,ygaus,zgaus,theta,phi,w0,k0,ss,pp,E0,Eder,nstop
     $                 ,infostr)
                  wrk(j,4)=Eder(1,2)
                  wrk(j,5)=Eder(2,2)
                  wrk(j,6)=Eder(3,2)
                  wrk(j,7)=Eder(1,3)
                  wrk(j,8)=Eder(2,3)
                  wrk(j,9)=Eder(3,3)          
               elseif(beam(1:17).eq.'gparawavecircular') then
                  call  gaussianparacircd(xs(j),ys(j),zs(j),xgaus,ygaus
     $                 ,zgaus,theta,phi,w0,k0,ss,E0,Eder,nstop,infostr)
                  wrk(j,4)=Eder(1,2)
                  wrk(j,5)=Eder(2,2)
                  wrk(j,6)=Eder(3,2)
                  wrk(j,7)=Eder(1,3)
                  wrk(j,8)=Eder(2,3)
                  wrk(j,9)=Eder(3,3)
               elseif (beam(1:8).eq.'gfftwave') then
                  Eder(1,1)=wrk(j,1)
                  Eder(2,1)=wrk(j,2)
                  Eder(3,1)=wrk(j,3)
               elseif (beam(1:9).eq.'arbitrary') then
                  call incidentarbitraryder1(test,xs(j),ys(j),zs(j)
     $                 ,Eder)
               endif 
               wrk(j,1)=Eder(1,1)+FF0(ii-2)
               wrk(j,2)=Eder(2,1)+FF0(ii-1)
               wrk(j,3)=Eder(3,1)+FF0(ii)            

            endif
         enddo

         test=2
         call derivativefield2(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
     $        ,ntotal,ntotalm,nmax,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,test,FF0,planb,planf)
         do i=1,ndipole
            ii=3*i
            if (Tabdip(i).ne.0) then
               j=Tabdip(i)    
               if (beam(1:11).eq.'pwavelinear') then
                  call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
     $                 ,phi,test,Eder) 
               elseif(beam(1:13).eq.'pwavecircular') then
                  call  ondecirced(xs(j),ys(j),zs(j),k0,E0,ss,theta,phi
     $                 ,test,Eder)
               elseif (beam(1:15).eq.'wavelinearmulti') then
                  call ondeplanedmulti(xs(j),ys(j),zs(j),k0,E0m,ssm,ppm
     $                 ,thetam,phim,nbinc,test,Eder,nstop,infostr)
               elseif (beam(1:7).eq.'antenna') then
                  call dipoleincder(xdip,ydip,zdip,theta,phi,xs(j),ys(j)
     $                 ,zs(j),aretecube,k0,E0,Eder,test,nstop,infostr)
               elseif(beam(1:11).eq.'gwavelinear') then
                  Eder(1,2)=wrk(j,4)
                  Eder(2,2)=wrk(j,5)
                  Eder(3,2)=wrk(j,6)         
               elseif(beam(1:13).eq.'gwavecircular') then
                  Eder(1,2)=wrk(j,4)
                  Eder(2,2)=wrk(j,5)
                  Eder(3,2)=wrk(j,6)    
               elseif(beam(1:15).eq.'gparawavelinear') then
                  Eder(1,2)=wrk(j,4)
                  Eder(2,2)=wrk(j,5)
                  Eder(3,2)=wrk(j,6) 
               elseif(beam(1:17).eq.'gparawavecircular') then
                  Eder(1,2)=wrk(j,4)
                  Eder(2,2)=wrk(j,5)
                  Eder(3,2)=wrk(j,6)
               elseif (beam(1:8).eq.'gfftwave') then
                  Eder(1,2)=wrk(j,4)
                  Eder(2,2)=wrk(j,5)
                  Eder(3,2)=wrk(j,6)
               elseif (beam(1:9).eq.'arbitrary') then
                  call incidentarbitraryder1(test,xs(j),ys(j),zs(j)
     $                 ,Eder)
               endif 
               j=j+nmax
               wrk(j,1)=Eder(1,2)+FF0(ii-2)
               wrk(j,2)=Eder(2,2)+FF0(ii-1)
               wrk(j,3)=Eder(3,2)+FF0(ii)
            endif
         enddo
c     Compute the z derivative of the local field
         
         test=3
         call derivativefield2(aretecube,k0,nx,ny,nz,nx2,ny2,nz2,nxy2
     $        ,ntotal,ntotalm,nmax,FFTTENSORxx,FFTTENSORxy,FFTTENSORxz
     $        ,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty,vectz
     $        ,test,FF0,planb,planf)
         do i=1,ndipole
            ii=3*i
            if (Tabdip(i).ne.0) then
               j=Tabdip(i) 
               if (beam(1:11).eq.'pwavelinear') then
                  call  ondeplaned(xs(j),ys(j),zs(j),k0,E0,ss,pp,theta
     $                 ,phi,test,Eder)   
c     write(99,*) 'coucou',xs(j),ys(j),zs(j),k0,E0,ss,pp
c     $                 ,theta,phi,test,Eder,j,i,ndipole
               elseif(beam(1:13).eq.'pwavecircular') then
                  call  ondecirced(xs(j),ys(j),zs(j),k0,E0,ss,theta,phi
     $                 ,test,Eder)
               elseif (beam(1:15).eq.'wavelinearmulti') then
                  call ondeplanedmulti(xs(j),ys(j),zs(j),k0,E0m,ssm,ppm
     $                 ,thetam,phim,nbinc,test,Eder,nstop,infostr)
               elseif (beam(1:7).eq.'antenna') then
                  call dipoleincder(xdip,ydip,zdip,theta,phi,xs(j),ys(j)
     $                 ,zs(j),aretecube,k0,E0,Eder,test,nstop,infostr)
               elseif(beam(1:11).eq.'gwavelinear') then
                  Eder(1,3)=wrk(j,7)
                  Eder(2,3)=wrk(j,8)
                  Eder(3,3)=wrk(j,9)
               elseif(beam(1:13).eq.'gwavecircular') then
                  Eder(1,3)=wrk(j,7)
                  Eder(2,3)=wrk(j,8)
                  Eder(3,3)=wrk(j,9)
               elseif(beam(1:15).eq.'gparawavelinear') then
                  Eder(1,3)=wrk(j,7)
                  Eder(2,3)=wrk(j,8)
                  Eder(3,3)=wrk(j,9)       
               elseif(beam(1:17).eq.'gparawavecircular') then
                  Eder(1,3)=wrk(j,7)
                  Eder(2,3)=wrk(j,8)
                  Eder(3,3)=wrk(j,9)
               elseif (beam(1:8).eq.'gfftwave') then
                  Eder(1,3)=wrk(j,7)
                  Eder(2,3)=wrk(j,8)
                  Eder(3,3)=wrk(j,9)
               elseif (beam(1:9).eq.'arbitrary') then
c     call  ondeplaned(xs(j),ys(j),zs(j),k0,uncomp,0.d0,1.d0
c     $                 ,0.d0,0.d0,test,Eder) 
c     write(*,*) 'Ederz1',Eder
                  call incidentarbitraryder1(test,xs(j),ys(j),zs(j)
     $                 ,Eder)
c     write(*,*) 'Ederz2',Eder
               endif 
               j=j+2*nmax
               wrk(j,1)=Eder(1,3)+FF0(ii-2)
               wrk(j,2)=Eder(2,3)+FF0(ii-1)
               wrk(j,3)=Eder(3,3)+FF0(ii)
            endif
         enddo

c     comptute the optical force from the derivative of the local field
c     and the dipole
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) 
         do i=1,nbsphere3
            FF0(i)=dconjg(FF(i))
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         forcet(1)=0.d0
         forcet(2)=0.d0
         forcet(3)=0.d0

         do i=1,nbsphere 
c     write(*,*) 'ooooooooooooo',i,nbsphere
            k=3*(i-1)
            forcex(i)=0.5d0*dreal(FF0(k+1)*wrk(i,1)+FF0(k+2) *wrk(i,2)
     $           +FF0(k+3)*wrk(i,3))*quatpieps0       
c            write(*,*) 'x',wrk(i,1),wrk(i,2),wrk(i,3),i
            j=i+nmax
            forcey(i)=0.5d0*dreal(FF0(k+1)*wrk(j,1)+FF0(k+2) *wrk(j,2)
     $           +FF0(k+3)*wrk(j,3))*quatpieps0       
c            write(*,*) 'y',wrk(j,1),wrk(j,2),wrk(j,3)
            j=i+2*nmax
            forcez(i)=0.5d0*dreal(FF0(k+1)*wrk(j,1)+FF0(k+2) *wrk(j,2)
     $           +FF0(k+3)*wrk(j,3))*quatpieps0 
c            write(*,*) 'z',wrk(j,1),wrk(j,2),wrk(j,3)
            forcet(1)=forcet(1)+forcex(i)
            forcet(2)=forcet(2)+forcey(i)
            forcet(3)=forcet(3)+forcez(i)
         enddo
  
c     write(*,*) 'coucou'
         if (nstop == -1) then
            infostr = 'Calculation cancelled at the optical force!'
            return
         endif
c     save the density of force
c     subunit=0
         if (nforced.eq.1.and.nmat.eq.0) then
            do i=1,ndipole
c     subunit= subunit+1
               k=tabdip(i)
               if (k.ne.0) then
                  write(60,*) forcex(k)
                  write(61,*) forcey(k)
                  write(62,*) forcez(k)
               else
                  write(60,*) 0.d0
                  write(61,*) 0.d0
                  write(62,*) 0.d0
               endif
            enddo
c     save the density of the optical force
            close(60)
            close(61)
            close(62)            
         endif

         forcem=dsqrt(forcet(1)*forcet(1)+forcet(2)*forcet(2)+forcet(3)
     $        *forcet(3))
         write(99,*) 'optical force x',forcet(1)
         write(99,*) 'optical force y',forcet(2)
         write(99,*) 'optical force z',forcet(3)
         write(*,*) 'Optical force x      : ',forcet(1)
         write(*,*) 'Optical force y      : ',forcet(2)
         write(*,*) 'Optical force z      : ',forcet(3)
         write(*,*) 'Modulus of the force : ',forcem
c     write(*,*) 'optical force x',forcet(1)
c     write(*,*) 'optical force y',forcet(2)
c     write(*,*) 'optical force z',forcet(3)
c     forcemie=(MIECext-GSCA*MIECsca)/8.d0/pi

         if (nsection.eq.1.or.nsectionsca.eq.1) then
            if (object(1:6).eq.'sphere') then
               forcemie=(MIECEXT-GSCA*MIECSCA)/8.d0 /pi*I0*quatpieps0
               write(*,*) 'Force with Mie',forcemie
               write(*,*) 'Error with Mie in %',100.d0*(forcemie-forcem)
     $              /forcemie
            endif
         endif

         
         write(99,*) 'Modulus of the force',forcem
         write(*,*) ' ********** End optical force *********'
         write(*,*) ' '
      endif 

      if (ntorque*nforce.eq.1) then
         write(99,*) ' Compute the optical torque '
         write(*,*) '********** Compute the optical torque *******'
c     Compute the optical torque on the particle 
         couplet(1)=0.d0
         couplet(2)=0.d0
         couplet(3)=0.d0
c     compute the center of gravity of the object
         xg=0.d0
         yg=0.d0
         zg=0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:xg,yg,zg)         
         do i=1,nbsphere
            xg=xg+xs(i)
            yg=yg+ys(i)
            zg=zg+zs(i)
         enddo
!$OMP ENDDO 
!$OMP END PARALLEL
         
         xg=xg/dble(nbsphere)
         yg=yg/dble(nbsphere)
         zg=zg/dble(nbsphere)
         write(99,*) 'position of the center of gravity',xg,yg,zg
         write(*,*) 'Position of the center of gravity :'
         write(*,*) 'xg=',xg
         write(*,*) 'yg=',yg
         write(*,*) 'zg=',zg
c     gamma=r x F+1/2 Re(p^* x p / alpha_0)

         do i=1,nbsphere
            kk=3*(i-1)

            do ii=1,3
               Em(ii)=2.d0/3.d0*icomp*k03*FF(kk+ii)+FFloc(kk+ii)
            enddo

            torquex(i)=0.5d0*dreal(FF0(kk+2)*Em(3)-FF0(kk+3)*Em(2))
     $           *quatpieps0+((ys(i)-yg)*forcez(i)-(zs(i)-zg)*forcey(i))
            torquey(i)=0.5d0*dreal(-FF0(kk+1)*Em(3)+FF0(kk+3)*Em(1))
     $           *quatpieps0+((zs(i)-zg)*forcex(i)-(xs(i)-xg)*forcez(i))
            torquez(i)=0.5d0*dreal(FF0(kk+1)*Em(2)-FF0(kk+2)*Em(1))
     $           *quatpieps0+((xs(i)-xg)*forcey(i)-(ys(i)-yg)*forcex(i))
            couplet(1)=couplet(1)+torquex(i)
            couplet(2)=couplet(2)+torquey(i)
            couplet(3)=couplet(3)+torquez(i)               
         enddo

c     save the density of torque
         if (nstop == -1) then
            infostr = 'Calculation cancelled at the optical torque!'
            return
         endif
c     subunit=0
         if (ntorqued.eq.1.and.nmat.eq.0) then
            do i=1,ndipole
c     subunit= subunit+1
               k=tabdip(i)
               if (k.ne.0) then
                  write(63,*) torquex(k)
                  write(64,*) torquey(k)
                  write(65,*) torquez(k)
               else
                  write(63,*) 0.d0
                  write(64,*) 0.d0
                  write(65,*) 0.d0                                    
               endif
            enddo
c     save the density of optical torque
            close(63)
            close(64)
            close(65)            
         endif

         couplem=dsqrt(couplet(1)*couplet(1)+couplet(2)*couplet(2)
     $        +couplet(3)*couplet(3))
         write(99,*) 'Optical torque x',couplet(1)
         write(99,*) 'Optical torque y',couplet(2)
         write(99,*) 'Optical torque z',couplet(3)
         write(*,*) 'optical torque x',couplet(1)
         write(*,*) 'optical torque y',couplet(2)
         write(*,*) 'optical torque z',couplet(3)
         write(99,*) 'Modulus of the optical torque',couplem
         write(*,*) 'Modulus of the optical torque',couplem

         if (nsection.eq.1.or.nsectionsca.eq.1) then
            if (object(1:6).eq.'sphere' .and. beam(1:13).eq
     $           .'pwavecircular') then
               couplemie=MIECABS/8.d0/k0/pi*I0 *quatpieps0
               write(*,*) 'Torque with Mie',couplemie
               write(*,*) 'Error with Mie in %',100.d0*(couplemie
     $              -couplem)/couplemie
            endif
         endif
         
         write(*,*) '******* End optical torque ********'
         write(*,*) ' '
      endif
c     calcul sur des forces et couple sur differents objets si presents
      if (numberobjet.ne.1) then
         forcexmulti=0.d0
         forceymulti=0.d0
         forcezmulti=0.d0
         torquexmulti=0.d0
         torqueymulti=0.d0
         torquezmulti=0.d0    
         do i=1,nbsphere
            is=tabmulti(i)
            forcexmulti(is)=forcexmulti(is)+forcex(i)
            forceymulti(is)=forceymulti(is)+forcey(i)
            forcezmulti(is)=forcezmulti(is)+forcez(i)
            torquexmulti(is)=torquexmulti(is)+torquex(i)
            torqueymulti(is)=torqueymulti(is)+torquey(i)
            torquezmulti(is)=torquezmulti(is)+torquez(i)            
         enddo
      endif
      if (numberobjet.ge.2) then
         do is=1,numberobjet
            write(*,*) 'Optical force for the',is,'object'
            write(*,*) 'optical force x',forcexmulti(is)
            write(*,*) 'optical force y',forceymulti(is)
            write(*,*) 'optical force z',forcezmulti(is)
            write(*,*) 'Optical torque for the',is,'object'
            write(*,*) 'optical torque x',torquexmulti(is)
            write(*,*) 'optical torque y',torqueymulti(is)
            write(*,*) 'optical torque z',torquezmulti(is)
         enddo
      endif

c     output file
      close(99)
 9999 format(512(d22.15,1x))
c$$$  write(*,*) 'prout',  lambda,beam,object,trope, materiaumulti,nnnr
c$$$  $     ,tolinit,methodeit,polarizability, nquad,nlecture,filereread,
c$$$  c     output file cdm.out
c$$$  $     nlocal,nmacro,nsection,nsectionsca,nquickdiffracte,nrig,
c$$$  $     nforce,nforced,ntorque,ntorqued,nproche,nlentille,nquicklens
c$$$  $     ,nobjet,
c$$$  c     cube, sphere (includes multiple)
c$$$  $     side, sidex, sidey, sidez, hauteur,
c$$$  $     numberobjet, rayonmulti, xgmulti, ygmulti, zgmulti,
c$$$  $     epsmulti, epsanimulti,lc,hc,ng,
c$$$  c     ellipsoid+arbitrary
c$$$  $     demiaxea,demiaxeb,demiaxec,thetaobj,phiobj,psiobj,
c$$$  $     namefileobj,
c$$$  c     planewavecircular.in / planewavelinear.in files
c$$$  $     theta, phi, pp, ss, P0, w0, xgaus, ygaus, zgaus,namefileinc,
c$$$  c     return info stringf
c$$$  $     infostr, nstop,
c$$$  c     return scalar results
c$$$  $     nbsphere, ndipole, aretecube,
c$$$  $     lambda10n, k0, tol1, ncompte, nloop,
c$$$  $     Cext, Cabs, Csca, Cscai, gasym, irra, E0,
c$$$  $     forcet, forcem,
c$$$  $     couplet, couplem,
c$$$  $     nxm, nym, nzm
c$$$  write(*,*) 'prout1'
c$$$  do i=1,nxm*nym*nzm
c$$$  write(*,*) 'i',i,incidentfield(i), localfield(i),
c$$$  $        macroscopicfield(i),xs(i), ys(i), zs(i), xswf(i), yswf(i),
c$$$  $        zswf(i),forcex(i),forcey(i),forcez(i), incidentfieldx(i),
c$$$  $        incidentfieldy(i), incidentfieldz(i), localfieldx(i),
c$$$  $        localfieldy(i), localfieldz(i), macroscopicfieldx(i),
c$$$  $        macroscopicfieldy(i), macroscopicfieldz(i), torquex(i)
c$$$  $        ,torquey(i),torquez(i)
c$$$  
c$$$  enddo
c$$$  write(*,*) 'prout2',ntheta, nphi, forcexmulti ,forceymulti
c$$$  $     ,forcezmulti,torquexmulti,torqueymulti, to rquezmulti, nfft2d
c$$$  $     ,numaper ,gross
      if (nstop == -1) then
         infostr = 'Calculation cancelled at the end!'
         return
      endif

      call date_and_time(date,time,zone,valuesf)
      write(*,*) 'Real time used by the code',(valuesf(5) -valuesi(5))
     $     *3600.d0+(valuesf(6)-valuesi(6))*60.d0+(valuesf(7)
     $     -valuesi(7)) +(valuesf(8)-valuesi(8))/1000.d0
      
      write(*,*) 'COMPLETED!'
      infostr='COMPLETED!'
      write(*,*) 'END'
      end

