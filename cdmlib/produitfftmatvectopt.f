c     FFB=FFB-A*FFX==> FFB=XI et FFX=XI*pola==> FFB=(I-A*pola)*XI
c     routine optimiser pour le cube.
      subroutine produitfftmatvectopt(FFTTENSORxx,FFTTENSORxy
     $     ,FFTTENSORxz,FFTTENSORyy,FFTTENSORyz,FFTTENSORzz,vectx,vecty
     $     ,vectz,ntotalm,ntotal,nmax,nxm,nym,nzm,nx,ny,nz,nx2,ny2,nxy2
     $     ,nz2,FFX,FFB)
      implicit none
      integer jj,i,j,k,indice,ntotal,ntotalm,nxm,nym,nzm,nx,ny,nz
     $     ,nmax,nx2,ny2,nxy2,nz2
      double complex FFTTENSORxx(ntotalm),FFTTENSORxy(ntotalm)
     $     ,FFTTENSORxz(ntotalm),FFTTENSORyy(ntotalm)
     $     ,FFTTENSORyz(ntotalm),FFTTENSORzz(ntotalm),FFX(3*nmax),FFB(3
     $     *nmax),vectx(ntotalm) ,vecty(ntotalm),vectz(ntotalm)
      double complex ctmpx,ctmpy,ctmpz
      integer nxy
c      double precision t1,t2

      nxy=nx*ny
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(indice) 
!$OMP DO SCHEDULE(STATIC)      
      do indice=1,ntotal
         vectx(indice)=0.d0
         vecty(indice)=0.d0
         vectz(indice)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(jj,i,j,k,indice) 
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)       
      do k=1,nz
         do j=1,ny
            do i=1,nx
c     position du dipole
               indice=i+nx2*(j-1)+nxy2*(k-1)
               jj=3*(i+nx*(j-1)+nxy*(k-1))
               vectx(indice)=FFX(jj-2)
               vecty(indice)=FFX(jj-1)
               vectz(indice)=FFX(jj)
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
!$OMP PARALLEL  DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION   
      CALL ZFFT3D(vectx,NX2,NY2,NZ2,1)
!$OMP SECTION      
      CALL ZFFT3D(vecty,NX2,NY2,NZ2,1)
!$OMP SECTION      
      CALL ZFFT3D(vectz,NX2,NY2,NZ2,1)
!$OMP END SECTIONS
!$OMP END PARALLEL 
      
!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(indice,ctmpx,ctmpy,ctmpz)   
!$OMP DO SCHEDULE(STATIC)   
      do indice=1,ntotal
         ctmpx=vectx(indice)
         ctmpy=vecty(indice)
         ctmpz=vectz(indice)
         vectx(indice)=FFTTENSORxx(indice)*ctmpx+FFTTENSORxy(indice)
     $        *ctmpy+ FFTTENSORxz(indice)*ctmpz
         vecty(indice)=FFTTENSORxy(indice)*ctmpx+FFTTENSORyy(indice)
     $        *ctmpy+ FFTTENSORyz(indice)*ctmpz         
         vectz(indice)=FFTTENSORxz(indice)*ctmpx+FFTTENSORyz(indice)
     $        *ctmpy+ FFTTENSORzz(indice)*ctmpz
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL

      
c     FFT inverse (deconvolution)
!$OMP PARALLEL  DEFAULT(SHARED)
!$OMP SECTIONS 
!$OMP SECTION        
      CALL ZFFT3D(vectx,nx2,ny2,nz2,-1)
!$OMP SECTION      
      CALL ZFFT3D(vecty,nx2,ny2,nz2,-1)
!$OMP SECTION      
      CALL ZFFT3D(vectz,nx2,ny2,nz2,-1)
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(indice,jj,i,j,k)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)               
      do k=1,nz
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)+nxy2*(k-1)
               jj=3*(i+nx*(j-1)+nxy*(k-1))
               FFB(jj-2)=FFB(jj-2)-vectx(indice)
               FFB(jj-1)=FFB(jj-1)-vecty(indice)
               FFB(jj)=FFB(jj)-vectz(indice)
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      end
   

c*************************************************
