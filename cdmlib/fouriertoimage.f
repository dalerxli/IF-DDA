      subroutine fouriertoimage(deltakx,deltaky,Eimagex,Eimagey,Eimagez
     $     ,Eimageincx,Eimageincy,Eimageincz,nfft2D,nfft2d2)
      
      implicit none
      integer nfft2D,nfft2d2
      double precision deltakx,deltaky
      double complex Eimagex(nfft2d*nfft2d),Eimagey(nfft2d*nfft2d)
     $     ,Eimagez(nfft2d*nfft2d),Eimageincx(nfft2d*nfft2d)
     $     ,Eimageincy(nfft2d *nfft2d),Eimageincz(nfft2d*nfft2d)

      integer i,j,indicex,indicey,indice,kk
      double precision tmp
      double complex ctmp
      
      tmp=deltakx*deltaky

c     Compute the FFT inverse to get the image through the microscope
!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP SECTIONS 
!$OMP SECTION            
      CALL ZFFT2D(Eimageincx,nfft2D,nfft2D,1)
!$OMP SECTION         
      CALL ZFFT2D(Eimageincy,nfft2D,nfft2D,1)
!$OMP SECTION          
      CALL ZFFT2D(Eimageincz,nfft2D,nfft2D,1)
!$OMP SECTION            
      CALL ZFFT2D(Eimagex,nfft2D,nfft2D,1)
!$OMP SECTION         
      CALL ZFFT2D(Eimagey,nfft2D,nfft2D,1)
!$OMP SECTION          
      CALL ZFFT2D(Eimagez,nfft2D,nfft2D,1)         
!$OMP END SECTIONS
!$OMP END PARALLEL
      
      
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,indicex,indicey,indice,kk,ctmp)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)           
      do i=-nfft2d2,nfft2d2-1
         do j=-nfft2d2,-1
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
            kk=i+nfft2d2+1+nfft2d*(j+nfft2d2)
            
            ctmp=Eimagex(kk)
            Eimagex(kk)=Eimagex(indice)*tmp
            Eimagex(indice)=ctmp*tmp
            ctmp=Eimagey(kk)
            Eimagey(kk)=Eimagey(indice)*tmp
            Eimagey(indice)=ctmp*tmp
            ctmp=Eimagez(kk)
            Eimagez(kk)=Eimagez(indice)*tmp
            Eimagez(indice)=ctmp*tmp

            ctmp=Eimageincx(kk)
            Eimageincx(kk)=Eimageincx(indice)*tmp
            Eimageincx(indice)=ctmp*tmp
            ctmp=Eimageincy(kk)
            Eimageincy(kk)=Eimageincy(indice)*tmp
            Eimageincy(indice)=ctmp*tmp
            ctmp=Eimageincz(kk)
            Eimageincz(kk)=Eimageincz(indice)*tmp
            Eimageincz(indice)=ctmp*tmp
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL


      end
