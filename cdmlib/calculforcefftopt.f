      subroutine calculforcefftopt(nx,ny,nz,nx2,ny2,nz2,nxy2,nxm,nym,nzm
     $     ,vectx,vecty,vectz,FF,planb)
      implicit none
      integer nx,ny,nz,nx2,ny2,nz2,nxy2,nxm,nym,nzm
      double complex, dimension(8*nxm*nym*nzm) :: vectx,vecty,vectz
      double complex, dimension(3*nxm*nym*nzm) :: FF
      
      integer i,j,k,ii,jj,indice
      integer*8 planb

      
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO  SCHEDULE(STATIC)
      do i=1,8*nx*ny*nz
         vectx(i)=0.d0
         vecty(i)=0.d0
         vectz(i)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,indice,jj)
!$OMP DO  SCHEDULE(STATIC) COLLAPSE(3)      
      do k=1,nz
         do j=1,ny
            do i=1,nx
               indice=i+nx2*(j-1)+nxy2*(k-1)               
               jj=3*(i+nx*(j-1)+nx*ny*(k-1))
               vectx(indice)=FF(jj-2)
               vecty(indice)=FF(jj-1)
               vectz(indice)=FF(jj)
            enddo
         enddo
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL
      
      call dfftw_execute_dft(planb, vectx, vectx)   
      call dfftw_execute_dft(planb, vecty, vecty)     
      call dfftw_execute_dft(planb, vectz, vectz)

      end
