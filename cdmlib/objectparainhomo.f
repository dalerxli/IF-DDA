      subroutine objetparainhomo(eps,eps0,xs,ys,zs,k0,aretecube
     $     ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym,nzm
     $     ,methode,epsilon ,polarisa,sidex,sidey,sidez,xg,yg,zg,lc,hc
     $     ,ng,epsb,na,nmat,infostr,nstop)
      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,i,j,k ,test
     $     ,IP(3),nnnr,dddis,inv,na,nstop,nxm,nym,nzm,ng,ngraine,nk,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,xg,yg,zg,x,y,z,eps0
     $     ,aretecube,sidex,sidey,sidez,side,pi,lc,hc
      double complex icomp,eps,polarisa(nmax,3 ,3) ,epsilon(nmax,3,3)
     $     ,ctmp,epsb(nmax)

      double precision kx,ky,kz,dkx,dky,dkz,coeff1,coeff2,coeff3 ,phase
     $     ,spec,moyenne,ecartype,lx,ly,lz,sunif
      
      character*2 methode
      character(64) infostr
      write(*,*) 'object para inhomo'
c     Initialization
      nbsphere=0
      ndipole=0 
      Tabdip=0
      polarisa=0.d0
      epsilon=0.d0
      dddis=1
      inv=1
      pi=dacos(-1.d0)
      icomp=(0.d0,1.d0)
      ngraine=4*ng+1

      
 100  call factor235(nnnr,IP,test)        
      if (test.eq.1) then
         nnnr=nnnr-1
         goto 100
      endif

c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      sidex=sidex*1.d-9
      sidey=sidey*1.d-9
      sidez=sidez*1.d-9
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      lc=lc*1.d-9

      write(*,*) 'object4',sidex,sidey,sidez
      if (sidex.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidex=0!'
         return
      elseif (sidey.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidey=0!'
         return
      elseif (sidez.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidez=0!'
         return
      endif
      write(*,*) 'object'


      
      side=max(sidex,sidey,sidez)
      aretecube=side/dble(nnnr)


c     boite au plus pres de l'objet
      nx=int(sidex/aretecube)
      ny=int(sidey/aretecube)
      nz=int(sidez/aretecube)
      write(*,*) 'nnn',nx,ny,nz
 200  call factor235(nx,IP,test)        
      if (test.eq.1) then
         nx=nx+1
         goto 200
      endif
 300  call factor235(ny,IP,test)        
      if (test.eq.1) then
         ny=ny+1
         goto 300
      endif
 400  call factor235(nz,IP,test)        
      if (test.eq.1) then
         nz=nz+1
         goto 400
      endif
      write(*,*) 'nnn',nx,ny,nz

      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box : Box too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif



      lx=dble(nx)*aretecube
      ly=dble(ny)*aretecube
      lz=dble(nz)*aretecube
      dkx=2.d0*pi/lx
      dky=2.d0*pi/ly
      dkz=2.d0*pi/lz
      nk=0
      
      coeff1=0.125d0*hc*hc*lc*lc*lc/(pi*dsqrt(pi))
      coeff2=lc*lc
      coeff3=4.d0*dkx*dky*dkz
      do i=1,nz
         do j=1,ny
            do k=1,nx
               nk=nk+1
               if (k.ge.1 .and. k.le.(nx/2+1)) then
                  kx=dble(k-1)*dkx
               else
                  kx=dble(k-1-nx)*dkx
               endif
               if (j.ge.1 .and. j.le.(ny/2+1)) then
                  ky=dble(j-1)*dky
               else
                  ky=dble(j-1-ny)*dky
               endif
               if (i.ge.1 .and. i.le.(nz/2+1)) then
                  kz=dble(i-1)*dkz
               else
                  kz=dble(i-1-nz)*dkz
               endif
c     calcul du spec
               spec=coeff1*dexp(-0.25d0*coeff2*(kx*kx+ky*ky+kz*kz))
c     phase aleatoire
               phase=2.d0*pi*SUNIF(ngraine)
c     Amplitude complexe
               epsb(nk)=dsqrt(spec*coeff3)*cdexp(icomp *phase)
c     Spectre symetrique pour diffraction harmonique
               if (nk.eq.1) epsb(nk)=0.d0
               if (nk.ge.nx/2+2.and.nk.le.nx) epsb(nk)=0.d0
               if (nk.ge.(ny/2+1)*nx+1.and.nk.le.nx*ny) epsb(nk)=0.d0
               if (nk.ge.(nz/2+1)*nx*ny+1.and.nk.le.nx*ny*nz) epsb(nk)
     $              =0.d0

            enddo
         enddo
      enddo     
      
c     Profil des hauteurs
      call ZFFT3D(epsb,NX,NY,NZ,1)
      moyenne=0.d0
      ecartype=0.d0
      do i=1,nk
         epsb(i)=dreal(epsb(i))+eps
         moyenne=moyenne+dreal(epsb(i))
         ecartype=ecartype+cdabs(epsb(i))**2.d0      
      enddo
      
      
      moyenne=moyenne/dble(nk)
      ecartype=ecartype/dble(nk)
      
      if (na.eq.-1) then

         do i=1,nz
            do j=1,ny
               do k=1,nx

                  x=-sidex/2.d0+aretecube*(dble(k)-0.5d0)
                  y=-sidey/2.d0+aretecube*(dble(j)-0.5d0)
                  z=-sidez/2.d0+aretecube*(dble(i)-0.5d0)                 

                  if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(20,*) x+xg

                  ndipole=ndipole+1                   
                  
                  if (dabs(x).le.sidex/2.d0.and.dabs(y).le.sidey
     $                 /2.d0.and.dabs(z).le.sidez/2.d0) then  

                     nbsphere=nbsphere+1
                     Tabdip(ndipole)=nbsphere
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg                   
                     eps=epsb(ndipole)
                     
                     call poladiff(aretecube,eps,eps0,k0,dddis ,methode
     $                    ,ctmp)  
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps
                     if (nmat.eq.0) then
                        write(10,*) xs(nbsphere)
                        write(11,*) ys(nbsphere)
                        write(12,*) zs(nbsphere)
                     endif
                  endif
               enddo
            enddo
         enddo 
      elseif (na.eq.0) then
         

c     fin essai
         do i=1,nz
            do j=1,ny
               do k=1,nx

                  x=-sidex/2.d0+aretecube*(dble(k)-0.5d0)
                  y=-sidey/2.d0+aretecube*(dble(j)-0.5d0)
                  z=-sidez/2.d0+aretecube*(dble(i)-0.5d0)                 

                  

                  if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(20,*) x+xg

                  ndipole=ndipole+1  
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  

                  if (dabs(x).le.sidex/2.d0.and.dabs(y).le.sidey
     $                 /2.d0.and.dabs(z).le.sidez/2.d0) then
     $                 
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     eps=epsb(ndipole)
                     call poladiff(aretecube,eps,eps0,k0,dddis
     $                    ,methode,ctmp)  
                     polarisa(nbsphere,1,1)=ctmp
                     polarisa(nbsphere,2,2)=ctmp
                     polarisa(nbsphere,3,3)=ctmp
                     epsilon(nbsphere,1,1)=eps
                     epsilon(nbsphere,2,2)=eps
                     epsilon(nbsphere,3,3)=eps

                     if (nmat.eq.0) then
                        write(10,*) xs(nbsphere)
                        write(11,*) ys(nbsphere)
                        write(12,*) zs(nbsphere)
                     endif
                  else
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg     
                     epsilon(nbsphere,1,1)=(1.d0,0.d0)
                     epsilon(nbsphere,2,2)=(1.d0,0.d0)
                     epsilon(nbsphere,3,3)=(1.d0,0.d0)
                     if (nmat.eq.0) then
                        write(10,*) xs(nbsphere)
                        write(11,*) ys(nbsphere)
                        write(12,*) zs(nbsphere)
                     endif
                  endif
               enddo
            enddo
         enddo
      else 
         infostr='na should be equal to -1 or 0'
         nstop=1
         return
      endif

      if (ndipole.gt.nmax) then
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif

      close(10)
      close(11)
      close(12)
      close(20)
      close(21)
      close(22)

      write(*,*) 'Cuboid inhomogeneous::nbsphere=',nbsphere
      end
