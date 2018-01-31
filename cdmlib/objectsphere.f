      subroutine objetsphere(trope,eps,epsani,eps0,xs,ys,zs,k0,aretecube
     $     ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,methode,na
     $     ,epsilon,polarisa,rayon,xg,yg,zg,nmat,infostr,nstop)

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,na,nstop,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,xg,yg,zg,x,y,z,eps0
     $     ,aretecube
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp
      
      double precision rayon
      character*2 methode
      character*3 trope
      character(64) infostr

c     Initialization
      nbsphere=0
      ndipole=0 
      Tabdip=0
      polarisa=0.d0
      epsilon=0.d0
      dddis=1
      inv=1
      
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      rayon=rayon*1.d-9
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      if (rayon.eq.0.d0) then
         nstop=1
         infostr='object sphere: radius=0!'
         return
      endif

c     verfie si on est bien multiple de 2 3 5 pour la discretisation,
c     car ma FFT est basee sur une decomposition en nombre premier de 2
c     3 5. Si on utilisait FFTW ceci disparaitrait.
 100  call factor235(nnnr,IP,test)        
      if (test.eq.1) then
         nnnr=nnnr-1
         goto 100
      endif
      nx=nnnr
      ny=nnnr
      nz=nnnr
c     size of the subunit
      aretecube=2.d0*rayon/dble(nnnr)

      if (na.eq.-1) then
         do i=1,nnnr         
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(20,*) x+xg

                  ndipole=ndipole+1
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then
                     nbsphere=nbsphere+1
                     Tabdip(ndipole)=nbsphere
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     if (trope.eq.'iso') then
                        call poladiff(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else
                        call polaepstens(aretecube,epsani,eps0,k0,dddis
     $                       ,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
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
         do i=1,nnnr         
            do j=1,nnnr              
               do k=1,nnnr                
                  x=-rayon+aretecube*(dble(k)-0.5d0)
                  y=-rayon+aretecube*(dble(j)-0.5d0)
                  z=-rayon+aretecube*(dble(i)-0.5d0)

                  if (j.eq.1.and.k.eq.1.and.nmat.eq.0) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmat.eq.0) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1.and.nmat.eq.0) write(20,*) x+xg

                  ndipole=ndipole+1
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  if (dsqrt(x*x+y*y+z*z).lt.rayon) then                    
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     if (trope.eq.'iso') then
                        call poladiff(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
c                        write (*,*) 'nbsphere = ', nbsphere
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else
                        call polaepstens(aretecube,epsani,eps0,k0,dddis
     $                       ,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
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
      write (*,*) ' OBJECT SPHERE FINISHED'
      end
