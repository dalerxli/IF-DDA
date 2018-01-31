      subroutine objetsphererandom(trope,eps,epsani,eps0,xs,ys,zs,k0
     $     ,aretecube,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym
     $     ,nzm,methode,epsilon ,polarisa,sidex,sidey,sidez,xg,yg,zg
     $     ,xsphe,ysphe,zsphe,radius,density,nseed,na,nmat,infostr
     $     ,nstop)
      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k,ks
     $     ,test,IP(3),nnnr,dddis,inv,na,nstop,nxm,nym,nzm,nseed,ngraine
     $     ,is,ns,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),xsphe(nmax)
     $     ,ysphe(nmax),zsphe(nmax),k0,xg,yg,zg,x,y,z,xt,yt,zt,eps0
     $     ,sunif,aretecube,sidex,sidey,sidez,side,pi,radius,density
     $     ,dist ,volsphere
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp

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
      pi=dacos(-1.d0)
      nseed=0
      ngraine=4*nseed+1

      sidex=sidex*1.d-9
      sidey=sidey*1.d-9
      sidez=sidez*1.d-9
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      radius=radius*1.d-9
      
c     diminue nnnr si pas bon facteur
 100  call factor235(nnnr,IP,test)        
      if (test.eq.1) then
         nnnr=nnnr-1
         goto 100
      endif
      side=max(sidex,sidey,sidez)
      aretecube=side/dble(nnnr)
      write(*,*) 'a',aretecube
c     calcul nx,ny,nz optimum
      nx=int(sidex/aretecube)
      ny=int(sidey/aretecube)
      nz=int(sidez/aretecube)
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

      
c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

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
      
      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box : Box too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif

c     calcul nb de sphere d'apres density
      if (density.ge.0.5d0.and.density.le.2.d0) then
         nstop=1
         infostr='density too large!'
         return
      endif
      volsphere=4.d0*pi/3.d0*radius*radius*radius
      ns=int(density*sidex*sidey*sidez/volsphere)

      if (ns.le.1) then 
         nstop=1
         infostr='density too small!'
         return
      endif
      if (density.gt.2.d0) then
         ns=int(density)
         density=dble(ns)*volsphere/sidex/sidey/sidez
      endif
      
      ks=0
      do i=1,ns*10
c     tirage aleatoire dans le boite               
         
         xt=-sidex/2.d0+radius+(sidex-2.d0*radius)*SUNIF(ngraine)
         yt=-sidey/2.d0+radius+(sidey-2.d0*radius)*SUNIF(ngraine)
         zt=-sidez/2.d0+radius+(sidez-2.d0*radius)*SUNIF(ngraine)
         
c     test si intersection
         if (ks.ge.1) then
            do j=1,ks
               dist=dsqrt((xt-xsphe(j))*(xt-xsphe(j))+(yt-ysphe(j))*(yt
     $              -ysphe(j))+(zt-zsphe(j))*(zt-zsphe(j)))
               if (dist.lt.radius*2.d0) goto 10
            enddo
         endif
         ks=ks+1
         xsphe(ks)=xt
         ysphe(ks)=yt
         zsphe(ks)=zt
         if (ks.eq.ns) goto 20
         
 10   enddo

 20   write(*,*) 'Number of sphere',ks,na,nx,ny,nz
      
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
                  
                  test=0
                  do is=1,ks
                     if ((x-xsphe(is))**2+(y-ysphe(is))**2+(z-zsphe(is))
     $                    **2.le.radius**2) then
                        test=is
                     endif
                  enddo
                  if (test.ne.0) then
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
                        call polaepstens(aretecube,epsani,eps0,k0
     $                       ,dddis,methode,inv,polaeps)
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

                  test=0
                  do is=1,ks
                     if ((x-xsphe(is))**2+(y-ysphe(is))**2+(z-zsphe(is))
     $                    **2.le.radius**2) then
                        test=is
                     endif
                  enddo
                  
                  if (test.ne.0) then                       
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
c     write(*,*) 'cubo',nbsphere,ctmp,eps
                     else
                        call polaepstens(aretecube,epsani,eps0,k0 ,dddis
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

      write(*,*) 'Random sphere::nbsphere    =',nbsphere,ndipole
      write(*,*) 'density of sphere obtained =',dble(ks)*volsphere/sidex
     $     /sidey/sidez
     
      end
