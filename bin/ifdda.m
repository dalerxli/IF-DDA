close all
clear all

load inputmatlab.mat -ascii

nproche=inputmatlab(1);         % Defined box of computation for the near field
nlocal=inputmatlab(2);          % Compute the local field
nmacro=inputmatlab(3);          % Compute the macroscopic field
nsection=inputmatlab(4);        % Compute the cross section
nsectionsca=inputmatlab(5);     % Compute  C_sca, Poynting and g
nquickdiffracte=inputmatlab(6); % Compute  C_sca, Poynting and g with FFT
nforce=inputmatlab(7);          % Compute the optical force
nforced=inputmatlab(8);         % Compute the optical force
ntorque=inputmatlab(9);         % Compute the optical force
ntorqued=inputmatlab(10);       % Compute the optical force
nlentille=inputmatlab(11);      % Compute the object through the microscope
nquicklens=inputmatlab(12);     % Compte the lens with FFT
nphi=inputmatlab(13);           % nphi
ntheta=inputmatlab(14);         % ntheta
niso=inputmatlab(15);           % 0 isotrope, 1 anisotrope
nfft=inputmatlab(16);           % size of the FFT
k0=inputmatlab(17);             % Wavenumber
numaper=inputmatlab(18);        % Numerical aperture
nprochefft=inputmatlab(19);     % =1 if wide field
nobjet=inputmatlab(20);         % =1 do ony the objet

  icomp=complex(0,1);
  
load x.mat -ascii
load y.mat -ascii
load z.mat -ascii

nx=max(size(x));
ny=max(size(y));
nz=max(size(z));

%%%%%%%%%%%%%%%% Begin plot dipole %%%%%%%%%%%%%%%%%%%%%%%%%%

load xc.mat -ascii
load yc.mat -ascii
load zc.mat -ascii

  
load epsilon.mat -ascii

if (niso == 0);
if (nproche == -1 );

figure(1)
set(1,'DefaultAxesFontName','Times')
set(1,'DefaultAxesFontSize',12)
set(1,'DefaultAxesFontWeight','Bold')
set(1,'DefaultTextfontName','Times')
set(1,'DefaultTextfontSize',12)
set(1,'DefaultTextfontWeight','Bold')
set(1,'Position',[0 600 500 500])


m=0;n=max(size( epsilon)); 
for i=1:n; if epsilon(i,1)~= 1 || epsilon(i,2)~= 0;
m=m+1;epsilonbg(m)=epsilon(i,1);end;end;

  
scatter3(xc,yc,zc,10,epsilonbg)
axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title('Position of the dipoles')

if nobjet == 1; return ;end;

 else

figure(1)
set(1,'DefaultAxesFontName','Times')
set(1,'DefaultAxesFontSize',12)
set(1,'DefaultAxesFontWeight','Bold')
set(1,'DefaultTextfontName','Times')
set(1,'DefaultTextfontSize',12)
set(1,'DefaultTextfontWeight','Bold')
set(1,'Position',[0 600 1200 500])

  subplot(1,2,1)
     
scatter3(xc,yc,zc,10,epsilon(:,1))
axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title('Position of the dipoles')
subplot(1,2,2)

m=0;n=max(size( epsilon)); 
for i=1:n; if epsilon(i,1)~= 1 || epsilon(i,2)~= 0 ; m=m+1;
xcc(m)=xc(i);ycc(m)=yc(i);zcc(m)=zc(i);epsilonbg(m)=epsilon(i,1);
end;end;


  
scatter3(xcc,ycc,zcc,10,epsilonbg)
axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title({'Position of the dipoles' 'without background'})
end;

else

figure(1)
set(1,'DefaultAxesFontName','Times')
set(1,'DefaultAxesFontSize',12)
set(1,'DefaultAxesFontWeight','Bold')
set(1,'DefaultTextfontName','Times')
set(1,'DefaultTextfontSize',12)
set(1,'DefaultTextfontWeight','Bold')
set(1,'Position',[0 600 1200 500])

 
     
  scatter3(xc,yc,zc,10)
  axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title('Position of the dipoles (anisotropic object)')
if nobjet == 1; return ;end;
end;
%%%%%%%%%%%%%%%% End plot dipole %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin plot epsilon %%%%%%%%%%%%%%%%%%%%%%%%%%



if niso == 0;


matxyepsilonr=reshape(epsilon(:,1),nx,ny,nz);
matxyepsiloni=reshape(epsilon(:,2),nx,ny,nz);

clear epsilon

figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style','text','Fontsize',16,'Fontweight','bold',...
'Position',[380 440 300 50],'String','Plot relative permittivity');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 400 170 50],...
'Callback',{@plotepsilon,nx,ny,nz,x,y,z,matxyepsilonr,matxyepsiloni});

else

for i=1:nz;for j=1:ny;for k=1:nx;kk=9*(k+nx*(j-1)+nx*ny*(i-1)-1);
matepsilonrxx(i,j,k)=epsilon(kk+1,1);
matepsilonixx(i,j,k)=epsilon(kk+1,2);
matepsilonrxy(i,j,k)=epsilon(kk+2,1);
matepsilonixy(i,j,k)=epsilon(kk+2,2);
matepsilonrxz(i,j,k)=epsilon(kk+3,1);
matepsilonixz(i,j,k)=epsilon(kk+3,2);
matepsilonryx(i,j,k)=epsilon(kk+4,1);
matepsiloniyx(i,j,k)=epsilon(kk+4,2);
matepsilonryy(i,j,k)=epsilon(kk+5,1);
matepsiloniyy(i,j,k)=epsilon(kk+5,2);
matepsilonryz(i,j,k)=epsilon(kk+6,1);
matepsiloniyz(i,j,k)=epsilon(kk+6,2);
matepsilonrzx(i,j,k)=epsilon(kk+7,1);
matepsilonizx(i,j,k)=epsilon(kk+7,2);
matepsilonrzy(i,j,k)=epsilon(kk+8,1);
matepsilonizy(i,j,k)=epsilon(kk+8,2);
matepsilonrzz(i,j,k)=epsilon(kk+9,1);
matepsilonizz(i,j,k)=epsilon(kk+9,2);
end;end;end;
clear epsilon

  
figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[300 440 380 50],'String','Plot tensor of permittivity in (x,y)-plane for        component');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
	  {'xx','xy','xz','yx','yy','yz','zx','zy','zz'},...
'Position',[490 425 40 40],...
	  'Callback',{@plotepsilonani,nx,ny,nz,x,y,z,matepsilonrxx,matepsilonixx,matepsilonrxy,matepsilonixy,matepsilonrxz,matepsilonixz,matepsilonryx,matepsiloniyx,matepsilonryy,matepsiloniyy,matepsilonryz,matepsiloniyz,matepsilonrzx,matepsilonizx,matepsilonrzy,matepsilonizy,matepsilonrzz,matepsilonizz});
		   
end;
%%%%%%%%%%%%%%%% End plot epsilon %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin incident field %%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nprochefft == 0) 

load incidentfield.mat -ascii
load incidentfieldx.mat -ascii
load incidentfieldy.mat -ascii
load incidentfieldz.mat -ascii


matxyincifield=reshape(incidentfield,nx,ny,nz);
matxyincifieldx=reshape(incidentfieldx(:,1)+icomp*incidentfieldx(:,2),nx,ny,nz);
matxyincifieldy=reshape(incidentfieldy(:,1)+icomp*incidentfieldy(:,2),nx,ny,nz);
matxyincifieldz=reshape(incidentfieldz(:,1)+icomp*incidentfieldz(:,2),nx,ny,nz);

clear incidentfieldx
clear incidentfieldy
clear incidentfieldz


figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])




uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 210 20],'String','Plot incident field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
'Callback',{@plotincifield,nx,ny,nz,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz});

 else

   
load xwf.mat -ascii
load ywf.mat -ascii
load zwf.mat -ascii

nxm=max(size(xwf));
nym=max(size(ywf));
nzm=max(size(zwf));

load incidentfieldwf.mat -ascii
load incidentfieldxwf.mat -ascii
load incidentfieldywf.mat -ascii
load incidentfieldzwf.mat -ascii


matxyincifield=reshape(incidentfieldwf,nxm,nym,nzm);
matxyincifieldx=reshape(incidentfieldxwf(:,1)+icomp*incidentfieldxwf(:,2),nxm,nym,nzm);
matxyincifieldy=reshape(incidentfieldywf(:,1)+icomp*incidentfieldywf(:,2),nxm,nym,nzm);
matxyincifieldz=reshape(incidentfieldzwf(:,1)+icomp*incidentfieldzwf(:,2),nxm,nym,nzm);

clear incidentfieldxwf
clear incidentfieldywf
clear incidentfieldzwf

figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])




uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 210 20],'String','Plot incident field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
'Callback',{@plotincifield,nx,ny,nz,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz});

end;



%%%%%%%%%%%%%%%% End incident field %%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Plot local field %%%%%%%%%%%%%%%%%%%%%%%%%%

  
if nlocal == 1;

if (nprochefft ==0) 

load localfield.mat -ascii
load localfieldx.mat -ascii
load localfieldy.mat -ascii
load localfieldz.mat -ascii


matxylocalfield=reshape(localfield,nx,ny,nz);
matxylocalfieldx=reshape(localfieldx(:,1)+icomp*localfieldx(:,2),nx,ny,nz);
matxylocalfieldy=reshape(localfieldy(:,1)+icomp*localfieldy(:,2),nx,ny,nz);
matxylocalfieldz=reshape(localfieldz(:,1)+icomp*localfieldz(:,2),nx,ny,nz);

clear localfieldx
clear localfieldy
clear localfieldz

figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 200 20],'String','Plot local field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
'Callback',{@plotlocalfield,nx,ny,nz,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz});


 else
   
load xwf.mat -ascii
load ywf.mat -ascii
load zwf.mat -ascii

nxm=max(size(xwf));
nym=max(size(ywf));
nzm=max(size(zwf));

load localfieldwf.mat -ascii
load localfieldxwf.mat -ascii
load localfieldywf.mat -ascii
load localfieldzwf.mat -ascii


  
matxylocalfield=reshape(localfieldwf,nxm,nym,nzm);
matxylocalfieldx=reshape(localfieldxwf(:,1)+icomp*localfieldxwf(:,2),nxm,nym,nzm);
matxylocalfieldy=reshape(localfieldywf(:,1)+icomp*localfieldywf(:,2),nxm,nym,nzm);
matxylocalfieldz=reshape(localfieldzwf(:,1)+icomp*localfieldzwf(:,2),nxm,nym,nzm);

clear localfieldxwf
clear localfieldywf
clear localfieldzwf

figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 200 20],'String','Plot local field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
'Callback',{@plotlocalfield,nxm,nym,nzm,xwf,ywf,zwf,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz});


end;


end;

%%%%%%%%%%%%%%%% End local field %%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%% Begin macrocopic field %%%%%%%%%%%%%%%%%%%%%%%%%%
if nmacro == 1;


if (nprochefft ==0) 

load macroscopicfield.mat -ascii
load macroscopicfieldx.mat -ascii
load macroscopicfieldy.mat -ascii
load macroscopicfieldz.mat -ascii


matxymacrofield=reshape(macroscopicfield,nx,ny,nz);
matxymacrofieldx=reshape(macroscopicfieldx(:,1)+icomp*macroscopicfieldx(:,2),nx,ny,nz);
matxymacrofieldy=reshape(macroscopicfieldy(:,1)+icomp*macroscopicfieldy(:,2),nx,ny,nz);
matxymacrofieldz=reshape(macroscopicfieldz(:,1)+icomp*macroscopicfieldz(:,2),nx,ny,nz);

clear macroscopicfieldx
clear macroscopicfieldy
clear macroscopicfieldz

figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 550 200 50],'String','Plot macroscopic field')

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
'Callback',{@plotmacrofield,nx,ny,nz,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz});

 else
load xwf.mat -ascii
load ywf.mat -ascii
load zwf.mat -ascii

nxm=max(size(xwf));
nym=max(size(ywf));
nzm=max(size(zwf));

load macroscopicfieldwf.mat -ascii
load macroscopicfieldxwf.mat -ascii
load macroscopicfieldywf.mat -ascii
load macroscopicfieldzwf.mat -ascii


  
matxymacrofield=reshape(macroscopicfieldwf,nxm,nym,nzm);
matxymacrofieldx=reshape(macroscopicfieldxwf(:,1)+icomp*macroscopicfieldxwf(:,2),nxm,nym,nzm);
matxymacrofieldy=reshape(macroscopicfieldywf(:,1)+icomp*macroscopicfieldywf(:,2),nxm,nym,nzm);
matxymacrofieldz=reshape(macroscopicfieldzwf(:,1)+icomp*macroscopicfieldzwf(:,2),nxm,nym,nzm);

clear macrofieldxwf
clear macrofieldywf
clear macrofieldzwf

figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 550 200 50],'String','Plot macroscopic field')

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
'Callback',{@plotmacrofield,nxm,nym,nzm,xwf,ywf,zwf,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz});

end;



end;
%%%%%%%%%%%%%%%% End macrocopic field %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin Poynting vector %%%%%%%%%%%%%%%%%%%%%%%%%%
if nsectionsca == 1;
if nquickdiffracte == 0;

load poynting.mat -ascii

ray=reshape(poynting,nphi,ntheta);
ray(nphi+1,:)=ray(1,:);
  
for i=1:ntheta;for j=1:nphi+1;
theta=pi*(i-1)/(ntheta-1);
phi=2*pi*(j-1)/(nphi);
xpo(j,i)=ray(j,i)*sin(theta)*cos(phi);
ypo(j,i)=ray(j,i)*sin(theta)*sin(phi);
zpo(j,i)=ray(j,i)*cos(theta);
end;
end;


figure(100)
set(100,'DefaultAxesFontName','Times')
set(100,'DefaultAxesFontSize',12)
set(100,'DefaultAxesFontWeight','Bold')
set(100,'DefaultTextfontName','Times')
set(100,'DefaultTextfontSize',12)
set(100,'DefaultTextfontWeight','Bold')
set(100,'Position',[1000 0 600 600])
surf(xpo,ypo,zpo,ray)
shading interp

title('Poynting vector')

xlabel('x')
ylabel('y')
zlabel('z')


else


load poynting.mat -ascii

ray=reshape(poynting,nphi,ntheta);
ray(nphi+1,:)=ray(1,:);
  
for i=1:ntheta;for j=1:nphi+1;
theta=pi*(i-1)/(ntheta-1);
phi=2*pi*(j-1)/(nphi);
xpo(j,i)=ray(j,i)*sin(theta)*cos(phi);
ypo(j,i)=ray(j,i)*sin(theta)*sin(phi);
zpo(j,i)=ray(j,i)*cos(theta);
end;
end;


figure(100)
set(100,'DefaultAxesFontName','Times')
set(100,'DefaultAxesFontSize',12)
set(100,'DefaultAxesFontWeight','Bold')
set(100,'DefaultTextfontName','Times')
set(100,'DefaultTextfontSize',12)
set(100,'DefaultTextfontWeight','Bold')
set(100,'Position',[1000 0 600 600])

surf(xpo,ypo,zpo,ray)
shading interp

title('Poynting vector interpolated in 3D')

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)
zlabel('$z$','Interpreter','latex','Fontsize',18)



load poyntingpos.mat -ascii
load poyntingneg.mat -ascii
load kx.mat -ascii
load ky.mat -ascii


nxp=max(size(kx));
nyp=max(size(ky));



figure(101)
set(101,'DefaultAxesFontName','Times')
set(101,'DefaultAxesFontSize',12)
set(101,'DefaultAxesFontWeight','Bold')
set(101,'DefaultTextfontName','Times')
set(101,'DefaultTextfontSize',12)
set(101,'DefaultTextfontWeight','Bold')
set(101,'Position',[1000 0 1000 600])

suptitle('Poynting Modulus in $k_x$ and $k_y$ plane','Interpreter','latex','Fontsize',18)

subplot(1,2,1)

imagesc(kx/k0,ky/k0,reshape(poyntingpos,nxp,nyp)')
axis xy

hold on
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

axis image
axis equal
title('$k_z>0$','Interpreter','latex','Fontsize',18)

xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
colorbar

subplot(1,2,2)

imagesc(kx/k0,ky/k0,reshape(poyntingneg,nxp,nyp)')
axis xy

hold on
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

axis image
axis equal
title('$k_z<0$','Interpreter','latex','Fontsize',18)

xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

colorbar
  
end;

end;
%%%%%%%%%%%%%%%% End Poynting vector %%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%% Begin optical force %%%%%%%%%%%%%%%%%%%%%%%%%%
if nforced==1;

load forcex.mat -ascii
load forcey.mat -ascii
load forcez.mat -ascii



matxyforcex=reshape(forcex,nx,ny,nz);
matxyforcey=reshape(forcey,nx,ny,nz);
matxyforcez=reshape(forcez,nx,ny,nz);

clear forcex
clear forcey
clear forcez

for i=1:nz;for j=1:ny;for k=1:nx;
xx(k,j,i)=x(k);yy(k,j,i)=y(j);zz(k,j,i)=z(i);
end;end;end;

figure(200)
set(200,'DefaultAxesFontName','Times')
set(200,'DefaultAxesFontSize',12)
set(200,'DefaultAxesFontWeight','Bold')
set(200,'DefaultTextfontName','Times')
set(200,'DefaultTextfontSize',12)
set(200,'DefaultTextfontWeight','Bold')
set(200,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[350 570 300 25],'String','Plot Optical force');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotforce,nx,ny,nz,x,y,z,xx,yy,zz,matxyforcex,matxyforcey,matxyforcez});



figure(201)
set(201,'DefaultAxesFontName','Times')
set(201,'DefaultAxesFontSize',12)
set(201,'DefaultAxesFontWeight','Bold')
set(201,'DefaultTextfontName','Times')
set(201,'DefaultTextfontSize',12)
set(201,'DefaultTextfontWeight','Bold')
set(201,'Position',[0 0 1000 600])

quiver3(xx,yy,zz,matxyforcex,matxyforcey,matxyforcez)


xlabel('x')
ylabel('y')
zlabel('z')
title('Density of optical force')

end;  
%%%%%%%%%%%%%%%% End  optical force %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin optical torque %%%%%%%%%%%%%%%%%%%%%%%%%%
if ntorqued==1;


load torquex.mat -ascii
load torquey.mat -ascii
load torquez.mat -ascii



matxytorquex=reshape(torquex,nx,ny,nz);
matxytorquey=reshape(torquey,nx,ny,nz);
matxytorquez=reshape(torquez,nx,ny,nz);

clear torquex
clear torquey
clear torquez

for i=1:nz;for j=1:ny;for k=1:nx;
xx(k,j,i)=x(k);yy(k,j,i)=y(j);zz(k,j,i)=z(i);
end;end;end;

figure(300)
set(300,'DefaultAxesFontName','Times')
set(300,'DefaultAxesFontSize',12)
set(300,'DefaultAxesFontWeight','Bold')
set(300,'DefaultTextfontName','Times')
set(300,'DefaultTextfontSize',12)
set(300,'DefaultTextfontWeight','Bold')
set(300,'Position',[0 0 1000 600])


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[350 570 300 25],'String','Plot Optical torque');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plottorque,nx,ny,nz,x,y,z,xx,yy,zz,matxytorquex,matxytorquey,matxytorquez});



figure(301)
set(301,'DefaultAxesFontName','Times')
set(301,'DefaultAxesFontSize',12)
set(301,'DefaultAxesFontWeight','Bold')
set(301,'DefaultTextfontName','Times')
set(301,'DefaultTextfontSize',12)
set(301,'DefaultTextfontWeight','Bold')
set(301,'Position',[0 0 1000 600])

quiver3(xx,yy,zz,matxytorquex,matxytorquey,matxytorquez)


xlabel('x')
ylabel('y')
zlabel('z')
title('Density of optical torque')

end;  
%%%%%%%%%%%%%%%%%%%%% End  optical torque %%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%% Microscopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nlentille == 1;

load fourier.mat -ascii
load fourierx.mat -ascii
load fouriery.mat -ascii
load fourierz.mat -ascii
load kxfourier.mat -ascii

nnfft=max(size(kxfourier));

fourierm=reshape(fourier,nnfft,nnfft);
fourierxc=reshape(fourierx(:,1)+icomp*fourierx(:,2),nnfft,nnfft);
fourieryc=reshape(fouriery(:,1)+icomp*fouriery(:,2),nnfft,nnfft);
fourierzc=reshape(fourierz(:,1)+icomp*fourierz(:,2),nnfft,nnfft);

clear fourier
clear fourierx
clear fouriery
clear fourierz

load fourierinc.mat -ascii
load fourierincx.mat -ascii
load fourierincy.mat -ascii
load fourierincz.mat -ascii

fourierincm=reshape(fourierinc,nnfft,nnfft);
fourierincxc=reshape(fourierincx(:,1)+icomp*fourierincx(:,2),nnfft,nnfft);
fourierincyc=reshape(fourierincy(:,1)+icomp*fourierincy(:,2),nnfft,nnfft);
fourierinczc=reshape(fourierincz(:,1)+icomp*fourierincz(:,2),nnfft,nnfft);

clear fourierinc
clear fourierincx
clear fourierincy
clear fourierincz


load image.mat -ascii
load imagex.mat -ascii
load imagey.mat -ascii
load imagez.mat -ascii

imagem=reshape(image,nfft,nfft);
imagexc=reshape(imagex(:,1)+icomp*imagex(:,2),nfft,nfft);
imageyc=reshape(imagey(:,1)+icomp*imagey(:,2),nfft,nfft);
imagezc=reshape(imagez(:,1)+icomp*imagez(:,2),nfft,nfft);

clear image
clear imagex
clear imagey
clear imagez


load imageinc.mat -ascii
load imageincx.mat -ascii
load imageincy.mat -ascii
load imageincz.mat -ascii

imageincm=reshape(imageinc,nfft,nfft);
imageincxc=reshape(imageincx(:,1)+icomp*imageincx(:,2),nfft,nfft);
imageincyc=reshape(imageincy(:,1)+icomp*imageincy(:,2),nfft,nfft);
imageinczc=reshape(imageincz(:,1)+icomp*imageincz(:,2),nfft,nfft);

clear imageinc
clear imageincx
clear imageincy
clear imageincz

load ximage.mat -ascii

%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier %%%%%%%%%%%%%%%%%%%%%%%%%

figure(400)

set(400,'DefaultAxesFontName','Times')
set(400,'DefaultAxesFontSize',12)
set(400,'DefaultAxesFontWeight','Bold')
set(400,'DefaultTextfontName','Times')
set(400,'DefaultTextfontSize',12)
set(400,'DefaultTextfontWeight','Bold')
set(400,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])
  
imagesc(kxfourier,kxfourier,fourierm')
axis xy  
axis equal
axis image

hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

caxis([ 0 max(max(fourierm))])
shading interp
axis equal
axis image
colorbar
  
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[350 575 100 20],'String','Fourier:')
uicontrol('Style','text','Fontsize',12,'Fontweight','bold','Position',[350 545 200 18],'String','Numerical aperture:')
uicontrol('Style', 'text','Fontsize',12,'Fontweight','bold', 'String', num2str(numaper),'Position', [540 545 40 18]);

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Modulus','x-component','y-component','z-component'},...
'Position', [460 566 150 30],...
'Callback',{@plotfourier,numaper,kxfourier,fourierm,fourierxc,fourieryc,fourierzc});



%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier +incident %%%%%%%%%%%%%%%%%%%%%%%%%

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])
  
imagesc(kxfourier,kxfourier,fourierincm')
axis xy  
axis equal
axis image

hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

caxis([ 0 max(max(fourierincm))])
shading interp
axis equal
axis image
colorbar
  
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[240 575 220 20],'String','Fourier+incident:')
uicontrol('Style','text','Fontsize',12,'Fontweight','bold','Position',[350 545 200 18],'String','Numerical aperture:')
uicontrol('Style', 'text','Fontsize',12,'Fontweight','bold', 'String', num2str(numaper),'Position', [540 545 40 18]);

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Modulus','x-component','y-component','z-component'},...
'Position', [460 566 150 30],...
'Callback',{@plotfourierinc,numaper,kxfourier,fourierincm,fourierincxc,fourierincyc,fourierinczc});



%%%%%%%%%%%%%%%%% Image  %%%%%%%%%%%%%%%%%%%%%%


figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,imagem')
axis xy
caxis([ 0 max(max(imagem))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[350 560 100 30],'String','Image:')

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
'Callback',{@plotimage,ximage,imagem,imagexc,imageyc,imagezc});

%%%%%%%%%%%%%%%%% Image + incident %%%%%%%%%%%%%%%%%%%%%%



figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,imageincm')
axis xy
caxis([min(min(imageincm)) max(max(imageincm))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[250 560 200 30],'String','Image+incident:')

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
'Callback',{@plotimageinc,ximage,imageincm,imageincxc,imageincyc,imageinczc});

  
 end;
