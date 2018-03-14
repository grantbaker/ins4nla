function void=ins();
% This program solves the incompressible Navier-Stokes equations on
% the domain [x,y] \in [0,Lx] \times [0,Ly] using a grid spacing
% h=Lx/Nx, Ly/Ny.

Lx=1;
Ly=1;
Nx=200;
Ny=200;
Nsteps=1000;
k=0.01;

alpha = 0.1/k;
nu    = 0.001;
mms   = 3;
timemet = 3;
pl = 1;

t  = 0;
hx = Lx/Nx;
hy = Ly/Ny;
% Set up the grid, we include the ghost points.
x=(-hx:hx:Lx+hx)';
y=(-hy:hy:Ly+hy)';
[Y,X]=meshgrid(y,x);

% We store u,v,p as two dimensional arrays
u=zeros(Nx+3,Ny+3);
v=zeros(Nx+3,Ny+3);
p=zeros(Nx+3,Ny+3);
% arrays for forcing.
Fu=zeros(Nx+3,Ny+3);
Fv=zeros(Nx+3,Ny+3);
% setup pressure laplacian
disp('Setting up Laplacians')
tt=cputime;
[LapP]= setupLapP(Nx,Ny,hx,hy);

r = sparse(ones((Nx+1)*(Ny+1),1));
LapPBig = sparse([[LapP r];[r' 0]]);
A = LapPBig;
A(end,end) = 1;
[M1,M2] = ilu(A);
[LapUV]= setupLapUV(Nx,Ny,hx,hy);
LapUV=-0.5*nu*LapUV+speye(size(LapUV))/k;
fprintf('... done setting up Laplacians it took %g seconds\n',cputime-tt)
% These are the arrays for the boundary forcing u=gu, v=gv on the
% boundary. We order them
% gux(1,:) = "left"
% gux(2,:) = "right"
% guy(1,:) = "bottom"
% guy(2,:) = "top"

gux=zeros(2,Ny+1);
guy=zeros(2,Nx+1);
gvx=zeros(2,Ny+1);
gvy=zeros(2,Nx+1);
gtux=zeros(2,Ny+1);
gtuy=zeros(2,Nx+1);
gtvx=zeros(2,Ny+1);
gtvy=zeros(2,Nx+1);

% These hold boundary forcings
pbx=zeros(2,Ny+1);
pby=zeros(2,Nx+1);

% mms mode checks various parts of the discretization
if (mms ==1)
  t=0.;
  u=mmsfun(X,Y,t,1,0,0,0);
  v=mmsfun(X,Y,t,2,0,0,0);
  p=mmsfun(X,Y,t,3,0,0,0);
  u(1:2,:)=0;u(end-1:end,:)=0;u(:,1:2)=0;u(:,end-1:end)=0;
  v(1:2,:)=0;v(end-1:end,:)=0;v(:,1:2)=0;v(:,end-1:end)=0;
  p(1,:)=0;p(end,:)=0;p(:,1)=0;p(:,end)=0;
  [Fu,Fv]                               = computePDEForcing(X,Y,t,mms,nu);
  [gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy] = computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t,mms);
  [u,v]                                 = updateBCforU(u,v,X,Y,t,gux,guy,gvx,gvy,hx,hy);
  fprintf('Interior ERROR u and v: %e %e \n',...
          max(max(abs(u(2:end-1,2:end-1)-mmsfun(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),t,1,0,0,0)))),...
          max(max(abs(v(2:end-1,2:end-1)-mmsfun(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),t,2,0,0,0)))))
  [u,v]                     = computeAndUpdateGPforU(u,v,X,Y,hx,hy);
  [pbx,pby]         = computeGPforP(pbx,pby,u,v,X,Y,t,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fu,Fv);
  [LapP]            = setupLapP(Nx,Ny,hx,hy);
  r=sparse(ones((Nx+1)*(Ny+1),1));
  LapPBig = [[LapP r];[r' 0]];    
  [bP]              = setupRhsideP(Nx,Ny,hx,hy,u,v,Fu,Fv,alpha,pbx,pby,mms,X,Y,t,nu);
  [p]               = computeP(LapP,LapPBig,bP,p,Nx,Ny,mms,t,M1,M2);
  xt=X(2:end-1,2:end-1);
  yt=Y(2:end-1,2:end-1);
  % this tests the boundary conditions and the pressure solver
  e1=u-mmsfun(X,Y,t,1,0,0,0);
  e2=v-mmsfun(X,Y,t,2,0,0,0);
  e1(1,1)=0;e1(1,end)=0;e1(end,1)=0;e1(end,end)=0;
  e2(1,1)=0;e2(1,end)=0;e2(end,1)=0;e2(end,end)=0;
  fprintf('MMS ERROR u, v and p: %e %e %e \n',...
          max(max(abs(e1))),...
          max(max(abs(e2))),...
          max(max(abs(p(2:end-1,2:end-1)-mmsfun(xt,yt,t,3,0,0,0)))))
  % this tests the computation of the leplicit and implicit part of
  % the right hand side.
  [Leu,Lev]=computeLE(u,v,p,X,Y,t,hx,hy);
  [Liu,Liv]=computeLI(u,v,X,Y,t,hx,hy,nu);
  Leumms=-mmsfun(xt,yt,t,1,0,0,0).*mmsfun(xt,yt,t,1,0,1,0)...
         -mmsfun(xt,yt,t,2,0,0,0).*mmsfun(xt,yt,t,1,0,0,1)...
         -mmsfun(xt,yt,t,3,0,1,0);
  Levmms=-mmsfun(xt,yt,t,1,0,0,0).*mmsfun(xt,yt,t,2,0,1,0)...
         -mmsfun(xt,yt,t,2,0,0,0).*mmsfun(xt,yt,t,2,0,0,1)...
         -mmsfun(xt,yt,t,3,0,0,1);
  Liumms=nu*(mmsfun(xt,yt,t,1,0,2,0)+mmsfun(xt,yt,t,1,0,0,2));
  Livmms=nu*(mmsfun(xt,yt,t,2,0,2,0)+mmsfun(xt,yt,t,2,0,0,2));
  e1=(Leu-Leumms);
  e2=(Lev-Levmms);
  e3=(Liu-Liumms);
  e4=(Liv-Livmms);
  e1(1,:)=0;e1(end,:)=0;e1(:,1)=0;e1(:,end)=0;
  e2(1,:)=0;e2(end,:)=0;e2(:,1)=0;e2(:,end)=0;
  e3(1,:)=0;e3(end,:)=0;e3(:,1)=0;e3(:,end)=0;
  e4(1,:)=0;e4(end,:)=0;e4(:,1)=0;e4(:,end)=0;
  fprintf('MMS ERROR Leu, Lev Liu Liv: %e %e %e %e \n',...
          max(max(abs(e1))),max(max(abs(e2))),...
          max(max(abs(e3))),max(max(abs(e4))))
end
if (mms == 2)
  u=rand(size(u))-0.5;
  v=rand(size(v))-0.5;
end

if (mms == 3)
  nv = 500;
  xv = 0.25+0.5*rand(nv,1);
  yv = 0.25+0.5*rand(nv,1);
  ga = 0.5+rand(nv,1);
  rr = 0.02;
  i = 1;
  [ut,vt]=taylor(X,Y,xv(i),yv(i),rr,ga(i));
  u=ut;
  v=vt;
  for i = 2:nv
  [ut,vt]=taylor(X,Y,xv(i),yv(i),rr,ga(i));
  u=u+ut;
  v=v+vt;
  end
end

if (mms == 5)
  u=(Y-0.5).*exp(-((Y-0.5)/0.2).^2);
  v=0*Y;
end


% start up the computation with a single Euler step
uold=0*u;
vold=0*v;
pold=0*p;
up=0*u;
vp=0*v;
[Leu,Lev]=computeLE(u,v,p,X,Y,t,hx,hy);
[Liu,Liv]=computeLI(u,v,X,Y,t,hx,hy,nu);
uold(2:end-1,2:end-1)=u(2:end-1,2:end-1)-k*(Leu+Liu+Fu(2:end-1,2:end-1));
vold(2:end-1,2:end-1)=v(2:end-1,2:end-1)-k*(Lev+Liv+Fv(2:end-1,2:end-1));
[gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy] = computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t-k,mms);
[uold,vold]                           = updateBCforU(uold,vold,X,Y,t-k,gux,guy,gvx,gvy,hx,hy);
[uold,vold]                           = computeAndUpdateGPforU(uold,vold,X,Y,hx,hy);
[Fuold,Fvold]                         = computePDEForcing(X,Y,t-k,mms,nu);
[pbx,pby] = computeGPforP(pbx,pby,uold,vold,X,Y,t-k,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fuold,Fvold);
[bP]      = setupRhsideP(Nx,Ny,hx,hy,uold,vold,Fuold,Fvold,alpha,pbx,pby,mms,X,Y,t-k,nu);
[pold]       = computeP(LapP,LapPBig,bP,pold,Nx,Ny,mms,t-k,M1,M2);

if (mms == 1)
  e1=uold-mmsfun(X,Y,t-k,1,0,0,0);
  e2=vold-mmsfun(X,Y,t-k,2,0,0,0);
  e1(1,1)=0;e1(1,end)=0;e1(end,1)=0;e1(end,end)=0;
  e2(1,1)=0;e2(1,end)=0;e2(end,1)=0;e2(end,end)=0;

  fprintf('Euler startup ERROR u, v and p: %e %e %e \n',max(max(abs(e1))),max(max(abs(e2))),...
          max(max(abs(pold(2:end-1,2:end-1)-mmsfun(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),t-k,3,0,0,0)))))
end
[Leuold,Levold]         = computeLE(uold,vold,pold,X,Y,t-k,hx,hy);
[Liuold,Livold]         = computeLI(uold,vold,X,Y,t-k,hx,hy,nu);

if (mms == 1)
  [Fu,Fv] = computePDEForcing(X,Y,t,mms,nu);
  [Fuold,Fvold] = computePDEForcing(X,Y,t-k,mms,nu);
  uold=mmsfun(X,Y,t-k,1,0,0,0);
  vold=mmsfun(X,Y,t-k,2,0,0,0);
  pold=mmsfun(X,Y,t-k,3,0,0,0);
  u=mmsfun(X,Y,t,1,0,0,0);
  v=mmsfun(X,Y,t,2,0,0,0);
  p=mmsfun(X,Y,t,3,0,0,0);
  [Leuold,Levold]         = computeLE(uold,vold,pold,X,Y,t-k,hx,hy);
  [Liuold,Livold]         = computeLI(uold,vold,X,Y,t-k,hx,hy,nu);
end

disp('Starting Time Loop....')
for nt=1:Nsteps
  tt=cputime;
  if (timemet == 1)
    % Adams-Bashfort;
    % Compute Forcing and right hand side
    [Leu,Lev]=computeLE(u,v,p,X,Y,t,hx,hy);
    [Liu,Liv]=computeLI(u,v,X,Y,t,hx,hy,nu);
    % Compute u^(t+k), v^(t+k)
    u(2:end-1,2:end-1)=u(2:end-1,2:end-1)+1.5*k*(Leu+Liu+Fu(2:end-1,2:end-1))...
        -0.5*k*(Leuold+Liuold+Fuold(2:end-1,2:end-1));
    v(2:end-1,2:end-1)=v(2:end-1,2:end-1)+1.5*k*(Lev+Liv+Fv(2:end-1,2:end-1))...
        -0.5*k*(Levold+Livold+Fvold(2:end-1,2:end-1));
    Liuold=Liu; Livold=Liv; Leuold=Leu; Levold=Lev; Fuold=Fu; Fvold=Fv;
    t=t+k;
    [gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy] = computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t,mms);
    [u,v]                                 = updateBCforU(u,v,X,Y,t,gux,guy,gvx,gvy,hx,hy);
    [u,v]                     = computeAndUpdateGPforU(u,v,X,Y,hx,hy);

    [Fu,Fv]           = computePDEForcing(X,Y,t,mms,nu);
    [pbx,pby]         = computeGPforP(pbx,pby,u,v,X,Y,t,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fu,Fv);
    [bP]              = setupRhsideP(Nx,Ny,hx,hy,u,v,Fu,Fv,alpha,pbx,pby,mms,X,Y,t,nu);
    [p]               = computeP(LapP,LapPBig,bP,p,Nx,Ny,mms,t,M1,M2);

  elseif (timemet == 2)
    % Hyman method
    [Leu,Lev]=computeLE(u,v,p,X,Y,t,hx,hy);
    [Liu,Liv]=computeLI(u,v,X,Y,t,hx,hy,nu);
    Liuold=Liu; Livold=Liv; Leuold=Leu; Levold=Lev; Fuold=Fu; Fvold=Fv;
    % Predict with leapfrog to get up and vp.
    up(2:end-1,2:end-1)=uold(2:end-1,2:end-1)+2*k*(Leu+Liu+Fu(2:end-1,2:end-1));
    vp(2:end-1,2:end-1)=vold(2:end-1,2:end-1)+2*k*(Lev+Liv+Fv(2:end-1,2:end-1));
    t=t+k;
    [gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy] = computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t,mms);
    [up,vp]                               = updateBCforU(up,vp,X,Y,t,gux,guy,gvx,gvy,hx,hy);
    [up,vp]                               = computeAndUpdateGPforU(up,vp,X,Y,hx,hy);

    [Fu,Fv]           = computePDEForcing(X,Y,t,mms,nu);
    [pbx,pby]         = computeGPforP(pbx,pby,up,vp,X,Y,t,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fu,Fv);
    [bP]              = setupRhsideP(Nx,Ny,hx,hy,up,vp,Fu,Fv,alpha,pbx,pby,mms,X,Y,t,nu);
    [p]               = computeP(LapP,LapPBig,bP,p,Nx,Ny,mms,t,M1,M2);

    % These are for the predictor right hand side
    [Leu,Lev]=computeLE(up,vp,p,X,Y,t,hx,hy);
    [Liu,Liv]=computeLI(up,vp,X,Y,t,hx,hy,nu);
    % Store u and v at time n to stick into uold and vold below
    up=u;
    vp=v;
    % Corrector
    u(2:end-1,2:end-1)=(4/5)*u(2:end-1,2:end-1)...
        +(1/5)*uold(2:end-1,2:end-1)...
        +(2/5)*k*(Leu+Liu+Fu(2:end-1,2:end-1)...
                  +2*(Leuold+Liuold+Fuold(2:end-1,2:end-1)));
    v(2:end-1,2:end-1)=(4/5)*v(2:end-1,2:end-1)...
        +(1/5)*vold(2:end-1,2:end-1)...
        +(2/5)*k*(Lev+Liv+Fv(2:end-1,2:end-1)...
                  +2*(Levold+Livold+Fvold(2:end-1,2:end-1)));
    % update pressure anb boundary conditions at time n+1
    [gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy] = computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t,mms);
    [u,v]                                 = updateBCforU(u,v,X,Y,t,gux,guy,gvx,gvy,hx,hy);
    [u,v]                                 = computeAndUpdateGPforU(u,v,X,Y,hx,hy);
    [Fu,Fv]           = computePDEForcing(X,Y,t,mms,nu);
    [pbx,pby]         = computeGPforP(pbx,pby,u,v,X,Y,t,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fu,Fv);
    [bP]              = setupRhsideP(Nx,Ny,hx,hy,u,v,Fu,Fv,alpha,pbx,pby,mms,X,Y,t,nu);
    [p]               = computeP(LapP,LapPBig,bP,p,Nx,Ny,mms,t,M1,M2);
    % update old values
    uold=up;
    vold=vp;
  elseif (timemet == 3)
    % semi-implicit method
    % Get boundary conditions at time n+1 ;
    t=t+k;
    [gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy] = computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t,mms);
    [bU,bV] = setupRhsideUV(Nx,Ny,hx,hy,gux,guy,gvx,gvy);
    [u,v,uold,vold] = computeUV(Nx,Ny,bU,bV,Leu,Lev,Fu,Fv,...
                                Leuold,Levold,Liu,Liv,Fuold,Fvold,u,v,uold,vold,k,nu,LapUV);
    % update boundary conditions and compute pressure
    [u,v]                                 = updateBCforU(u,v,X,Y,t,gux,guy,gvx,gvy,hx,hy);
    [u,v]                                 = computeAndUpdateGPforU(u,v,X,Y,hx,hy);
    Fuold=Fu; Fvold=Fv;
    [Fu,Fv]           = computePDEForcing(X,Y,t,mms,nu);
    [pbx,pby]         = computeGPforP(pbx,pby,u,v,X,Y,t,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fu,Fv);
    [bP]              = setupRhsideP(Nx,Ny,hx,hy,u,v,Fu,Fv,alpha,pbx,pby,mms,X,Y,t,nu);
    [p]               = computeP(LapP,LapPBig,bP,p,Nx,Ny,mms,t,M1,M2);

    Liuold=Liu; Livold=Liv; Leuold=Leu; Levold=Lev; 
    [Leu,Lev]=computeLE(u,v,p,X,Y,t,hx,hy);
    [Liu,Liv]=computeLI(u,v,X,Y,t,hx,hy,nu);

  end
  vorz=(v(3:end,2:end-1)-v(1:end-2,2:end-1))/(2*hx)-(u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2*hy);
  divr=(u(3:end,2:end-1)-u(1:end-2,2:end-1))/(2*hx)+(v(2:end-1,3:end)-v(2:end-1,1:end-2))/(2*hy);
  if (pl == 1 && nt > 5)
      if (mod(nt,10) == 0 )
          subplot(2,2,1)
          contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),vorz,30)
          colorbar
          subplot(2,2,2)
          contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),divr,30)
          colorbar
          subplot(2,2,3)
          contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),u(2:end-1,2:end-1),30)
          colorbar
          subplot(2,2,4)
          contour(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),v(2:end-1,2:end-1),30)
          colorbar
          drawnow
      end
  end
  if (mms == 1)
    e1=u-mmsfun(X,Y,t,1,0,0,0);
    e2=v-mmsfun(X,Y,t,2,0,0,0);
    e1(1,1)=0;e1(1,end)=0;e1(end,1)=0;e1(end,end)=0;
    e2(1,1)=0;e2(1,end)=0;e2(end,1)=0;e2(end,end)=0;

    fprintf('MMS ERROR u, v and p: %e %e %e \n',max(max(abs(e1))),max(max(abs(e2))),...
            max(max(abs(p(2:end-1,2:end-1)-mmsfun(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),t,3,0,0,0)))))
    mesh(e1)
  end
  fprintf('Time step %i took %g seconds\n',nt,cputime-tt)

end

function [u,v]=taylor(X,Y,Lx,Ly,r0,gamma);
R=sqrt((X-Lx).^2+(Y-Ly).^2);
R=gamma*exp(-(R/r0).^2);
u=-(Y-Ly).*R;
v=(X-Lx).*R;

function [p] = computeP(LapP,LapPBig,bP,p,Nx,Ny,mms,t,M1,M2)
%r=sparse(ones((Nx+1)*(Ny+1),1));
%pInner = [[LapP r];[r' 0]]\[bP;0];
TOL = 1e-2;
MAXIT = 500;
bPp = [bP;0];
pInner = LapPBig\bPp;

%p0 = [reshape(p(2:end-1,2:end-1),(Nx+1)*(Ny+1),1);0];
%pInner = qmr(LapPBig,bPp,TOL,MAXIT,M1,M2,p0);
%pInner = gmres(LapPBig,bPp,40,TOL,MAXIT,M1,M2,p0);

%pInner = cgs(LapPBig,bPp,TOL,MAXIT,M1,M2,p0);
%pInner = bicg(LapPBig,bPp,TOL,MAXIT,M1,M2,p0);
beta   = pInner(end);
if (mms == 1)
  pmean=mmsfun(0,0,t,3,0,0,0)-pInner(1);
else
  pmean = 0 ;
end
pInner = pmean+pInner(1:end-1);
[p]    = updateP(p,pInner,Nx,Ny);

function [p] = updateP(p,pInner,Nx,Ny);
nx=(Nx+1);
ny=(Ny+1);
for j=1:ny
  for i=1:nx
    k=i+(j-1)*nx;
    p(i+1,j+1) = pInner(k);
  end
end

function [u,v,uold,vold] = computeUV(Nx,Ny,bU,bV,Leu,Lev,Fu,Fv,...
                                     Leuold,Levold,Liu,Liv,Fuold,Fvold,u,v,uold,vold,k,nu,Lap)


bU=-0.5*bU*nu;
bV=-0.5*bV*nu;
uold=u;
vold=v;
nx=(Nx-1);
ny=(Ny-1);

FFu=u(3:end-2,3:end-2)/k+(1.5*(Leu(2:end-1,2:end-1)+Fu(3:end-2,3:end-2)) ...
                          -0.5*(Leuold(2:end-1,2:end-1)+Fuold(3:end-2,3:end-2))...
                          +0.5*Liu(2:end-1,2:end-1));
FFv=v(3:end-2,3:end-2)/k+(1.5*(Lev(2:end-1,2:end-1)+Fv(3:end-2,3:end-2)) ...
                          -0.5*(Levold(2:end-1,2:end-1)+Fvold(3:end-2,3:end-2))...
                          +0.5*Liv(2:end-1,2:end-1));

for j=1:ny
  for i=1:nx
    k=i+(j-1)*nx;
    bU(k)=bU(k)+FFu(i,j);
    bV(k)=bV(k)+FFv(i,j);
  end
end

U=Lap\bU;
V=Lap\bV;

for j=1:ny
  for i=1:nx
    k=i+(j-1)*nx;
    u(i+2,j+2) = U(k);
    v(i+2,j+2) = V(k);
  end
end

function [LapP] = setupLapP(Nx,Ny,hx,hy);
LapP=sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
hxi2=1/hx^2;
hyi2=1/hy^2;
nx=(Nx+1);
ny=(Ny+1);
for j=1:ny
  for i=1:nx
    k=i+(j-1)*nx;
    kp=k+1;
    km=k-1;
    kb=k-nx;
    kt=k+nx;
    LapP(k,k)=-2*(hxi2+hyi2);
    if (i<nx)
      LapP(k,kp)=hxi2;
      % special stencil at the boundary
      if (i==1)
        LapP(k,kp)=2*hxi2;
      end
    end
    if (i>1)
      LapP(k,km)=hxi2;
      % special stencil at the boundary
      if (i==nx)
        LapP(k,km)=2*hxi2;
      end
    end
    if (j<ny)
      LapP(k,kt)=hyi2;
      % special stencil at the boundary
      if (j==1)
        LapP(k,kt)=2*hyi2;
      end
    end
    if (j>1)
      LapP(k,kb)=hyi2;
      % special stencil at the boundary
      if (j==ny)
        LapP(k,kb)=2*hyi2;
      end
    end
  end
end

function [LapUV] = setupLapUV(Nx,Ny,hx,hy);
LapUV=sparse((Nx-1)*(Ny-1));
hxi2=1/hx^2;
hyi2=1/hy^2;
nx=(Nx-1);
ny=(Ny-1);
for j=1:ny
  for i=1:nx
    k=i+(j-1)*nx;
    kp=k+1;
    km=k-1;
    kb=k-nx;
    kt=k+nx;
    LapUV(k,k)=-2*(hxi2+hyi2);
    if (i<nx)
      LapUV(k,kp)=hxi2;
    end
    if (i>1)
      LapUV(k,km)=hxi2;
    end
    if (j<ny)
      LapUV(k,kt)=hyi2;
    end
    if (j>1)
      LapUV(k,kb)=hyi2;
    end
  end
end

function [bP] = setupRhsideP(Nx,Ny,hx,hy,u,v,Fu,Fv,alpha,pbx,pby,mms,X,Y,t,nu);
bP=zeros((Nx+1)*(Ny+1),1);
hxi2=1/hx^2;
hyi2=1/hy^2;
nx=(Nx+1);
ny=(Ny+1);
D0xu = (u(3:end,2:end-1)-u(1:end-2,2:end-1))/(2*hx);
D0xv = (v(3:end,2:end-1)-v(1:end-2,2:end-1))/(2*hx);
D0yu = (u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2*hy);
D0yv = (v(2:end-1,3:end)-v(2:end-1,1:end-2))/(2*hy);

if (mms == 1)
  xt=X(2:end-1,2:end-1);
  yt=Y(2:end-1,2:end-1);
  Fux=mmsfun(xt,yt,t,1,1,1,0)...
      +mmsfun(xt,yt,t,1,0,1,0).*mmsfun(xt,yt,t,1,0,1,0)...
      +mmsfun(xt,yt,t,1,0,0,0).*mmsfun(xt,yt,t,1,0,2,0)...
      +mmsfun(xt,yt,t,2,0,1,0).*mmsfun(xt,yt,t,1,0,0,1)...
      +mmsfun(xt,yt,t,2,0,0,0).*mmsfun(xt,yt,t,1,0,1,1)...
      +mmsfun(xt,yt,t,3,0,2,0)...
      -nu*(mmsfun(xt,yt,t,1,0,3,0)+mmsfun(xt,yt,t,1,0,1,2));
  Fvy=mmsfun(xt,yt,t,2,1,0,1)...
      +mmsfun(xt,yt,t,1,0,0,1).*mmsfun(xt,yt,t,2,0,1,0)...
      +mmsfun(xt,yt,t,1,0,0,0).*mmsfun(xt,yt,t,2,0,1,1)...
      +mmsfun(xt,yt,t,2,0,0,1).*mmsfun(xt,yt,t,2,0,0,1)...
      +mmsfun(xt,yt,t,2,0,0,0).*mmsfun(xt,yt,t,2,0,0,2)...
      +mmsfun(xt,yt,t,3,0,0,2)...
      -nu*(mmsfun(xt,yt,t,2,0,2,1)+mmsfun(xt,yt,t,2,0,0,3));
else
  Fux  = (Fu(3:end,2:end-1)-Fu(1:end-2,2:end-1))/(2*hx);
  Fvy  = (Fv(2:end-1,3:end)-Fv(2:end-1,1:end-2))/(2*hy);
end
F=alpha*(D0xu+D0yv)-(D0xu.^2+D0yv.^2+2*D0xv.*D0yu)+Fux+Fvy;

% inner points
for j=2:ny-1
  for i=2:nx-1
    k=i+(j-1)*nx;
    bP(k)=F(i,j);
  end
end

% side
for j=1:ny-1:ny
  for i=2:nx-1
    k=i+(j-1)*nx;
    bP(k)=F(i,j);
    if (j==1)
      bP(k)=bP(k)-pby(1,i)*hyi2;
    end
    if (j==ny)
      bP(k)=bP(k)-pby(2,i)*hyi2;
    end
  end
end

% side
for j=2:ny-1
  for i=1:nx-1:nx
    k=i+(j-1)*nx;
    bP(k)=F(i,j);
    if (i==1)
      bP(k)=bP(k)-pbx(1,j)*hxi2;
    end
    if (i==nx)
      bP(k)=bP(k)-pbx(2,j)*hxi2;
    end
  end
end
% corner
for j=1:ny-1:ny
  for i=1:nx-1:nx
    k=i+(j-1)*nx;
    bP(k)=F(i,j);
    if (i==1)
      bP(k)=bP(k)-pbx(1,j)*hxi2;
    end
    if (i==nx)
      bP(k)=bP(k)-pbx(2,j)*hxi2;
    end
    if (j==1)
      bP(k)=bP(k)-pby(1,i)*hyi2;
    end
    if (j==ny)
      bP(k)=bP(k)-pby(2,i)*hyi2;
    end
  end
end

function [bU,bV] = setupRhsideUV(Nx,Ny,hx,hy,gux,guy,gvx,gvy);
bU=zeros((Nx-1)*(Ny-1),1);
bV=zeros((Nx-1)*(Ny-1),1);
hxi2=1/hx^2;
hyi2=1/hy^2;
nx=(Nx-1);
ny=(Ny-1);

% loop over side
for j=2:ny-1
  for i=1:nx-1:nx
    k=i+(j-1)*nx;
    if (i==1)
      bU(k)=bU(k)-gux(1,j+1)*hxi2;
      bV(k)=bV(k)-gvx(1,j+1)*hxi2;
    end
    if (i==nx)
      bU(k)=bU(k)-gux(2,j+1)*hxi2;
      bV(k)=bV(k)-gvx(2,j+1)*hxi2;
    end
  end
end
% loop over side
for j=1:ny-1:ny
  for i=2:nx-1
    k=i+(j-1)*nx;
    if (j==1)
      bU(k)=bU(k)-guy(1,i+1)*hyi2;
      bV(k)=bV(k)-gvy(1,i+1)*hyi2;
    end
    if (j==ny)
      bU(k)=bU(k)-guy(2,i+1)*hyi2;
      bV(k)=bV(k)-gvy(2,i+1)*hyi2;
    end
  end
end
% loop over corners
for j=1:ny-1:ny
  for i=1:nx-1:nx
    k=i+(j-1)*nx;
    if (i==1)
      bU(k)=bU(k)-gux(1,j+1)*hxi2;
      bV(k)=bV(k)-gvx(1,j+1)*hxi2;
    end
    if (i==nx)
      bU(k)=bU(k)-gux(2,j+1)*hxi2;
      bV(k)=bV(k)-gvx(2,j+1)*hxi2;
    end
    if (j==1)
      bU(k)=bU(k)-guy(1,i+1)*hyi2;
      bV(k)=bV(k)-gvy(1,i+1)*hyi2;
    end
    if (j==ny)
      bU(k)=bU(k)-guy(2,i+1)*hyi2;
      bV(k)=bV(k)-gvy(2,i+1)*hyi2;
    end
  end
end

function [Leu,Lev]=computeLE(u,v,p,x,y,t,hx,hy)
% this routine computes the explicit part of the right hand side

Leu=-(+u(2:end-1,2:end-1).*(u(3:end,2:end-1)-u(1:end-2,2:end-1))/(2*hx)...
      +v(2:end-1,2:end-1).*(u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2*hy)...
      +(p(3:end,2:end-1)-p(1:end-2,2:end-1))/(2*hx));
Lev=-(+u(2:end-1,2:end-1).*(v(3:end,2:end-1)-v(1:end-2,2:end-1))/(2*hx)...
      +v(2:end-1,2:end-1).*(v(2:end-1,3:end)-v(2:end-1,1:end-2))/(2*hy)...
      +(p(2:end-1,3:end)-p(2:end-1,1:end-2))/(2*hy));

function [Liu,Liv]=computeLI(u,v,x,y,t,hx,hy,nu)
% this routine computes the implicit part of the right hand side

Liu=nu*((u(3:end,2:end-1)-2*u(2:end-1,2:end-1)+u(1:end-2,2:end-1))/(hx^2)...
        +(u(2:end-1,3:end)-2*u(2:end-1,2:end-1)+u(2:end-1,1:end-2))/(hy^2));
Liv=nu*((v(3:end,2:end-1)-2*v(2:end-1,2:end-1)+v(1:end-2,2:end-1))/(hx^2)...
        +(v(2:end-1,3:end)-2*v(2:end-1,2:end-1)+v(2:end-1,1:end-2))/(hy^2));


function [pbx,pby] = computeGPforP(pbx,pby,u,v,x,y,t,gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,hx,hy,nu,Fu,Fv);
% We compute "ghost values" for p
% pbx will be used to the right hand side in the Poisson eq.
% This approximation uses the curl-curl condition
pbx(1,:)=-2*hx*(-gtux(1,:)...
                -gux(1,:).*(u(3,2:end-1)-u(1,2:end-1))/(2*hx)...
                -gvx(1,:).*(u(2,3:end)  -u(2,1:end-2))/(2*hy)...
                +nu*(-(v(3,3:end)-v(1,3:end)-v(3,1:end-2)+v(1,1:end-2))/(4*hx*hy)...
                     +(u(2,3:end)-2*u(2,2:end-1)+u(2,1:end-2))/(hy^2))...
                +Fu(2,2:end-1));
% right
pbx(2,:)=+2*hx*(-gtux(2,:)...
                -gux(2,:).*(u(end,2:end-1)-u(end-2,2:end-1))/(2*hx)...
                -gvx(2,:).*(u(end-1,3:end)-u(end-1,1:end-2))/(2*hy)...
                +nu*(-(v(end,3:end)-v(end-2,3:end)-v(end,1:end-2)+v(end-2,1:end-2))/(4*hx*hy)...
                     +(u(end-1,3:end)-2*u(end-1,2:end-1)+u(end-1,1:end-2))/(hy^2))...
                +Fu(end-1,2:end-1));
% bottom
pby(1,:)=-2*hy*(-gtvy(1,:)...
                -gvy(1,:).*((v(2:end-1,3)-v(2:end-1,1))')/(2*hy)...
                -guy(1,:).*((v(3:end,2)  -v(1:end-2,2))')/(2*hx)...
                +nu*(((v(3:end,3)-v(1:end-2,3)-v(3:end,1)+v(1:end-2,1))/(4*hy*hx)...
                      +(v(3:end,2) -2*v(2:end-1,2)+v(1:end-2,2))/(hx^2))')...
                +Fv(2:end-1,2)');
% top
pby(2,:)=+2*hy*(-gtvy(2,:)...
                -gvy(2,:).*((v(2:end-1,end)-v(2:end-1,end-2))')/(2*hy)...
                -guy(2,:).*((v(3:end,end-1)-v(1:end-2,end-1))')/(2*hx)...
                +nu*(((v(3:end,end)-v(3:end,end-2)-v(1:end-2,end)+v(1:end-2,end-2))/(4*hy*hx)...
                      +(v(3:end,end-1)-2*v(2:end-1,end-1)+v(1:end-2,end-1))/(hx^2))')...
                +Fv(2:end-1,end-1)');


function [u,v]=computeAndUpdateGPforU(u,v,x,y,hx,hy);
% This routine assumes the interior point and the physical boundary
% conditions have been updated.
% We extrapolate v to the left and right
vbx(1,:)=3*v(2,2:end-1)-3*v(3,2:end-1)+v(4,2:end-1);
vbx(2,:)=3*v(end-1,2:end-1)-3*v(end-2,2:end-1)+v(end-3,2:end-1);
% We extrapolate u to the bottom and top
uby(1,:)=3*u(2:end-1,2)-3*u(2:end-1,3)+u(2:end-1,4);
uby(2,:)=3*u(2:end-1,end-1)-3*u(2:end-1,end-2)+u(2:end-1,end-3);
% update ghost point values
v(1,2:end-1)   = vbx(1,:);
v(end,2:end-1) = vbx(2,:);
u(2:end-1,1)   = uby(1,:)';
u(2:end-1,end) = uby(2,:)';
% we also need to extrapolate data for the stencil
u(1,2)=3*u(2,2)-3*u(3,2)+u(4,2);
u(end,2)=3*u(end-1,2)-3*u(end-2,2)+u(end-3,2);
u(1,end-1)=3*u(2,end-1)-3*u(3,end-1)+u(4,end-1);
u(end,end-1)=3*u(end-1,end-1)-3*u(end-2,end-1)+u(end-3,end-1);
v(2,1)=3*v(2,2)-3*v(2,3)+v(2,4);
v(2,end)=3*v(2,end-1)-3*v(2,end-2)+v(2,end-3);
v(end-1,1)=3*v(end-1,2)-3*v(end-1,3)+v(end-1,4);
v(end-1,end)=3*v(end-1,end-1)-3*v(end-1,end-2)+v(end-1,end-3);
% To the left and right we get u from the zero divergence condition
% to the left
ubx(1,:)=u(3,2:end-1)+(hx/hy)*(v(2,3:end)-v(2,1:end-2));
% to the right
ubx(2,:)=u(end-2,2:end-1)-(hx/hy)*(v(end-1,3:end)-v(end-1,1:end-2));
% At the bottom and top we get v from the zero divergence condition
% at the bottom
vby(1,:)=v(2:end-1,3)+(hy/hx)*(u(3:end,2)-u(1:end-2,2));
vby(2,:)=v(2:end-1,end-2)-(hy/hx)*(u(3:end,end-1)-u(1:end-2,end-1));
% update ghost point values
u(1,2:end-1)   = ubx(1,:);
u(end,2:end-1) = ubx(2,:);
v(2:end-1,1)   = vby(1,:)';
v(2:end-1,end) = vby(2,:)';
% Finally we extrapolate to the corners
u(1,1)=3*u(2,1)-3*u(3,1)+u(4,1);
u(end,1)=3*u(end-1,1)-3*u(end-2,1)+u(end-3,1);
u(1,end)=3*u(2,end)-3*u(3,end)+u(4,end);
u(end,end)=3*u(end-1,end)-3*u(end-2,end)+u(end-3,end);
v(1,1)=3*v(1,2)-3*v(1,3)+v(1,4);
v(1,end)=3*v(1,end-1)-3*v(1,end-2)+v(1,end-3);
v(end,1)=3*v(end,2)-3*v(end,3)+v(end,4);
v(end,end)=3*v(end,end-1)-3*v(end,end-2)+v(end,end-3);

function [u,v]=updateBCforU(u,v,x,y,t,gux,guy,gvx,gvy,hx,hy);
% we update u and v to the left and right.
u(2,2:end-1)     = gux(1,:);
u(end-1,2:end-1) = gux(2,:);
v(2,2:end-1)     = gvx(1,:);
v(end-1,2:end-1) = gvx(2,:);
% Then we do the top and bottom.
u(2:end-1,2)     = guy(1,:)';
u(2:end-1,end-1) = guy(2,:)';
v(2:end-1,2)     = gvy(1,:)';
v(2:end-1,end-1) = gvy(2,:)';

function [Fu,Fv]=computePDEForcing(X,Y,t,mms,nu);
if (mms == 1)
  Fu=mmsfun(X,Y,t,1,1,0,0)...
     +mmsfun(X,Y,t,1,0,0,0).*mmsfun(X,Y,t,1,0,1,0)...
     +mmsfun(X,Y,t,2,0,0,0).*mmsfun(X,Y,t,1,0,0,1)...
     +mmsfun(X,Y,t,3,0,1,0)...
     -nu*(mmsfun(X,Y,t,1,0,2,0)+mmsfun(X,Y,t,1,0,0,2));
  Fv=mmsfun(X,Y,t,2,1,0,0)...
     +mmsfun(X,Y,t,1,0,0,0).*mmsfun(X,Y,t,2,0,1,0)...
     +mmsfun(X,Y,t,2,0,0,0).*mmsfun(X,Y,t,2,0,0,1)...
     +mmsfun(X,Y,t,3,0,0,1)...
     -nu*(mmsfun(X,Y,t,2,0,2,0)+mmsfun(X,Y,t,2,0,0,2));
else
  Fu=0*(X+Y);
  Fv=0*(X+Y);
end

function [gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy]=computeBCForcing(gux,guy,gvx,gvy,gtux,gtuy,gtvx,gtvy,x,y,t,mms);
if (mms == 1)
  gux(1,:)=mmsfun(x(2),y(2:end-1),t,1,0,0,0);
  gux(2,:)=mmsfun(x(end-1),y(2:end-1),t,1,0,0,0);
  guy(1,:)=mmsfun(x(2:end-1),y(2),t,1,0,0,0);
  guy(2,:)=mmsfun(x(2:end-1),y(end-1),t,1,0,0,0);
  gvx(1,:)=mmsfun(x(2),y(2:end-1),t,2,0,0,0);
  gvx(2,:)=mmsfun(x(end-1),y(2:end-1),t,2,0,0,0);
  gvy(1,:)=mmsfun(x(2:end-1),y(2),t,2,0,0,0);
  gvy(2,:)=mmsfun(x(2:end-1),y(end-1),t,2,0,0,0);
  % time derivaives of boundary conditions used to update p.
  gtux(1,:)=mmsfun(x(2),y(2:end-1),t,    1,1,0,0);
  gtux(2,:)=mmsfun(x(end-1),y(2:end-1),t,1,1,0,0);
  gtuy(1,:)=mmsfun(x(2:end-1),y(2),t,    1,1,0,0);
  gtuy(2,:)=mmsfun(x(2:end-1),y(end-1),t,1,1,0,0);
  gtvx(1,:)=mmsfun(x(2),y(2:end-1),t,    2,1,0,0);
  gtvx(2,:)=mmsfun(x(end-1),y(2:end-1),t,2,1,0,0);
  gtvy(1,:)=mmsfun(x(2:end-1),y(2),t, 2,1,0,0);
  gtvy(2,:)=mmsfun(x(2:end-1),y(end-1),t,2,1,0,0);
elseif (mms == 3 || mms == 4)
  gux(1,:)=0;
  gux(2,:)=0;
  guy(1,:)=0;
  guy(2,:)=0;
  gvx(1,:)=1*exp(-((y(2:end-1)-2.5)/0.2).^2);
  gvx(2,:)=-1*exp(-((y(2:end-1)-2.5)/0.2).^2);
  gvy(1,:)=0;
  gvy(2,:)=0;
    % time derivaives of boundary conditions used to update p.
    gtux(1,:)=0;
    gtux(2,:)=0;
    gtuy(1,:)=0;
    gtuy(2,:)=0;
    gtvx(1,:)=0;
    gtvx(2,:)=0;
    gtvy(1,:)=0;
    gtvy(2,:)=0;
elseif (mms == 5)
  gux(1,:)=+(y(2:end-1)-0.5).*exp(-((y(2:end-1)-0.5)/0.2).^2);
  gux(2,:)=+(y(2:end-1)-0.5).*exp(-((y(2:end-1)-0.5)/0.2).^2);
  guy(1,:)=0;
  guy(2,:)=0;
  gvx(1,:)=0;
  gvx(2,:)=0;
  gvy(1,:)=0;
  gvy(2,:)=0;
    % time derivaives of boundary conditions used to update p.
    gtux(1,:)=0;
    gtux(2,:)=0;
    gtuy(1,:)=0;
    gtuy(2,:)=0;
    gtvx(1,:)=0;
    gtvx(2,:)=0;
    gtvy(1,:)=0;
    gtvy(2,:)=0;
else
  gux(1,:)=0;
  gux(2,:)=0;
  guy(1,:)=0;
  guy(2,:)=0;
  gvx(1,:)=0;
  gvx(2,:)=0;
  gvy(1,:)=0;
  gvy(2,:)=0;
  % time derivaives of boundary conditions used to update p.
  gtux(1,:)=0;
  gtux(2,:)=0;
  gtuy(1,:)=0;
  gtuy(2,:)=0;
  gtvx(1,:)=0;
  gtvx(2,:)=0;
  gtvy(1,:)=0;
  gtvy(2,:)=0;
end

function f=mmsfun(x,y,t,field,td,xd,yd);
% Manufactured solution
if (td == 0)
  T=1+t/2;%+(t.^2)/3;
elseif(td == 1)
  T=1/2;%+2*t/3;
elseif(td == 2)
  T=0;%2/3;
else
  T=0;
end

if (xd == 0  && yd == 0)
  if (field == 1)
    f=(x.^2+2*x.*y+y.^2).*T;
  elseif (field == 2)
    f=(x.^2-2*x.*y-y.^2).*T;
  elseif (field == 3)
    f=(x.^2+0.5*x.*y+y.^2-1).*T;
  end
end

if (xd == 1  && yd == 0)
  if (field == 1)
    f=(2*x+2*y).*T;
  elseif (field == 2)
    f=(2*x-2*y).*T;
  elseif (field == 3)
    f=(2*x+0.5*y).*T;
  end
end

if (xd == 0  && yd == 1)
  if (field == 1)
    f=(2*x+2*y).*T;
  elseif (field == 2)
    f=(-2*x-2*y).*T;
  elseif (field == 3)
    f=(0.5*x+2*y).*T;
  end
end


if (xd == 1  && yd == 1)
  if (field == 1)
    f=(2).*T;
  elseif (field == 2)
    f=(-2).*T;
  elseif (field == 3)
    f=(0.5).*T;
  end
end

if (xd == 2  && yd == 0)
  if (field == 1)
    f=(2).*T;
  elseif (field == 2)
    f=(2).*T;
  elseif (field == 3)
    f=(2).*T;
  end
end

if (xd == 0  && yd == 2)
  if (field == 1)
    f=(2).*T;
  elseif (field == 2)
    f=(-2).*T;
  elseif (field == 3)
    f=(2).*T;
  end
end

if (xd+yd > 2)
  if (field == 1)
    f=(0*(x+y)).*T;
  elseif (field == 2)
    f=(0*(x+y)).*T;
  elseif (field == 3)
    f=(0*(x+y)).*T;
  end
end

