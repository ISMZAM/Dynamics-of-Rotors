function batch_rotode
%
% Rotating blade model for Dynamics of Rotors / UoB
%
% Properties:
% + Euler-Bernoulli rotating flapping beam
% + constant mass and stiffness distribution
% + root BC: hinge with torsional spring or cantilevered root
% + tip BC: free tip with point or no mass 
%

% Brano Titurus, 25022013
% Last revision: 20112013


% -------------------------------------------------------------------------
% PROBLEM PARAMETERS
% -------------------------------------------------------------------------
% ... blade and rotor parameters
%%ROTOR
width=25.5e-3;
h=12.7e-3;
t=1.4e-3;
L=1.003;
%%Steel thing
Offseth=4.9e-3; %thickness of stee;
Offsetw=22.2e-3; %width of steel
OffsetL=24.5e-3; %length of steel
G=78e9; %G of steel
OffsetE=210e9;
%Calc
Volume=(width*h-(width-2*t)*(h-2*t))*L;
density=2710; %kg/m3
Mass=Volume*density;

m=Mass/L ;% [kg/m], mass density
E=70e9; %Aluminium
J=(1/12)*OffsetL*Offseth*(OffsetL^2+Offseth^2); %J of steel (side)

I=2*(t*h^3)/(12)+2*(((width-2*t)*t^3)/12+(t*(width-2*t)*(0.5*(h-t))^2));
OffsetI=(Offsetw*Offseth^3)/12;

EIf=E*I;              % [N*m^2], flapping stiffness
ef=OffsetL/L;                % [-], non-dimensional flapping hinge offset
R=L+OffsetL;                  % [m], rotor radius
OmegaN=33;              % [rad/s], rotor speed
% ... mode parameters
im=1;                   % [-], mode number
% ... BC parameters
is_spring=true;  % blade root BC status (cantilevered vs. hinge)

 %kf=G*J/Offsetw  ; % [N*m/rad], flapping hinge spring
kf=(OffsetE*OffsetI)/OffsetL;
is_tipmass=false;         % blade tip BC status (free end vs. tip mass)
mt=30e-3;                   % [kg], blade tip mass
% ... auxiliary parameters
Nomegas=20;             % [-], numer of omegas for frequency diagram
% -------------------------------------------------------------------------


% secondary parameters
efR=ef*R;                            % [m], flapping hinge offset
Mbld=m*(R-efR)+mt;                   % [kg], blade mass
Omvec=linspace(0,OmegaN,Nomegas);    % range of rotor speeds
if is_tipmass,mtT=mt;else mtT=0;end  % T(x) change due to tip mass


% -----------------
% Initial estimates
% -----------------
T=4*(R-efR)/(2*im-1);   % [m], wave period
is_init_shape=true;
% fine grid for the first estimated shape
nfine=1e3; xfine=linspace(0,R-efR,nfine);
% calculate estimates
[w_est,W_est]=fcn_init_est(xfine);
% first estimates
% ... eigenvalue
lambda_guess=(w_est)^2;    % [rad/s]^2
% ... mode function and its derivatives
D0wj=W_est;                xfine0=xfine+efR; dx=xfine(2)-xfine(1);
D1wj=diff(D0wj,1)./(dx^1); xfine1=linspace(efR,R,nfine-1);
D2wj=diff(D0wj,2)./(dx^2); xfine2=linspace(efR,R,nfine-2);
D3wj=diff(D0wj,3)./(dx^3); xfine3=linspace(efR,R,nfine-3);


% ----------------------
% Boundary Value Problem
% ----------------------
% bvp4c inputs
io=1;                      % [-], loop counter
Nest=35;                   % [-], number of initial bvp4c grid points
wsvec=zeros(1,Nomegas);    % Omega range

% ... solution pair estimate for bvp4c
solinit = bvpinit(linspace(efR,R,Nest),@fcn_init,lambda_guess);
% ... solver parameters
options = bvpset('stats','off','RelTol',1e-3,'AbsTol',1e-6);%, ...
% ... field equation in bvp4c format
field_eq=@fcn_ode;
% ... boundary conditions in bvp4c format
boundary_cond=@fcn_bc;

for io=1:Nomegas
   fprintf('.')
   % define new rotational speed from the range
   Omega=Omvec(io);
   % BVP solver execution
   sol=bvp4c(field_eq,boundary_cond,solinit,options);
   % solution pair after convergence
   ws=sqrt(sol.parameters); % natural frequency
   xbld=sol.x;              % blade stations
   Wi=sol.y(1,:);           % shape
   dWi=sol.y(2,:);          % diff shape
   ddWi=sol.y(3,:);         % diff diff shape
   dddWi=sol.y(4,:);        % diff diff diff shape
   % update initial estimate
   solinit=bvpinit(xbld,@fcn_init,ws^2);
   % store natural frequency for later use
   wsvec(io)=ws;
end,disp(' ')


% --------------
% Postprocessing
% --------------
% frequency change with Omega
figure, grid on, box on, hold on
plot(Omvec,Omvec/2/pi,'b--')                        % 1R line
text(Omvec(end),Omvec(end)/2/pi,' 1R')              % 1R text
plot(Omvec,wsvec/2/pi,'r.-')                        % frequency line
text(Omvec(end),wsvec(end)/2/pi,sprintf(' F%d',im)) % frequency line text
xlabel('rotor speed [rad/s]'), ylabel('frequency [Hz]')
set(gca,'xlim',[-0.1,1.1*Omvec(end)])

% show last solution pair
figure, grid on, box on, hold on, axis([-efR*1.1,R*1.05,[-1,1]*1.05])
% ... plot hub and hinge
plot([1,1]*efR,get(gca,'ylim'),'--k')
plot([0,0,-efR,efR],[-0.5*efR,0,0,0],'r-','linewidth',3)
plot([efR],[0],'ro','linewidth',2)
% ... plot mode and its second derivative
plot(xbld,ddWi./max(abs(ddWi)),'linewidth',2,'color',[1,1,1]*0.8)
plot(xfine0,D0wj/D0wj(end),':')
plot(xbld,Wi,'linewidth',2)
% ... prepare text for legend, title, etc.
if is_spring,plot(efR,0,'gx','linewidth',2),end
if is_tipmass,plot(R,Wi(end),'ko','markersize',10,'markerfacecolor','m'),end
title(sprintf('Mode %d @ %.1frad/s: InitEst=%.2fHz (%.2fR), FinSol=%.2fHz (%.2fR)', ...
   im,OmegaN,w_est/2/pi,w_est/OmegaN,ws/2/pi,ws/OmegaN))
xlabel('blade span [m]'), ylabel('')
rottxt={'blade root line','rotor hub','Flapping hinge','bending moments', ...
   'initial mode estimate','mode shape @ \Omega_{max}'};
if is_spring, rottxt=[rottxt,{sprintf('k_F=%.1f [N*m/rad]',kf)}];end
if is_tipmass, rottxt=[rottxt,{sprintf('m_T=%.3f [kg]',mt)}];end
hll=legend(rottxt{:},'location','nw');
set(hll,'fontsize',8)

% end-of-batch



%=======================================================
% Local nested functions:
%  fcn_ode
%  fcn_bc
%  fcn_init
%  fcn_init_est
%=======================================================

   % ----------------------------------------------------
   function dydx=fcn_ode(x,y,lambda)
      % FIELD EQUATION: rotating flapping blade
      z=y(1);                      % shape
      u=y(2);                      % diff shape
      v=y(3);                      % diff diff shape
      w=y(4);                      % diff diff diff shape
      EI=(168/L^2)*x^2-(336/L)*x+200;
      EId=(336/L^2)*x-(336/L);
      EIdd=(336/L^2);
      Tx=0.5*m*(L^2-x^2);
      Txd=m*x;
      % model for bvp4c
      dydx = [ ...
         u ; ...
         v ; ...
         w ; ...                   % Euler-Bernoulli flapping model
         (Omega^2*Tx*v+Omega^2*Txd*u+lambda*m*z-2*EId*w-EIdd*v)/(EI)];
         %(lambda*m*z+Omega^2*((0.5*m*(R^2-x^2)+mtT*R)*v-m*x*u))/EIf];
   end
   % ----------------------------------------------------
   function res=fcn_bc(ya,yb,lambda)
      % BOUNDARY CONDITIONS: spring root, mass tip
      
      % left       right   end      deformation
      za=ya(1);   zb=yb(1);        % displacement
      ua=ya(2);   ub=yb(2);        % D_displacement
      va=ya(3);   vb=yb(3);        % DD_displacement
      wa=ya(4);   wb=yb(4);        % DDD_displacement
      % variable BC root
      if is_spring, cond1=va-(kf/200)*ua; % moment continuity
      else          cond1=ua;      % zero slope
      end
      % variable BC tip
      if is_tipmass, cond2=wb+(mtT/EIf)*(lambda*zb-Omega^2*R*ub); % shear continuity
      else           cond2=wb;            % zero shear force
      end
      % boundary conditions and norm for bvp4c
      res = [ ...
         za ; ...                  % blade root BC 1
         cond1 ; ...               % blade root BC 2
         vb ; ...                  % blade tip BC 1
         cond2 ; ...               % blade tip BC 2
         ...
         zb - 1 ];                 % normalisation condition
   end
   % ----------------------------------------------------
   function v=fcn_init(x)
      % INITIAL SOLUTION GUESS
      if io==1  % first estimate
         if is_init_shape
            v = [ ...  % external estimate
               interp1(xfine0,D0wj,x) ; ...
               interp1(xfine1,D1wj,x) ; ...
               interp1(xfine2,D2wj,x) ; ...
               interp1(xfine3,D3wj,x) ];
         else
            v = [ ...  % very crude estimate harmonic functions
                sin(2*pi*x/T) ; ...
                cos(2*pi*x/T)*(2*pi/T) ; ...
               -sin(2*pi*x/T)*(2*pi/T).^2 ; ...
               -cos(2*pi*x/T)*(2*pi/T).^3 ];
         end
      else      % estimates for loop calculations
         v = [ ...
            interp1(xbld,Wi,x) ; ...
            interp1(xbld,dWi,x) ; ...
            interp1(xbld,ddWi,x) ; ...
            interp1(xbld,dddWi,x) ];
      end
   end
   % ----------------------------------------------------
   
   function [ww,Wj]=fcn_init_est(xfine)
      % FIRST INITIAL SOLUTION ESTIMATE
      % + non-rotating cantilevered and pinned beam
      n=4:1e2;
      
      % Fixed and hinged models, free tip
      % ... fixed BC (left)
      blF=[1.88,4.69,7.85,(n-1/2)*pi];
      knF=(cosh(blF)+cos(blF))./(sinh(blF)+sin(blF));
      omegaF=@(m,EI,L,ii)((blF(ii)./L).^2*sqrt(EI./m));
      WnF=@(x,L,ii)(cosh((blF(ii)/L).*x)-cos((blF(ii)/L).*x)-...
         knF(ii).*(sinh((blF(ii)/L).*x)-sin((blF(ii)/L).*x)));
      % ... hinged BC (right - to change to left in future!)
      blH=[3.927,7.069,10.21,(4*n+1)*pi/4];
      knH=(cosh(blH)-cos(blH))./(sinh(blH)-sin(blH));
      if im==1      % RBM hinged BC
         WnH=@(x,L,ii)(1-x./L);
         omegaH=@(m,EI,L,ii)(sqrt(kf./(m*(R-efR)^3/3)));
      else          % EBM hinged BC
         WnH=@(x,L,ii)(cosh((blH(ii-1)/L).*x)+cos((blH(ii-1)/L).*x)-...
            knH(ii-1).*(sinh((blH(ii-1)/L).*x)+sin((blH(ii-1)/L).*x)));
         omegaH=@(m,EI,L,ii)((blH(ii-1)./L).^2*sqrt(EI./m));
      end
      
      % Fixed and hinged calculations
      % ... fixed BC
      wgF=omegaF(m,EIf,R-efR,im);   % [rad/s]
      WnFN=WnF(xfine,R-efR,im)./WnF(xfine(end),R-efR,im);
      % ... hinged BC
      wgH=omegaH(m,EIf,R-efR,im);   % [rad/s]
      WnHN=WnH(xfine,R-efR,im)./WnH(xfine(1),R-efR,im);
      WnHN=WnHN(end:-1:1);
      
      % Fixed-hinged BC blending
      % ... blending parameter
      if is_spring, weigFH=min([kf*(R-efR)*0.05/EIf,1]);
      else          weigFH=1;
      end
      % ... blended modal properties
      ww = weigFH*wgF  + (1-weigFH)*wgH;
      Wj = weigFH*WnFN + (1-weigFH)*WnHN;
   end
end
