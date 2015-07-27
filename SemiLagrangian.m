%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Solving 1-D wave equation with Semi-Lagrangian advection scheme
%
%                 dq/dt + df/dx = 0,  for x \in [a,b]
%          where f = u*q :: linear flux, solving as Dq/Dt = 0
%
%              coded by Manuel Diaz, NTU, 2015.05.29
%
%       +--------+--------+--------+--------+--------o--------+ n+1 
%       |        |        |        |        |  o     |        |
%       |        |        |        |      o |        |        |
%       +--------+--------+--------+++o++++++--------+--------+ n 
%       |        |        |      o |        |        |        | 
%       |        |        | o      |        |        |        |
%       +--------++++++o+++--------+--------+--------+--------+ n-1 
%       |        |        |        |        |        |        | 
%       m-5     m-4       m-3      m-2      m-1      m        m+1
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
     u = -2.5;	% scalar velocity in x direction
   CFL = 13.51;  % CFL/courant condition (should be large for SL!)
  tEnd = 20.0;  % Final time
method = 2;     % Advection method: {1} Upwind or {2} Semi-Lagrangian;

%% NOTE:
% For the Semi-Lagrangian, if the CFL = 1 or any Natural number, the system
% corresponds to the exact shifting of information over the domain.
% Renderging the interpolation unnecesary. Otherwise, we require the CFL 
% to be real value such that CFL>1.

%% Preprocess

% Domain discretization
a=-1; b=1; Lx=(b-a); nx=200; dx=Lx/nx; x=a+dx/2:dx:b; 

% set IC
condition=1;
switch condition
    case 1
        q0=TestingIC(x);
    case 2
        q0=IC1d(x,9);
end

% Exact solution
qe=q0;

% set plot range
plotrange=[a,b,min(q0)-0.1,max(q0)+0.1];
    
%% Solver Loop 

% Time discretization
dt0=CFL*dx/abs(u);

% load initial conditions
qo=q0; dt=dt0; it=0; t=0;

while t < tEnd
    
    % evalute time step
    if t+dt>tEnd, dt=tEnd-t; end

    % Update time and iteration counter
    it=it+1; t=t+dt;
    
    switch method
        case 1	% Forward Euler Upwind
            q = qo-dt*residual1d(qo,u,dx); % cfl ~ 0.9
        case 2	% SemiLagrangian
            % x: are now the solutions points at new time step, t=t+dt.
            % The old solutions points at the previous time step are
            xo = x-u*dt; % or xo = x-sign(u)*CFL*dx;
            
            % Set periodic BCs
            i=find(xo<a); np=round((b-xo(i))/Lx); xo(i)=xo(i)+np*Lx;
            i=find(xo>b); np=round((xo(i)-a)/Lx); xo(i)=xo(i)-np*Lx;
            
            % Interpolate q data to the old solution points 
            q=interp1(x,qo,xo,'pchip');
  
        otherwise
            error('method not available');
    end
    qo=q;

    % plot
    if(rem(it,10)==0)
        plot(x,q0,x,q,'.'); colormap Copper; axis(plotrange);
        title(['Upwind, dx = ',num2str(dx),', dt = ',num2str(dt0,'%1.1e'),' time: ',num2str(t)])
        xlabel('x'); ylabel('q(x)'); drawnow
    end
end

% Post Process
plot(x,q0,x,q,'.'); colormap Copper; axis(plotrange);
title(['Upwind, dx = ',num2str(dx),', dt = ',num2str(dt0,'%1.1e'),' time: ',num2str(t)])
xlabel('x'); ylabel('q(x)');

% How good is the approximation?
er1 = norm(q-qe,Inf);               disp(['Inf norm: ',num2str(er1)])
er2 = (sum(abs(q-qe).^2)/nx)^0.5;   disp(['L2 norm:  ',num2str(er2)])
er3 = norm(q-qe, 1 );               disp(['L1 norm:  ',num2str(er3)])