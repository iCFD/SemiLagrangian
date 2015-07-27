%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Solving 2-D scalar wave equation with semi-lagrangian scheme
%
%            dq/dt + df/dx + dg/dy = 0,  for x,y \in [a,b;c,d]
%            where f = u*q  and  g = v*q, solving as Dq/Dt = 0
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
   CFL = 8.5;   % CFL condition
  tEnd = 10.0;  % Final time
method = 2;     % {1} Upwind, {2} Semi-Lagrangian
 vfunc = 2;     % see below ;)

% Set analytic velocity Functions
switch vfunc
    case 1  % Linear Periodic advection
        U = @(X,Y,t) -1.0;  % scalar velocity in x direction
        V = @(X,Y,t)  0.5;  % scalar velocity in y direction
    case 2; % Rigid body (washing machine problem ;D)
        U = @(X,Y,t)  2*pi/tEnd*(Y-0.5);
        V = @(X,Y,t) -2*pi/tEnd*(X-0.5);
    case 3; % Shear Cell
        U = @(X,Y,t) -sin(pi*X).*cos(pi*Y)*cos((2*pi/tEnd)*t);
        V = @(X,Y,t)  cos(pi*X).*sin(pi*Y)*cos((2*pi/tEnd)*t);
    case 4; % full shear Cell
        U = @(X,Y,t) -sin(pi*X).*cos(pi*Y);
        V = @(X,Y,t)  cos(pi*X).*sin(pi*Y);
    otherwise
        error('No such %s model in the list',vfunc);
end

%% Preprocess
% Domain discretization
a=0; b=1; c=0; d=1; nx=200; [x,dx,y,dy]=grid2d(a,b,nx,c,d,nx);

% Set velocity Functions
u=U(x,y,0); v=V(x,y,0);

% set IC
q0=IC2d(x,y,2); %{1} 4 Quadrants, {2} Square Jump, {3} Guassian

% Time discretization
dt0=CFL*min(dy,dx);%/max([abs(v(:));abs(u(:))]); % <--check this out!

% set plot range
plotrange=[a,b/dx,c,d/dy,min(min(q0)),max(max(q0))];

%% Solver Loop 
% load initial conditions
qo=q0; it=0; t=0; dt=dt0;

while t<tEnd
    
    % Set time and time step
    if t+dt>tEnd; dt=tEnd-t; end; t=t+dt;
    
    switch method
        case 1	% Forward Euler Upwind
            q = qo-dt*residual2d(qo,U(x,y,t),V(x,y,t),dx,dy); % cfl ~ 0.5
            
        case 2	% Semi-Lagrangian 
            % x is now assumed to be the solutions points in the next time
            % step, so we compute the previous solutions point at the
            % last time step by doing:
            t_half = t-dt/2;
            xm=x-U(x ,y ,t_half)*dt/2;	ym=y-V(x ,y ,t_half)*dt/2;   % step 1
            xs=x-U(xm,ym,t_half)*dt/2;	ys=y-V(xm,ym,t_half)*dt/2;   % step 2
            xp=x-U(xs,ys,t_half)*dt;	yp=y-V(xs,ys,t_half)*dt;     % final 
            
            % Interpolate the solution into the previous solutions points
            q = interp2(x,y,qo,xp,yp,'cubic',0);
            
        otherwise
            error('method not available');
    end
    
    % Update info
    qo=q; it=it+1;
    
    % plot
    if(rem(it,10)==0)
        figure(1); imagesc(x(1,:),y(:,1),q); axis image; axis xy;
        set(gca,'FontWeight','bold','FontSize',14); 
        title(['time: ',num2str(t,'%1.1g')]); drawnow
    end
end

%% Post Processsing

% Compute Error as L1-norm
err=q-q0; % q0 is also our exact solution!

% Calculate Error Norms
abs_err = (dx)*norm(err(:));
rel_err = norm(err(:))/norm(q0(:));

% Plot Final Result 
figure(1); imagesc(x(1,:),y(:,1),q); axis image; axis xy; 
set(gca,'FontWeight','bold','FontSize',14); xlabel('x points'); ylabel('y points');
title(['dx = ',num2str(dx,'%1.1e'),', dy = ',num2str(dy,'%1.1e'),', dt = ',num2str(dt0,'%1.1e'),', time: ',num2str(t,'%1.1g')]);

% Plot Error Map
figure(2); imagesc(x(1,:),y(:,1),err); axis image; axis xy;
set(gca,'FontWeight','bold','FontSize',14);
xlabel('x'); ylabel('y'); title(sprintf('Advection Error, CFL=%g, N=[%d]x[%d]',CFL,nx,nx))
colorbar; hold on; contour(x,y,q,'k'); contour(x,y,q0,'w--');
caxis([ min(min(err)) max(max(err)) ]);
