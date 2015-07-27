function dF = residual2d(q,ax,ay,dx,dy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 2d wave equation using WENO5
%
%                    Residual: RHS of dq/dt
%
%              coded by Manuel Diaz, NTU, 2012.12.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax,zy: scalar advection velocities in x and y directions respectively.

% Along x
dq=   Upwindflux(ax,q,[0,1])/dx;

% Along y
dq=dq+Upwindflux(ay,q,[1,0])/dy;

% The Residual
dF = dq;
