function dF = residual1d(q,ax,dx)
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
dq = Upwindflux(ax,q,[0,1])/dx;

% The Residual
dF = dq;
